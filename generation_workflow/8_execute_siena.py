import json
from pathlib import Path
import logging
import gzip
import os
import subprocess
import tempfile
import argparse


class SienaCaller:
    """Wrapper around the SIENA tool"""

    def __init__(
        self,
        pdb,
        ligand,
        mirror_path,
        ligands_path,
        output_dir,
        tool_path,
        db_path,
        white_list,
    ):
        self.pdb = pdb
        self.ligand = ligand
        self.mirror_path = mirror_path
        self.siena = tool_path
        self.siena_db = db_path
        self.sdf_path = ligands_path / f"{ligand}-{pdb}.sdf"
        self.temp_dir = None
        self.pdb_path = None
        self.final_dir = output_dir
        self.intermediate_dir = output_dir / "intermediate_results"
        self.white_list = white_list

    def run(self):
        """Run SIENA"""
        try:
            self.__make_temp_dir()
            logging.debug("processing in %s", self.temp_dir)
            self.pdb_path = self.__get_pdb()
            exit_code = self.__call_siena().returncode
            if exit_code == 0:
                self.__postprocess_result()
                return exit_code
            else:
                return exit_code
        except subprocess.CalledProcessError as cmdline_error:
            logging.error(cmdline_error)
            return 1
        finally:
            if self.temp_dir:
                self.temp_dir.cleanup()
            # delete intermediate results
            for content in self.intermediate_dir.iterdir():
                if content.is_dir():
                    for file in content.iterdir():
                        file.unlink()
                    content.rmdir()
                else:
                    content.unlink()
            self.intermediate_dir.rmdir()

    def __make_temp_dir(self):
        temp_path = "/tmp/"
        self.temp_dir = tempfile.TemporaryDirectory(dir=temp_path)

    def __get_pdb(self):
        """Gets the PDB file from the PDB mirror and unzips it"""
        pdb_path = Path(self.temp_dir.name) / f"{self.pdb}.pdb"
        with gzip.open(
            str(self.mirror_path / f"pdb{self.pdb}.ent.gz")
        ) as compressed_pdb_file:
            compressed_pdb = compressed_pdb_file.read()
        with pdb_path.open("w") as pdb_file:
            pdb_file.write(compressed_pdb.decode("utf8"))
        return pdb_path

    def __call_siena(self):

        # Create a temporaray subdir for unprocessed results
        if not self.intermediate_dir.is_dir():
            self.intermediate_dir.mkdir()
        command = [
            str(self.siena),
            "--radius",
            "6.5",
            "--identity",
            "1",
            "--bb_rmsd",
            "1",
            "--fragment_length",
            "10",
            "--fragment_distance",
            "4",
            "--resolution",
            "2.5",
            "--flexibility_senitivity",
            "0.6",
            "--holo_only",
            "--filter_unwanted_ligands",
            "false",
            "--database",
            str(self.siena_db),
            "--ligand",
            str(self.sdf_path),
            "--protein",
            str(self.pdb_path),
            "--output",
            str(self.intermediate_dir),
        ]
        logging.debug("calling %s", " ".join(command))
        return subprocess.run(command, check=True)

    def __postprocess_result(self):
        file_prefix = "-".join([self.ligand, self.pdb])
        align_file = self.final_dir / f"{file_prefix}-alignment.txt"
        stats_file = self.final_dir / f"{file_prefix}-resultStatistic.csv"

        ligands_dir = self.intermediate_dir / "ligands"
        if not ligands_dir.is_dir():
            logging.info(f"Siena results were empty for {file_prefix}")
            return

        # Rename files to keep and store them in final dir
        for file in self.intermediate_dir.iterdir():
            if "alignment" in file.name:
                file.rename(str(align_file))
            if "resultStatistic" in file.name:
                file.rename(str(stats_file))

        # Merge SDFs and tag ligands wrt their pdbs
        # > <pdb_id>
        final_ligands_path = self.final_dir / f"{file_prefix}-ligands.sdf"
        final_ligands_file = final_ligands_path.open("a")
        ligand_count = 0
        for file in ligands_dir.iterdir():
            pdb_id = file.name.split(".")[0].split("_")[0]
            lines = file.open("r").readlines()
            ligand = lines[0].strip()
            # make sure the ligand is in our list of accepted ligands after filtering
            # this list is what here is called "white list"
            if ligand in self.white_list[pdb_id]:
                ligand_count += 1
                new_lines = lines[0 : len(lines) - 1]
                new_lines.append("> <pdb_id>\n")
                new_lines.append(f"{pdb_id}\n\n")
                new_lines.append("$$$$\n")
                final_ligands_file.write("".join(new_lines))
            file.unlink()
        # should be empty after the loop
        ligands_dir.rmdir()

        # We do not keep alignments without ligands fitting our acceptance criteria
        if ligand_count < 1:
            os.remove(align_file)
            os.remove(stats_file)
            os.remove(final_ligands_path)
            logging.info(f"No ligand fitting our filters for {self.temp_dir}")


# ----------------------------------------------------------------------------------------


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(
        description="SIENA search for protein-ligand complexes."
    )
    arg_parser.add_argument(
        "--config",
        type=str,
        required=True,
        help="Path to config in JSON format including all file paths.",
    )
    args = arg_parser.parse_args()

    config = json.load(open(args.config, "r"))

    logging.basicConfig(level=logging.DEBUG)

    if not Path(config["6_finally_filtered_pdbs_ligands_json"]).is_file():
        print("Input file initial_csv specified in config does not exist")
        raise FileNotFoundError(config["6_finally_filtered_pdbs_ligands_json"])

    # files to handle
    pdb_dict = json.load(open(config["6_finally_filtered_pdbs_ligands_json"], "r"))
    tool_path = Path(config["tool_siena"])
    mirror_path = Path(config["PDB_file_mirror"])
    ligands_path = Path(config["4_removed_ligands_dir"]) / "sdfs"
    siena_db = Path(config["7_siena_db"])
    output_dir = Path(config["8_siena_results_dir"])

    # create result dir if necessary
    if not output_dir.is_dir():
        output_dir.mkdir()

    for pdb in pdb_dict.keys():
        for ligand in pdb_dict[pdb]:
            # watchout for the case of the pdbs because of the paths
            siena_call = SienaCaller(
                pdb.lower(),
                ligand,
                mirror_path,
                ligands_path,
                output_dir,
                tool_path,
                siena_db,
                pdb_dict,
            )
            siena_result = siena_call.run()
            logging.info(
                f"SIENA call for ligand {ligand} in PDB {pdb} finished"
                " with exit code {siena_result}"
            )

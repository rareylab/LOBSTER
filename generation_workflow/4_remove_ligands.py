import logging
import gzip
from pathlib import Path
import subprocess
import tempfile
import argparse, json
import pandas as pd


class LigandRemoval:
    """Wrapper around the remove_ligand NAOMI tool"""

    def __init__(self, pdb_path, ligand_name, tool_path, results_path):
        """Removes a specified ligand from a given PDB

        Returns an edf of the pocket, the ligand in an SD file
        and the pocket in binary format. Parameters of remove_ligand
        tool were kept as defaults.
        """
        self.pdb_path = pdb_path
        self.pdb_name = pdb_path.name.replace("pdb", "").replace(".ent.gz", "")
        self.ligand_name = ligand_name
        self.tool_path = tool_path
        self.results_dir = results_path
        self.temp_dir = None
        self.complex_path = None

    def run(self):
        """Run the ligand removal"""
        try:
            self.__make_temp_dir()
            logging.debug("processing in %s", self.temp_dir.name)
            self.complex_path = self.__get_pdb()
            exit_code = self.__remove_ligand().returncode
            return exit_code
        except subprocess.CalledProcessError as cmdline_error:
            logging.error(cmdline_error)
            return 1
        finally:
            if self.temp_dir:
                self.temp_dir.cleanup()

    def __make_temp_dir(self):
        temp_path = "/tmp/"
        self.temp_dir = tempfile.TemporaryDirectory(dir=temp_path)

    def __get_pdb(self):
        """Gets the PDB file from the PDB mirror and unzips it"""
        pdb_path = Path(self.temp_dir.name) / f"{self.pdb_name}.pdb"
        with gzip.open(self.pdb_path) as compressed_pdb_file:
            compressed_pdb = compressed_pdb_file.read()
        with pdb_path.open("w") as pdb_file:
            pdb_file.write(compressed_pdb.decode("utf8"))
        return pdb_path

    def __remove_ligand(self):
        """Calls remove_ligand"""
        command = [
            str(self.tool_path),
            "--protein",
            str(self.complex_path),
            "--ligand",
            str(self.ligand_name),
            "--output",
            f"{self.results_dir}/binaries/{self.ligand_name}-{self.pdb_name}.db",
            "--output_ligand",
            f"{self.results_dir}/sdfs/{self.ligand_name}-{self.pdb_name}.sdf",
            "--output_edf",
            f"{self.results_dir}/edfs/{self.ligand_name}-{self.pdb_name}.edf",
        ]
        logging.debug("calling %s", " ".join(command))
        return subprocess.run(command, check=True)


# ----------------------------------------------------------------------------------------


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(description="Remove ligands.")
    arg_parser.add_argument(
        "--config",
        type=str,
        required=True,
        help="Path to config in JSON format including all file paths.",
    )
    args = arg_parser.parse_args()

    config = json.load(open(args.config, "r"))

    args = arg_parser.parse_args()
    master_table_path = Path(config["3_activities_csv"])
    output_path = Path(config["4_removed_ligands_dir"])
    mirror_path = Path(config["PDB_file_mirror"])
    tool_path = Path(config["tool_remove_ligands"])
    # create result dirs if necessary
    if not output_path.is_dir():
        output_path.mkdir()
    if not (output_path / "edfs").is_dir():
        (output_path / "edfs").mkdir()
    if not (output_path / "sdfs").is_dir():
        (output_path / "sdfs").mkdir()
    if not (output_path / "binaries").is_dir():
        (output_path / "binaries").mkdir()

    logging.basicConfig(level=logging.DEBUG)
    zipped = list()
    df = pd.read_csv(master_table_path, delimiter=",")
    for idx, row in df.iterrows():
        zipped.append((row["pdb"], row["name"]))

    logging.info(f"Removing {len(zipped)} ligands...")
    for pdb, ligand in zipped:
        pdb = pdb.lower()
        pdb_mirror_path = mirror_path / f"pdb{pdb.lower()}.ent.gz"
        ligand_removal = LigandRemoval(pdb_mirror_path, ligand, tool_path, output_path)
        removal_result = ligand_removal.run()
        logging.info(
            f"Ligand removal for {ligand} in PDB {pdb} finished"
            f" with exit code {removal_result}"
        )

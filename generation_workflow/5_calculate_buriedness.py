import logging
import gzip
import json
import subprocess
import tempfile
from pathlib import Path
import argparse


class PropertyCalculation:
    """Wrapper around the complex_properties from the NAOMI library"""

    def __init__(self, tool_path, pdb_path, ligand_path, result_file):
        """Calculates SAS of ligand

        Complex Properties profiles SAS of free ligand and of
        bound ligand inside the complex. Console ouutput gives
        stats.
        """
        self.result_file = result_file
        self.ligand_path = ligand_path
        self.complex_properties_path = tool_path
        self.pdb_path = pdb_path
        self.pdb_name = pdb_path.name.replace("pdb", "").replace(".ent.gz", "")
        self.temp_dir = None
        self.complex_path = None

    def run(self):
        """Calculate the complex properties"""
        try:
            self.__make_temp_dir()
            logging.debug("processing in %s", self.temp_dir.name)
            self.complex_path = self.__get_pdb()
            calculation_result = self.__calculate_properties()
            exit_code = calculation_result.returncode
            processing_code = self.__process_stdout(calculation_result.stdout.decode())
            return exit_code, processing_code
        except subprocess.CalledProcessError as cmdline_error:
            log_file = "complex_properties.log"
            logging.error(cmdline_error)
            return 1, 1
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

    def __calculate_properties(self):
        """Calls the complex properties tool"""
        command = [
            str(self.complex_properties_path),
            "--complex",
            str(self.complex_path),
            "--ligand",
            str(self.ligand_path),
        ]
        logging.debug("calling %s", " ".join(command))
        return subprocess.run(command, check=True, stdout=subprocess.PIPE)

    def __process_stdout(self, outprints):
        buriedness = 0
        processing_code = 1
        for line in outprints.strip().split("\n"):
            if pdb.upper() in line:
                result_elements = line.split(";")
                bound_sasa = float(result_elements[4])
                free_sasa = float(result_elements[5])
                buriedness = 1 - (bound_sasa / free_sasa)
                with self.result_file.open("a") as f:
                    f.write(line)
                    f.write(f";{buriedness}\n")
                processing_code = 0
        return processing_code


# ----------------------------------------------------------------------------------------


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(description="Calculation of buriedness.")
    arg_parser.add_argument(
        "--config",
        type=str,
        required=True,
        help="Path to config in JSON format including all file paths.",
    )
    args = arg_parser.parse_args()

    config = json.load(open(args.config, "r"))

    mirror_path = Path(config["PDB_file_mirror"])
    ligands_path = Path(config["4_removed_ligands_dir"]) / "sdfs"
    tool_path = Path(config["tool_complex_properties"])
    property_results_path = Path(config["5_buriednes_csv"])

    # create result dir if necessary
    if not property_results_path.parent.is_dir():
        property_results_path.parent.mkdir()

    # check ligand, tool and mirror path
    if not mirror_path.is_dir():
        raise FileNotFoundError("PDB file mirror directory does not exist")
    if not ligands_path.is_dir():
        raise FileNotFoundError("Ligand SDF directory does not exist")
    if not tool_path.is_file():
        raise FileNotFoundError("Tool does not exist")

    with property_results_path.open("w") as f:
        f.write(
            "pdb;name;min_dist_protein;min_dist_any;real_ligand_sasa;free_ligand_sasa;"
            "atoms_without_coordinates;buriedness\n"
        )

    logging.basicConfig(level=logging.DEBUG)

    for ligand_path in ligands_path.iterdir():
        ligand, pdb = ligand_path.stem.split("-")
        pdb_mirror_path = mirror_path / f"pdb{pdb.lower()}.ent.gz"
        property_calculation = PropertyCalculation(
            tool_path, pdb_mirror_path, ligand_path, property_results_path
        )
        calculation_result = property_calculation.run()
        logging.info(
            f"Buriedness Calculation for {ligand} in PDB {pdb}"
            " finished with exit code {calculation_result}"
        )

import progressbar
import os
import json
import argparse
from pathlib import Path
import logging, subprocess


if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)

    arg_parser = argparse.ArgumentParser(
        description="Symlink PDBs stored in a json file"
        " and generate a SIENA database with these PDBs"
    )
    arg_parser.add_argument(
        "--config",
        type=str,
        required=True,
        help="Path to config in JSON format including all file paths.",
    )
    args = arg_parser.parse_args()

    config = json.load(open(args.config, "r"))

    if not Path(config["6_finally_filtered_pdbs_ligands_json"]).is_file():
        print("Input file initial_csv specified in config does not exist")
        raise FileNotFoundError(config["6_finally_filtered_pdbs_ligands_json"])

    # create symlinks
    pdb_dict = json.load(open(config["6_finally_filtered_pdbs_ligands_json"], "r"))
    output_dir = Path(config["7_pdb_symlinks_dir"])

    # create pdb dir if necessary
    if not output_dir.is_dir():
        output_dir.mkdir()
    for pdb in progressbar.progressbar(pdb_dict.keys()):
        file_name = f"pdb{pdb.lower()}.ent.gz"
        os.system(
            f"ln -s {Path(config['PDB_file_mirror']) / file_name} {output_dir / file_name}"
        )
    logging.info("Finished creating symlinks, generating SIENA DB")

    # -f 0 is for PDB file type ent.gz
    command = [
        config["tool_generate_siena_database"],
        "-d",
        str(output_dir),
        "-b",
        config["7_siena_db"],
        "-f",
        "0",
    ]
    logging.debug("calling %s", " ".join(command))
    try:
        result = subprocess.run(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
            check=True,
        )
        logging.debug(result.stdout)
        logging.error(result.stderr)
        logging.info(f"Finished generating SIENA DB with exit code {result.returncode}")
    except subprocess.CalledProcessError as cmdline_error:
        logging.error(cmdline_error)

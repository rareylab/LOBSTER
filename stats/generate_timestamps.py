from pathlib import Path
import gzip, tempfile
from progressbar import progressbar
import json
import pandas as pd
from argparse import ArgumentParser
from collections import defaultdict


arg_parser = ArgumentParser(
    description="Create statistics on timestamps of the molecules"
    " within the LOBSTER dataset"
)
arg_parser.add_argument(
    "--config",
    type=str,
    required=True,
    help="Path to config in JSON format including all file paths.",
)
args = arg_parser.parse_args()

config = json.load(open(args.config, "r"))

mirror_path = Path(config["PDB_file_mirror"])
all_benchmark_pdbs = Path(config["stats_all_pdbs_in_LOBSTER"])
if not all_benchmark_pdbs.is_file():
    print("Input directory 10_final_dataset_dir as specified in config does not exist")
    raise FileNotFoundError(str(all_benchmark_pdbs))
final_dataset_dir = Path(config["10_final_dataset_dir"])
if not final_dataset_dir.is_dir():
    print("Input directory 10_final_dataset_dir as specified in config does not exist")
    raise IOError(str(final_dataset_dir))

all_pairs_csv = final_dataset_dir / "subsets" / "all_pairs.csv"

temp_path = "/tmp/"
temp_file = tempfile.TemporaryFile(suffix=".pdb")
pdb_list = [x.strip() for x in all_benchmark_pdbs.open("r").readlines()]
pdb_timestamp_touple = []

for pdb in progressbar(pdb_list):
    # expects the PDB to be named pdbXXXX.ent.gz
    pdb_mirror_path = mirror_path / ("pdb" + pdb.lower() + ".ent.gz")
    with gzip.open(pdb_mirror_path) as compressed_pdb_file:
        compressed_pdb = compressed_pdb_file.read()
        uncompressed_pdb = compressed_pdb.decode("utf8")
        rev_dates = [x for x in uncompressed_pdb.split("\n") if "REVDAT" in x]
        for rev in rev_dates:
            split_rev = rev.split()
            revision = split_rev[1]
            if int(revision) < 1:
                # we expect the lowest revision to be 1
                print("REVISION 0:", pdb)
            if int(revision) == 1:
                timestamp = split_rev[2]
                pdb_timestamp_touple.append((pdb, timestamp))
                break

timestamps = pd.DataFrame.from_records(
    pdb_timestamp_touple, columns=["PDB code", "timestamp"]
)
all_pairs = defaultdict(
    lambda: {
        "count": 0,
        "pdb_resolutions": [],
        "edia_m_values": [],
        "mw_values": [],
        "hac_values": [],
        "rot_bonds": [],
        "pdb_ids": [],
    }
)
molecular_props = pd.read_csv(config["stats_molecules"], delimiter=";")

with all_pairs_csv.open("r") as f:
    for line in progressbar(f.readlines()):
        query, template, _ = line.split(";")
        ligand_1, pdb_1 = query.split(".")[0].split("-")
        ligand_2, pdb_2 = template.split(".")[0].split("-")
        ts_1 = int(
            timestamps.where(timestamps["PDB code"] == pdb_1)
            .dropna(subset=["PDB code"])["timestamp"]
            .values[0]
            .split("-")[2]
        )
        if ts_1 > 24:
            ts_1 = int("19" + str(ts_1))
        elif ts_1 < 10:
            ts_1 = int("200" + str(ts_1))
        else:
            ts_1 = int("20" + str(ts_1))
        ts_2 = int(
            timestamps.where(timestamps["PDB code"] == pdb_2)
            .dropna(subset=["PDB code"])["timestamp"]
            .values[0]
            .split("-")[2]
        )
        if ts_2 > 24:
            ts_2 = int("19" + str(ts_2))
        elif ts_2 < 10:
            ts_2 = int("200" + str(ts_2))
        else:
            ts_2 = int("20" + str(ts_2))

        if ts_1 > ts_2:
            all_pairs[ts_1]["count"] += 1  # this is counted in pairs
        else:
            all_pairs[ts_2]["count"] += 1

# add other properties besides the counts, not per pair but per molecule
for idx, row in progressbar(molecular_props.iterrows()):
    year = int(
        timestamps.where(timestamps["PDB code"] == row["pdb"])
        .dropna(subset=["PDB code"])["timestamp"]
        .values[0]
        .split("-")[2]
    )
    if year > 24:
        year = int("19" + str(year))
    elif year < 10:
        year = int("200" + str(year))
    else:
        year = int("20" + str(year))

    all_pairs[year]["pdb_resolutions"].append(row["pdb_resolution"])
    all_pairs[year]["edia_m_values"].append(row["edia_m"])
    all_pairs[year]["mw_values"].append(row["mw"])
    all_pairs[year]["hac_values"].append(row["hac"])
    all_pairs[year]["rot_bonds"].append(row["rot_bonds"])
    all_pairs[year]["pdb_ids"].append(row["pdb"])


timestamps_json = Path(config["timestamps_json"])
with timestamps_json.open("w") as f:
    json.dump(all_pairs, f)

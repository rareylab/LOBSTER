from pathlib import Path
from argparse import ArgumentParser
import json

import matplotlib.pyplot as plt
import pandas as pd

from matplotlib_venn import venn2
import plotnine as p9


EXTENSIVE = False


def plot_venn_diagram(
    a=0,
    b=0,
    ab=0,
    a_label="LOBSTER",
    b_label="AstraZeneca",
    title="test",
    output_name="try.png",
):
    plt.figure(figsize=(8, 8))
    figure = venn2(
        subsets=(a, b, ab),
        set_labels=[a_label, b_label],
        set_colors=["red", "cyan"],
        alpha=0.5,
    )
    for text in figure.set_labels:
        text.set_fontsize(20)
    for text in figure.subset_labels:
        text.set_fontsize(18)
    plt.title(title, fontsize=20)
    plt.savefig(plot_output / output_name, dpi=300)
    plt.close()


# ----------------------------------------------------------------------------------------


if __name__ == "__main__":
    arg_parser = ArgumentParser(
        description="Create diversity plots for the LOBSTER manuscript"
        " and the list of protein classes"
    )
    arg_parser.add_argument(
        "--config",
        type=str,
        required=True,
        help="Path to config in JSON format including all file paths.",
    )
    args = arg_parser.parse_args()

    config = json.load(open(args.config, "r"))

    az_broad_path = Path(config["az_broad_cluster_csv"])
    plot_output = Path(config["plots_dir"])
    pfam_families_csv = Path(config["pfam_classes_csv"])
    lobster_broad_path = Path(config["broad_cluster_csv"])

    pfam_lookup = {}
    az_protein_classes = {}
    lobster_protein_classes = {}
    pfam_families_data = pd.read_csv(pfam_families_csv, index_col=None)
    for index, row in pfam_families_data.iterrows():
        pfam = row["pfam"]
        protein_classes = row["class"]
        pfam_lookup[pfam] = protein_classes
        print(pfam, protein_classes)

    seperator = ";"
    az_dataset_data = pd.read_csv(az_broad_path, dtype=str)

    not_found_entries = []
    az_uniprot_ids = []
    az_pfam_ids = set()
    az_pfam_names = set()

    for index, row in az_dataset_data.iterrows():
        found = False
        cluster = row["Cluster"]
        original_uniprot = row["Cluster"].split(".")[0].split("_")[0]
        new_uniprots = (
            row["Uniprotids"]
            .replace("{", "")
            .replace("}", "")
            .replace("'", "")
            .replace("'", "")
        )
        PfamIDs = (
            row["PfamID"]
            .replace("{", "")
            .replace("}", "")
            .replace("'", "")
            .replace("'", "")
        )
        PfamNames = (
            row["PfamName"]
            .replace("{", "")
            .replace("}", "")
            .replace("'", "")
            .replace("'", "")
        )
        for uniprot in new_uniprots.split(","):
            uniprot = uniprot.strip()
            if uniprot == original_uniprot:
                az_uniprot_ids.append(uniprot)
                found = True
        for pfamid in PfamIDs.split(seperator):
            if pfamid != "failed":
                az_pfam_ids.add(pfamid)
                protein_class = pfam_lookup[pfamid]
                if cluster in az_protein_classes:
                    az_protein_classes[cluster].add(protein_class)
                else:
                    az_protein_classes[cluster] = set([protein_class])
        for pfamid in PfamNames.split(seperator):
            if pfamid != "failed":
                az_pfam_names.add(pfamid)
        if not found:
            print(original_uniprot, new_uniprots)
            not_found_entries.append(original_uniprot)
            az_uniprot_ids.append(new_uniprots)

    if EXTENSIVE:
        print(
            f"There are {len(not_found_entries)} targets in the az dataset that are not"
            f" reproduced from the original uniprot ids due to discontinuation by uniprot",
            not_found_entries,
        )
        print(
            f"There are {len(az_uniprot_ids)} targets in the az dataset with"
            f" {len(set(az_uniprot_ids))} unique uniprot ids",
            az_uniprot_ids,
        )
        print(
            f"There are {len(az_pfam_ids)} unique pfam ids in the az dataset",
            az_pfam_ids,
        )
        print(
            f"There are {len(az_pfam_names)} unique pfam names in the az dataset",
            az_pfam_names,
        )
    else:
        print(
            f"There are {len(not_found_entries)} targets in the az dataset that are not"
            f" reproduced from the original uniprot ids due to discontinuation by uniprot"
        )
        print(
            f"There are {len(az_uniprot_ids)} targets in the az dataset with"
            f" {len(set(az_uniprot_ids))} unique uniprot ids"
        )
        print(f"There are {len(az_pfam_ids)} unique pfam ids in the az dataset")
        print(f"There are {len(az_pfam_names)} unique pfam names in the az dataset")

    lobster_uniprot_ids = []
    lobster_pfam_ids = set()
    lobster_pfam_names = set()

    lobster_data = pd.read_csv(lobster_broad_path, dtype=str)
    for index, row in lobster_data.iterrows():
        cluster = row["Cluster"]
        new_uniprots = (
            row["Uniprotids"]
            .replace("{", "")
            .replace("}", "")
            .replace("'", "")
            .replace("'", "")
        )
        PfamIDs = (
            row["PfamID"]
            .replace("{", "")
            .replace("}", "")
            .replace("'", "")
            .replace("'", "")
        )
        PfamNames = (
            row["PfamName"]
            .replace("{", "")
            .replace("}", "")
            .replace("'", "")
            .replace("'", "")
        )
        for uniprot in new_uniprots.split(","):
            uniprot = uniprot.strip()
            if uniprot != "failed":
                lobster_uniprot_ids.append(uniprot)
        for pfamid in PfamIDs.split(seperator):
            if pfamid != "failed":
                lobster_pfam_ids.add(pfamid)
                protein_class = pfam_lookup[pfamid]
                if cluster in lobster_protein_classes:
                    lobster_protein_classes[cluster].add(protein_class)
                else:
                    lobster_protein_classes[cluster] = set([protein_class])
        for pfamid in PfamNames.split(seperator):
            if pfamid != "failed":
                lobster_pfam_names.add(pfamid)

    if EXTENSIVE:
        print(
            f"There are {len(lobster_uniprot_ids)} targets in the LOBSTER dataset with"
            f" {len(set(lobster_uniprot_ids))} unique uniprot ids",
            lobster_uniprot_ids,
        )
        print(
            f"There are {len(lobster_pfam_ids)} unique pfam ids in the LOBSTER dataset",
            lobster_pfam_ids,
        )
        print(
            f"There are {len(lobster_pfam_names)} unique pfam ids in the LOBSTER dataset",
            lobster_pfam_names,
        )
    else:
        print(
            f"There are {len(lobster_uniprot_ids)} targets in the LOBSTER dataset with"
            f" {len(set(lobster_uniprot_ids))} unique uniprot ids"
        )
        print(
            f"There are {len(lobster_pfam_ids)} unique pfam ids in the LOBSTER dataset"
        )
        print(
            f"There are {len(lobster_pfam_names)} unique pfam ids in the LOBSTER dataset"
        )

    lobster_uniprot_ids = set(lobster_uniprot_ids)
    az_uniprot_ids = set(az_uniprot_ids)

    print("Uniprot Analysis")
    intersection = lobster_uniprot_ids & az_uniprot_ids
    lobster_more = lobster_uniprot_ids.difference(az_uniprot_ids)
    az_more = az_uniprot_ids.difference(lobster_uniprot_ids)
    if EXTENSIVE:
        print(
            f"There are {len(intersection)} uniprot ids that both datasets have in common",
            intersection,
        )
        print(
            f"There are {len(lobster_more)} uniprot ids that are unique to LOBSTER dataset",
            lobster_more,
        )
        print(
            f"There are {len(az_more)} uniprot ids that are unique to az dataset",
            az_more,
        )
    else:
        print(
            f"There are {len(intersection)} uniprot ids that both datasets have in common"
        )
        print(
            f"There are {len(lobster_more)} uniprot ids that are unique to LOBSTER dataset"
        )
        print(f"There are {len(az_more)} uniprot ids that are unique to az dataset")

    plot_venn_diagram(
        a=len(lobster_more),
        b=len(az_more),
        ab=len(intersection),
        title="UniProt Accession Number Overlap",
        output_name="uniprot_venn.svg",
    )

    print("PfamID Analysis")
    intersection = lobster_pfam_ids & az_pfam_ids
    lobster_more = lobster_pfam_ids.difference(az_pfam_ids)
    az_more = az_pfam_ids.difference(lobster_pfam_ids)
    if EXTENSIVE:
        print(
            f"There are {len(intersection)} Pfam ids that both datasets have in common",
            intersection,
        )
        print(
            f"There are {len(lobster_more)} Pfam ids that are unique to LOBSTER dataset",
            lobster_more,
        )
        print(
            f"There are {len(az_more)} Pfam ids that are unique to az dataset", az_more
        )
    else:
        print(
            f"There are {len(intersection)} Pfam ids that both datasets have in common"
        )
        print(
            f"There are {len(lobster_more)} Pfam ids that are unique to LOBSTER dataset"
        )
        print(f"There are {len(az_more)} Pfam ids that are unique to az dataset")

    plot_venn_diagram(
        a=len(lobster_more),
        b=len(az_more),
        ab=len(intersection),
        title="Pfam ID Overlap",
        output_name="pfamid_venn.svg",
    )
    print(az_more)
    print("Pfam Name Analysis")
    intersection = lobster_pfam_names & az_pfam_names
    lobster_more = lobster_pfam_names.difference(az_pfam_names)
    az_more = az_pfam_names.difference(lobster_pfam_names)
    if EXTENSIVE:
        print(
            f"There are {len(intersection)} Pfam names that both datasets have in common",
            intersection,
        )
        print(
            f"There are {len(lobster_more)} Pfam names that are unique to LOBSTER dataset",
            lobster_more,
        )
        print(
            f"There are {len(az_more)} Pfam names that are unique to az dataset",
            az_more,
        )
    else:
        print(
            f"There are {len(intersection)} Pfam names that both datasets have in common"
        )
        print(
            f"There are {len(lobster_more)} Pfam names that are unique to LOBSTER dataset"
        )
        print(f"There are {len(az_more)} Pfam names that are unique to az dataset")

    plot_venn_diagram(
        a=len(lobster_more),
        b=len(az_more),
        ab=len(intersection),
        title="Pfam Name Overlap",
        output_name="pfamname_venn.svg",
    )
    print(az_more)

    protein_class_dict = {}
    total_az = 0
    total_lobster = 0
    for cluster in az_protein_classes:
        dataset_string = "az"
        total_az += 1
        for protein_class in az_protein_classes[cluster]:
            if protein_class in protein_class_dict:
                if dataset_string in protein_class_dict[protein_class]:
                    protein_class_dict[protein_class][dataset_string] += 1
                else:
                    protein_class_dict[protein_class][dataset_string] = 1
            else:
                protein_class_dict[protein_class] = {dataset_string: 1}

    for cluster in lobster_protein_classes:
        total_lobster += 1
        dataset_string = "lobster"
        for protein_class in lobster_protein_classes[cluster]:
            if protein_class in protein_class_dict:
                if dataset_string in protein_class_dict[protein_class]:
                    protein_class_dict[protein_class][dataset_string] += 1
                else:
                    protein_class_dict[protein_class][dataset_string] = 1
            else:
                protein_class_dict[protein_class] = {dataset_string: 1}

    print(protein_class_dict)
    final_data = []
    for protein_class in protein_class_dict:
        if "az" in protein_class_dict[protein_class]:
            n_az = protein_class_dict[protein_class]["az"]
        else:
            n_az = 0
        if "lobster" in protein_class_dict[protein_class]:
            n_lobster = protein_class_dict[protein_class]["lobster"]
        else:
            n_lobster = 0
        to_append_1 = [protein_class, n_az, n_az * 100 / total_az, "az"]
        to_append_2 = [
            protein_class,
            n_lobster,
            n_lobster * 100 / total_lobster,
            "lobster",
        ]
        final_data.append(to_append_1)
        final_data.append(to_append_2)

    final_data = pd.DataFrame(
        final_data, columns=["Protein Class", "Total", "Percentage", "Dataset"]
    )
    # For propper labels in plots
    final_data.replace("az", "AstraZeneca", inplace=True)
    final_data.replace("lobster", "LOBSTER", inplace=True)
    final_data.sort_values(by="Protein Class")
    final_data.to_csv(config["protein_class_csv"])

    plot = p9.ggplot(final_data)
    plot += p9.geom_col(
        p9.aes(x="Protein Class", y="Total", fill="Dataset"), position="dodge"
    )
    plot += p9.theme(axis_text_x=p9.element_text(rotation=45, hjust=1))
    plot.save(plot_output / "total_protein_class_diagram.png", dpi=300)

    plot = p9.ggplot(final_data)
    plot += p9.geom_col(
        p9.aes(x="Protein Class", y="Percentage", fill="Dataset"), position="dodge"
    )
    plot += p9.theme(axis_text_x=p9.element_text(rotation=45, hjust=1))
    plot.save(plot_output / "percentage_protein_class_diagram.png", dpi=300)

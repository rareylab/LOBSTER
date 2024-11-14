<img src="/assets/lobster_logo.png" width="350" align="right">

# LOBSTER

## Description
LOBSTER ("Ligand Overlays from Binding SiTe Ensemble Representatives") is a dataset of ligand overlays designed to evaluate small molecule superposition tools.
It is available for download at [Zenodo](https://doi.org/10.5281/zenodo.12658320).
This repository provides the scripts to generate and analyze the dataset.

## Workflow: Creating the LOBSTER Dataset
Please provide all paths in **config.json** so that the scripts will work properly.
Numbered path variables belong to the scripts with the respective number in the file name.
This number indicates the order in which the scripts are planned to be executed. This is detailed below.
Important note on the PDB file mirror: All PDB files in the mirror directory should have the format pdb<pdbIdLowerCase>.ent.gz.
Tool executables needed for this workflow are RemoveLigand (tool derived from the [NAOMI ChemBio Suite](https://software.zbh.uni-hamburg.de/)) and [SIENA](https://doi.org/10.1021/acs.jcim.5b00588).
The starting point of the workflow is the merged output of the [StructureProfiler](https://doi.org/10.1093/bioinformatics/bty692) and LigandFinder (tool derived from the [NAOMI ChemBio Suite](https://software.zbh.uni-hamburg.de/)), e.g. provided with the latest LOBSTER version as **ligand_structure_data.csv**. It can be used as **initial_csv** in the config.

All scripts required to execute the following steps can be found in the **generation_workflow** subdirectory.
Their successful execution will lead to the generation of the LOBSTER dataset.

### 1. Ligand and Structure Preparation
* Use **1_reduce_ligands_main_filter.py** to extract ligands fitting our filter criteria from a "master table" of structure and ligand data.
* The extracted master table can be further reduced by removing duplicates with **2_reduce_ligands_duplicates.py**.
 Here, we also consider duplicates as ligands with the same USMILES occurring in the same PDB. Here, only the ligand with the best EDIAm is kept.
* Use **3_reduce_ligands_fetch_activities_for_table.py** to fetch activity data from PDB, annotate them at the master table, and remove all ligands with a LE < 0.3.
* The new, filtered master table is used in the next step to extract ligands for calculating the buriedness.
  * To do so, prepare the ligands and run **4_remove_ligands.py**.
  * To calculate the buriedness for the ligands after the previous steps, **5_calculate_buriedness.py**. This will execute the complex_properties from NAOMI and write a file with the requested buriedness-factor if the respective executable is provided.
  * **6_reduce_ligands_buriedness.py** uses the buriedness factor to filter the master table and write a file of all PDB IDs for the creation of the SIENA database in the next step.


### 2. SIENA Search
* The new, filtered master table is used in the next step to extract binding pockets, serving as queries for SIENA.
* Besides the pockets, SIENA needs a database filled with the ligands surviving our filter criteria.
 For only including the PDBs with the surviving ligands in our database, we create a directory of symlinks to those PDBs from the mirror.
 This directory is used as input for the generate_siena_database tool of the NAOMI Chembio Suite.
 The script **7_generate_siena_db.py** performs these steps automatically.
* Having our pockets as queries and the database to search in, a call of **8_execute_siena.py** will perform the SIENA search.
 In this script, for each result, SIENA offers for a query, we post-process the result within the script to guarantee that only ligands occur in our results,
 that did pass the filter, even though there might be filtered-out ligands in the symlinked PDBs of the previous step.

### 3. Final LOBSTER Creation
* Having finished the SIENA search, the resulting ensembles might contain ligands with the same USMILES which were taken from different PDBs.
 Again we remove the duplicates and keep only one of these ligands. This is done via **9_unique_ligand_extractor.py**.
  * backbone RMSD determines for the duplicates which ligand is kept.
  * Additionally, duplicate ensembles are removed, i.e. if the query ligand for SIENA was the same but
 from a different PDB. Keep ensemble with lower average backbone RMSD
* To cluster the ensembles and retrieve the set of representatives, **10_clustering** outputs our benchmark.
  * There are two clustering options for picking a representative: by lowest mean backbone RMSD or by largest query heavy atom count.
  * For the initial LOBSTER generation, query heavy atom count was chosen.

## Post Processing: Pairs and Subsets
The following scripts can be found in the **post_processing** subdirectory and help to split ensembles into the pairs and subsets of the LOBSTER dataset.
* To create the dataset of all pairs within the LOBSTER directory, execute **11_create_pairs.py**.
* Next, to provide the statistics on the Shape Tversky Index, please execute **stats/csv_stats.py**.
* Now, a call of **12_create_subsets.py** will generate the subsets dataset within the LOBSTER directory.

## Evaluation

### 1. Statistics on the Dataset
The following scripts can be found in the **stats** subdirectory and provide statistics on the ensembles, pairs, or subsets of the LOBSTER dataset.
* **csv_stats.py** will generate statistics about the molecules, pairs, and ensembles.
* **size_stats.py** will generate violin plots within the plots directory specified in the config.
* **physchem_stats.py** will generate plots about cLogP, molecular weight, rotatable bonds, and H-bond donors/acceptors. It also writes a file containing all PDB IDs of molecules within LOBSTER.
* **tanimoto_boxplots.py** requires the stats generated by **csv_stats.py** and generates boxplots for the ECFP4 and Gobbi fingerprint Tanimoto similarities.
* **diversity_stats.py** generates statistics about median MACCS fingerprint similarity and Scaffold AUC as input for consensus diversity plots (CDPs). These are stored as CSV files for the subsets and ensembles of LOBSTER as well as for the AstraZeneca Overlays dataset if specified in the config. Important note: execute **12_create_subsets.py** first to provide subsets!
* **plot_diversity_stats.py** plots the results generated in **diversity_stats.py**.
* **generate_timestamps.py** generates a CSV that assigns the Revdat to each PDB ID according to the file of all PDB IDs generated in **physchem_stats.py**.
* **plot_timestamps.py** plots the timestamps stored in the CSV previously.
* **comparison_to_drugs_histograms.py** plots the histograms for the comparison of molecular properties from the LOBSTER compounds and a list of known orally available drugs from the FDA orange book as provided with the LOBSTER dataset.

### 2. Studies About Protein Diversity
The following scripts deal with the diversity of the molecules from LOBSTER and the AstraZeneca Overlays dataset and can be found in the subdirectory **protein_diversity**.
* **retrieve_chain_uniprot_ids.py** creates CSV files for all molecules in the LOBSTER dataset with information on the cluster, PDB ID, ligand hetcode, chain, Uniprot ID, PfamID, and Pfam Name.
* **retrieve_uniprot_az.py** creates an equal CSV for the AstraZeneca Overlays as it was created in the above-mentioned script for LOBSTER.
* **check_and_compare_uniprots.py** generates the plots and information provided in the manuscript about PFam and Uniprot IDs.
* **analyze_double_uniprots.py** creates a CSV of all proteins with double Uniprot IDs and all that fail during the analysis.

## Python Dependencies
* Python 3.10.11
* pandas=1.4.1
* rdkit=2022.0 9.1
* progressbar=4.0.0
* argparse=1.1
* json=2.0.9
* logging=0.5.1.2
* sklearn=1.1.3
* matplotlib=3.5.1
* plotnine=0.10.1
* seaborn=0.11.2
* matplotlib-venn=1.0.0
* requests=2.28.1
* re=2.2.1
* scipy=1.8.0

## Open Source License Information
Copyright (c) 2024

University of Hamburg, ZBH - Center for Bioinformatics

Sophia Hönig, Torben Gutermuth, Christiane Ehrt, Christian Lemmen, Matthias Rarey

All files of the workflow and its analysis are licensed by the BSD New license (see LICENSE file for more information).
# Starting off with search for genomes in a big file
import csv # Allows reading .tsv or .csv files.
from Bio import Entrez # Imported but not used. Could potentially be useful to query NCBI databases.
import time # Used to pause between genome downloads to avoid overwhelming the NCBI server.
import os # Enables shell commands (like downloading data via command line). Mostly used to move easily through folders here.


# Define target methanotrophic genera
target_genera = ["Methylobacter", "Crenothrix"] # "Methylobacterium" is also included in the search because the basis of the search is (in part) "Methylobacter" 

# Function to read TSV and extract accession numbers for target genera
def extract_accession_numbers(tsv_file):
	accession_numbers = []
	with open(tsv_file, 'r', encoding='utf-8') as file:
		reader = csv.reader(file, delimiter='\t')
		for row in reader:
			if len(row) >= 2:
				gtdb_accession = row[0]  # Accession number
				accession = "_".join(gtdb_accession.split("_")[1:])
				classification = row[1]  # Taxonomy string
				for genus in target_genera:
					if f"g__{genus}" in classification:
						accession_numbers.append(accession)
						print(accession)
						break
	return accession_numbers

acc_nr = extract_accession_numbers("bac120_taxonomy.tsv")

# Download datasets with genome information on all target genera
for acc in acc_nr:
	acc_filename = acc + ".zip"
	os.system(f".\datasets.exe download genome accession {acc} --filename {acc_filename} --include gff3,rna,cds,protein,genome,seq-report")
	time.sleep(1)

# If all downloaded data should contain .fna files then: --include gff3,rna,cds,protein,genome,seq-report

# If any problems arose during the download of a zipped files containing .fna file for the respective accession number, then
# the following line was executed in the terminal (not python): 
# Get specific zip file with accession number (example) GCF_030161815.1.zip:
.\datasets download genome accession GCF_030161815.1 --include gff3,rna,cds,protein,genome,seq-report 

### ALL ZIP-FILES WERE MOVED TO A FOLDER NAMED "all_genomes" ###


## List all .fna files
import os # Let me traverse directories 
import zipfile # To extract compressed .zip files

genome_folder = "C:\\Users\\emilf\\OneDrive\\Dokumenter\\UNI\\Bachelor\\all_genomes"

for file in os.listdir(genome_folder):
    if file.endswith(".zip"):
        zip_path = os.path.join(genome_folder, file)
        extract_path = os.path.join(genome_folder, "extracted_files", file.replace(".zip", "")) # Target folder to extract each genome into a subfolder with the same name
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            zip_ref.extractall(extract_path)
        print(f"Extracted: {file}") # confirmation message for each extracted file.


## Translate genome (.fna) to protein using prodigal
import subprocess # Enables to run external programs (like prodigal) from within Python

extracted_folder = "C:\\Users\\emilf\\OneDrive\\Dokumenter\\UNI\\Bachelor\\all_genomes\\extracted_files"

# Loop through each extracted genome folder 
for folder in os.listdir(extracted_folder):
    folder_path = os.path.join(extracted_folder, folder, "ncbi_dataset", "data", folder)
    
    # Look for a .fna file (genome sequence)
    for file in os.listdir(folder_path):
        if file.endswith(".fna"):  # For the files with .fna extension in the folders
            print(file) # To confirm the file
            genome_file = os.path.join(folder_path, file)
            protein_output = os.path.join("translated_proteins", folder + ".faa") # This defines where the translated protein FASTA file will be saved (in the "translated_proteins" folder). 
            # Each output file is named after the genome folder
            print(genome_file)
            print(protein_output)
            
            # Run Prodigal to translate genome to protein sequences
            cmd = f".\\prodigal.windows.exe -a {protein_output} -i {genome_file}"
            subprocess.run(cmd, shell=True)
            print(f"Translated: {file} -> {protein_output}") # Status message for user feedback

# Correct KB to acc_KB (kun i specifikt tilfælde)?
# This code corrects the protein header format for a single genome file that did not follow usual naming patterns
protein_file = "translated_proteins\\GCA_000333655.1.faa"
genome_name = protein_file.split("\\")[1].split(".")[0] 

# Modifies header lines to include the genome accession as a prefix.
for line in open(protein_file).readlines():
    if line.startswith(">"):
        newline = line.replace(">", ">" + genome_name + "_")
        print(newline)
    else:
        print(line) 

### Check the protein files ###

translated_protein_folder_path = "C:\\Users\\emilf\\OneDrive\\Dokumenter\\UNI\\Bachelor\\translated_proteins"

# Check length of translated protein folder path
len(os.listdir(translated_protein_folder_path))

# For loop iterating over every protein (.faa) file in the translated_proteins folder
# Make sure the files are named and saved properly.
for protein_file in os.listdir(translated_protein_folder_path):
    if protein_file.endswith(".faa"):
        print(protein_file)
# Check over

## Rename and merge all lines from all files to one file using os module
# This is done to combine all protein FASTA files into one, making sure each protein ID is prepended with the genome name, so that every protein is traceable.

# Rename
with open("C:\\Users\\emilf\\OneDrive\\Dokumenter\\UNI\\Bachelor\\all_proteins.faa", 'w') as outfile:
    for filename in os.listdir(translated_protein_folder_path):
        filename_fullpath = os.path.join(translated_protein_folder_path, filename)
        print(filename_fullpath)
        infile = open(filename_fullpath, 'r')
        genome_name = filename.split(".")[0]
        for line in infile.readlines():
            if line.startswith(">"):
                newline = line.replace(">", ">" + genome_name + "_")
                outfile.write(newline)
            else:
                outfile.write(line)


## Make a BLAST database and make a BLAST search for query proteins. This was executed in the terminal
& 'C:\Program Files\NCBI\blast-2.16.0+\bin\makeblastdb.exe' -dbtype prot -in all_proteins.faa

# Make BLAST query against all related proteins (related proteins were selected from table 2 in Castelle et al. 2008). This was also executed in the terminal
& 'C:\Program Files\NCBI\blast-2.16.0+\bin\blastp.exe' -query related_proteins\all_related_proteins.faa -db all_proteins.faa -outfmt 6 -out blast_output.txt
# I linux fra blastp -... (uden '')


## Read the blast output
pip install pandas # Executed in the terminal
import pandas # Enables manipulation of tables

# The actual BLAST output file
blast_output_file = "blast_output.txt"  

# E-value threshold (run with 1e-10)
e_value_cutoff = 1e-10

df = pandas.read_table(blast_output_file, sep = "\t")

# Creating table with columns containing taxa (rows with accession numbers) and columns with genes (rows with e-value and bit score)
col_names = [
    "query_id", "subject_id", "percent_identity", "alignment_length",
    "mismatches", "gap_opens", "q_start", "q_end", "s_start", "s_end",
    "e_value", "bit_score"
]


# Read BLAST output file (again). This re-loads the same file here but now with proper column names. The earlier df is replaced
df = pandas.read_csv(blast_output_file, sep="\t", names= col_names)

# Filter results based on E-value cutoff
df_filtered = df[df["e_value"] <= e_value_cutoff] # Filters for hits that have an e-value of significance 

# Save the filtered results to a CSV file for backup or external inspection.
table_output_csv = "blast_summary_table.csv"
df_filtered.to_csv(table_output_csv) 

print(df_filtered) # View/inspect the filtered df

# Returns a set (unique values only) of all query IDs that had at least one hit in the BLAST results. This confirms which query genes matched genomes
set(df_filtered["query_id"])

# Exports filtered table to Excel format for easy viewing
df_filtered.to_excel("blast_results.xlsx")

df_filtered # Execute line to view df_filtered in the console/terminal

# Creates list (initially empty) 
genome_list = list()

# Iterate over every genome in the "subject_id" column in the df_filtered list and append it to a list
for id in df_filtered["subject_id"]:
    genome = "_".join(id.split("_")[0:2])
    genome_list.append(genome)

genome_list # Execute to view the genome list
# Check length of genome list (some genomes occur more than once here though)
len(genome_list) # 1101

# Append the "genome_list" to a new column in df_filtered and call that column "genome"
df_filtered["genome"] = genome_list 

df_filtered.unstack() # Not used i think (maybe run this later)

# Iterate over every e-value in df_filtered
for value in df_filtered["e_value"]:
    print(value)

# Iterate over every bit-score in df_filtered
for score in df_filtered["bit_score"]:
    print(score)


# Pivot the table to have genomes as rows and query proteins as columns
df_pivot = df_filtered.pivot_table(index="genome", columns="query_id", values="e_value", aggfunc="min")

# Reset the index to make 'subject_id' a regular column
df_pivot.reset_index(inplace=True)

# Export to excel
df_pivot.to_excel("formatted_new_table.xlsx", index=False)

# Export to csv
df_pivot.to_csv("formatted_table.csv")

# Proof that genomes do not occur multiple times in the pivoted dataframe (df) 
len(df_pivot['genome']) # (check maybe)

g_list = list()

for data in df_pivot['genome']:
    if data not in g_list:
        g_list.append(data)

len(g_list) # Prints the same length 

### On the dnaseq1@.au.dk server
# THEN THE fegenie.sh FILE CONTENT WAS EXECUTED AND RESULTED IN A fegenie-summery FILE, CONTAINING HMMS FOR IRON-RELATED GENES IN OUR TARGET GENERA

### PHYLOGENY ###

from ete3 import Tree, NodeStyle, TextFace, TreeStyle, CircleFace, faces

# Initialize empty dictionary
species_dict = dict()

with open("bac120_taxonomy.tsv", "r") as file:
    for line in file:
        
        # Split lines by tab ("\t") character
        parts = line.strip().split("\t")
        
        genome_id = parts[0].split("_")
        genome_id = "_".join(genome_id[1:])
        
        species = parts[1].split(";")[-1]
        species = species.split("__")[1]
        # print(species)
        # print(genome_id)
        species_dict[genome_id] = species

species_dict # Check and prints the dictionary

# Collecting gene-presence sets as a test (OPTIONAL)
fegenie = open("genie\\FeGenie-geneSummary.csv")

# Prepares an empty set to store unique genome IDs that had genes in the "iron_oxidation" category
iron_oxidation_genomes = set()

# Iterates over every line in the file and checks if the line belong to the "iron_oxidation" category and if FeGenie outputs start with the category in each line
for line in fegenie:
    if line.startswith("iron_oxidation"):
        iron_oxidation_genomes.add(line.split(',')[1].strip(".fna")) # Strip the .fna file extension and add the genome ID to the set

iron_oxidation_genomes # Prints the final set of matched genomes to iron_oxidation related genes

### FIND SPECIFIC GENES ###
import csv

# Path to your FeGenie output
fe_file = "genie/FeGenie-geneSummary.csv"

## LIST RESPECTIVE GENOMES TO THEIR HMMS

genomes_to_hmm = {}

with open(fe_file, newline='') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        category = row["category"]
        gene   = row["HMM"] # The HMM (homologs) column
        genome = row["genome/assembly"].replace(".fna", "") # strips the file extension
        genome = genome.split("_")
        genome = "_".join(genome[0:2])
        
        # initialize a set if this is the first time we see this genome
        if genome not in genomes_to_hmm:
            genomes_to_hmm[genome] = set()
        
        # add the genes to that genomes set
        genomes_to_hmm[genome].add(gene)


genomes_to_hmm # Check which genomes have which HMM gene hits

### CHECK ###
# Check how many unique gene types are present (optional)
generinos = set()
for genes in genomes_to_hmm.values():
    for gene in genes:
        if gene not in generinos:
            generinos.add(gene)
    
len(generinos) # Same -79 genes

# ALSO GET RESPECTIVE GENES IN CERTAIN CATEGORIES (optional)

# Load the fegenie file
fegenie = open("genie\\FeGenie-geneSummary.csv")

irg = set()

for line in fegenie:
    if line.startswith("iron_oxidation"):
        irg.add(line.split(',')[3])
    if line.startswith("possible_iron_oxidation_and_possible_iron_reduction"):
        irg.add(line.split(',')[3])
    if line.startswith("probable_iron_reduction"):
        irg.add(line.split(',')[3])

irg # Prints the different genes within the 3 categories -encoded within all included genomes
# CHECK OVER 

## ANNOTATION

# Mapping from blast protein accession to a readable name
protein_name_map = {
    "CAA07038.1": "rusticyanin",
    "CAA07031.1": "Cyt c",
    "CAA07034.1": "Cyt ox II",
    "CAA07035.1": "Cyt ox I",
    "CAA07033.2": "hypothetical protein"
}

# INCLUDE BLAST DATA 

import pandas

# Load the Excel table
blast_df = pandas.read_excel("formatted_new_table.xlsx")

# Set the threshold again, to be sure of significant hits
threshold = 1e-10

# Create a dict: genome to list of protein IDs with significant hits
blast_hits = {}

for _, row in blast_df.iterrows():
    genome = row["genome"]
    hits = []
    for protein in blast_df.columns[1:]:
        try:
            evalue_str = str(row[protein]).replace(",", ".")  # convert comma to dot
            evalue = float(evalue_str)
            if evalue <= threshold:
                hits.append(protein)
        except:
            continue
    blast_hits[genome] = hits

blast_hits # Shows genome_accession as a key and the value pair is the iron-related genes of which there is a significant hit (e-value 1e-10) on.


## Creating the phylogenetic tree (following the protocol step 1-3, and 8-13)
the_tree = "fylogeni/gtdbtk_fasttree_phylogeny.tree"

# Loading the tree
t = Tree(the_tree)

# Annotating leaves
for leaf in t.iter_leaves():
    genome_accession = "_".join(leaf.name.split("_")[0:2])
    print(genome_accession)
    # 5.1 Add species name
    if genome_accession in species_dict.keys():
        leaf.add_face(TextFace(species_dict[genome_accession]), 1, position="aligned")
    # 5.2 Mark MHC-positive genomes
    if "Cyc2_repCluster2" in genomes_to_hmm[genome_accession]:
        leaf.add_face(TextFace("Cyc2", fgcolor="green"), 2, position="aligned")
    # 5.3 Mark MtrB-positive genomes
    if "MtrB_TIGR03509" in genomes_to_hmm[genome_accession]:
        leaf.add_face(TextFace("MtrB", fgcolor="red"), 3, position="aligned")
    # 5.4 Mark iron-oxidation (Cyc2) genomes
    if "MtoA" in genomes_to_hmm[genome_accession]:
        leaf.add_face(TextFace("MtoA", fgcolor="purple"), 4, position="aligned")
    if "Cyc1" in genomes_to_hmm[genome_accession]:
        leaf.add_face(TextFace("Cyc1", fgcolor="olive"), 5, position="aligned")
    if "MtrA" in genomes_to_hmm[genome_accession]:
        leaf.add_face(TextFace("MtrA", fgcolor="red"), 6, position="aligned")
    # Add blast hits
    short_genome = genome_accession.split(".")[0]  # e.g., GCA_000526475.1 → GCA_000526475. This is done because genome_acc... ends with '.1'
    if short_genome in blast_hits:
        hits = blast_hits[short_genome]
        nr = 0
        for hit in hits:
            nr += 1
            leaf.add_face(TextFace(f"Hit: {protein_name_map[hit]}   ", fgcolor="blue"), column=(6+nr), position="aligned")


# Tidying up
t.ladderize()          # Sorts the branches to make the tree look “neater”
t.show()               # Opens an interactive window (if you’re on a GUI-enabled machine)
t.render("second_hmm_and_blast_annotated_tree_v4.pdf", w=600, units="mm")  # Exports to a PDF file

### EXCLUDE METHYLOBACTERIUM ###

# For loop iterating over the name in every leaf and detaches leaves with the name "Methylobacterium" (optinal)
for leaf in t.iter_leaves():
    genome_accession = "_".join(leaf.name.split("_")[0:2])
    if genome_accession in species_dict:
        species_name = species_dict[genome_accession]
        if "Methylobacterium" in species_name:
            leaf.detach()

# Run again:
t.show()               # Opens an interactive window (if you’re on a GUI-enabled machine)
t.render("final_hmm_and_blast_annotated_tree_v4.pdf", w=600, units="mm")  # Exports to a PDF file


### THE END



# number of respective hits in HMM genomes

fegenie = open("genie\\FeGenie-geneSummary.csv")

hmm_gene_counts ={"iron_oxidation": 0, "possible_iron_oxidation_and_possible_iron_reduction": 0, "probable_iron_reduction": 0}

for line in fegenie:
    if line.startswith("iron_oxidation"):
        hmm_gene_counts["iron_oxidation"] += 1
    if line.startswith("possible_iron_oxidation_and_possible_iron_reduction"):
        hmm_gene_counts["possible_iron_oxidation_and_possible_iron_reduction"] += 1
    if line.startswith("probable_iron_reduction"):
        hmm_gene_counts["probable_iron_reduction"] += 1

hmm_gene_counts


respective_hmm_homolog_count = {'Cyc2_repCluster2': 0, 'MtoA': 0, 'MtrA': 0, 'Cyc1': 0, 'MtrB_TIGR03509': 0}

for line in fegenie:
    if line.startswith("iron_oxidation"):
        if line.split(',')[3] == "Cyc1":
            respective_hmm_homolog_count['Cyc1'] += 1
        elif line.split(',')[3] == "Cyc2_repCluster2":
            respective_hmm_homolog_count["Cyc2_repCluster2"] += 1
    if line.startswith("possible_iron_oxidation_and_possible_iron_reduction"):
        if line.split(',')[3] == "MtoA":
            respective_hmm_homolog_count['MtoA'] += 1
        elif line.split(',')[3] == "MtrB_TIGR03509":
            respective_hmm_homolog_count["MtrB_TIGR03509"] += 1
    if line.startswith("probable_iron_reduction"):
        if line.split(',')[3] == "MtrA":
                respective_hmm_homolog_count['MtrA'] += 1
        elif line.split(',')[3] == "MtrB_TIGR03509":
            respective_hmm_homolog_count["MtrB_TIGR03509"] += 1

# MtrB = 96, MtoA = 95, Cyc2 = 102, MtrA = 1, and Cyc1 = 3
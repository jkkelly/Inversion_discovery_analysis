"""
This script will automate the use of anchorwave per chromosome, once the files have been preprocessed

@Author: Luis J. Madrigal-Roca & John K. Kelly
@Date: 2024-08-7

"""

import os
import argparse

parser = argparse.ArgumentParser(description = "Anchorwave inversion extractor")

parser.add_argument("-rd", "--reference_directory", help = "Reference directory. This directory must contain reference fasta file and annotation file", required = True)
parser.add_argument("-qd", "--query_directory", help = "Query directory. This directory must contain query fasta files", required = True)
parser.add_argument("-ad", "--anchorwave_directory", help = "Anchorwave working directory. This directory will store temporary files", required = True)
parser.add_argument("-j", "--jobs_directory", help = "Jobs directory. This directory will store the bash scripts for each chromosome", required = True)
parser.add_argument("-r", "--results_directory", help = "Results directory. This directory will store the results of the anchorwave alignments", required = True)
parser.add_argument("-i", "--intermediate_directory", help = "Intermediate directory. This directory will store the intermediate files", required = True)
parser.add_argument("-a", "--anchorwave_executable", help = "Anchorwave executable", required = False, default = '/home/l338m483/bin/anchorwave/anchorwave')
parser.add_argument("-m", "--minimap2_executable", help = "Minimap2 executable", required = False, default = '/home/l338m483/bin/minimap2/minimap2')

args = parser.parse_args()

REF_DIR = args.reference_directory
QUE_DIR = args.query_directory
ANC_DIR = args.anchorwave_directory
JOBS = args.jobs_directory
RES = args.results_directory
INT_DIR = args.intermediate_directory
anchorwave = args.anchorwave_executable
minimap2 = args.minimap2_executable

ref_name = None
ref_gff = None

for file in os.listdir(REF_DIR):
    if file.endswith(".fasta") or file.endswith(".fa"):
        ref_name = file
    elif file.endswith(".gff3"):
        ref_gff = file

common_headers_files = []

for filename in os.listdir(INT_DIR):
    if "common" in filename:
        common_headers_files.append(filename)

for query in os.listdir(QUE_DIR):
    
    for file in common_headers_files:
        if query in file:
            current_header_file = file
    
    with open(f"{os.path.join(INT_DIR, current_header_file)}", "r") as file:
        headers = [header.strip().replace(">", "") for header in file.readlines()]

    for head in headers:

        # Create a bash script for each file
        with open(f"{os.path.join(JOBS,f'anchorwave_{query}_{head}_alig_job.sh')}", "w") as file:
            file.write("#!/bin/bash\n")
            file.write(f"#SBATCH --job-name=anchorwave_{query}_{head}_job\n")  # Job name
            file.write(f"#SBATCH --output=anchorwave_{query}_{head}_output\n")  # Output file name
            file.write("#SBATCH --partition=sixhour\n") # Work partition
            file.write("#SBATCH --nodes=1\n")  # Number of nodes
            file.write("#SBATCH --ntasks=1\n") # Number of tasks
            file.write("#SBATCH --constraint=avx512\n")
            file.write("#SBATCH --cpus-per-task=10\n") # Number of parallel processes 
            file.write("#SBATCH --time=6:00:00\n")  # Time
            file.write("#SBATCH --mail-user=l338m483@ku.edu\n")  # Mail user
            file.write("#SBATCH --mail-type=FAIL\n") # Mail type
            file.write("#SBATCH --mem-per-cpu=10g\n")  # Memory limit
            file.write("\n")
            file.write("module load samtools\n")
            file.write(f"cd {ANC_DIR}\n")
            file.write("\n")
            file.write(f"mkdir TEMP_{query}_{head}/\n")
            file.write(f"cp {os.path.join(REF_DIR, ref_name)} TEMP_{query}_{head}\n")
            file.write(f"cp {os.path.join(QUE_DIR, query)} TEMP_{query}_{head}\n")
            file.write("\n")
            file.write(f"grep '##' {os.path.join(REF_DIR, ref_gff)} > TEMP_{query}_{head}/Filtered_{ref_gff}\n")
            file.write(f"grep {head} {os.path.join(REF_DIR, ref_gff)} >> TEMP_{query}_{head}/Filtered_{ref_gff}\n")
            file.write("\n")
            file.write(f"echo {head} > TEMP_{query}_{head}/region.txt\n")
            file.write("\n")           
            file.write(f"samtools faidx -r TEMP_{query}_{head}/region.txt -o TEMP_{query}_{head}/Filtered_{ref_name} TEMP_{query}_{head}/{ref_name}\n")
            file.write(f"samtools faidx -r TEMP_{query}_{head}/region.txt -o TEMP_{query}_{head}/Filtered_{query} {QUE_DIR}{query}\n")
            file.write("\n")
            file.write(f"{anchorwave} gff2seq -i TEMP_{query}_{head}/Filtered_{ref_gff} -r TEMP_{query}_{head}/Filtered_{ref_name} -o TEMP_{query}_{head}/Filtered_{ref_name}.cds\n")
            file.write("\n")
            file.write(f"{minimap2} -x splice -t 10 -k 12 -a -p 0.4 -N 20 TEMP_{query}_{head}/Filtered_{ref_name} TEMP_{query}_{head}/Filtered_{ref_name}.cds > TEMP_{query}_{head}/{ref_name}.sam\n")
            file.write(f"{minimap2} -x splice -t 10 -k 12 -a -p 0.4 -N 20 TEMP_{query}_{head}/Filtered_{query} TEMP_{query}_{head}/Filtered_{ref_name}.cds > TEMP_{query}_{head}/{query}.sam\n")
            file.write("\n")
            file.write(f"{anchorwave} genoAli -t 10 -I 1 -i TEMP_{query}_{head}/Filtered_{ref_gff} -as TEMP_{query}_{head}/Filtered_{ref_name}.cds -r TEMP_{query}_{head}/Filtered_{ref_name} -a TEMP_{query}_{head}/{query}.sam -ar TEMP_{query}_{head}/{ref_name}.sam -s TEMP_{query}_{head}/Filtered_{query} -n {RES}{query}_{head}_anchors -o {RES}{query}_{head}_anchorwave.maf -f {RES}{query}_{head}_anchorwave.f.maf -IV\n")

        # Make the script executable
        os.system(f"chmod a+x {os.path.join(JOBS,f'anchorwave_{query}_{head}_alig_job.sh')}")
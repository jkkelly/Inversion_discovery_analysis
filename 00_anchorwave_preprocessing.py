"""
This script automates preprocessing of input files for Anchorwave implementations.

**Author:** Luis Javier Madrigal-Roca & John K. Kelly
**Date:** 2024-08-07

**Description:**
- Iterates through query chromosome files.
- Creates individual job scripts for each query to perform preprocessing tasks.

**Requirements:**
- samtools

**Notes:**
- Script assumes samtools is installed and available in the system path.
"""

import os
import argparse

parser = argparse.ArgumentParser(description = "Anchorwave preprocessing")

parser.add_argument("-rd", "--reference-directory", help="Reference directory where the fasta and gff3 for the reference are stored", required=True)
parser.add_argument("-qd", "--query-directory", help="Query directory where the fasta files for the query(ies) are stored", required=True)
parser.add_argument("-awd", "--anchorwave-working-directory", help="Directory where the anchorwave temporal files are stored", required=True)
parser.add_argument("-j", "--jobs-directory", help="Directory where the jobs are stored", required=True)
parser.add_argument("-id", "--intermediate-directory", help="Directory where the intermediate files will be stored stored. These are the outputs of this program", required=True)

args = parser.parse_args()

REF_DIR = args.reference_directory
QUE_DIR = args.query_directory
ANC_DIR = args.anchorwave_working_directory
JOBS = args.jobs_directory
INT_DIR = args.intermediate_directory

queries = os.listdir(QUE_DIR)

ref_name = None
ref_gff = None

for file in os.listdir(REF_DIR):
    if file.endswith(".fasta") or file.endswith(".fa"):
        ref_name = file
    elif file.endswith(".gff3"):
        ref_gff = file

if ref_name is None or ref_gff is None:
    raise ValueError("Reference fasta or gff3 file not found in the reference directory")

for query in queries: 

    # Create a bash script for each file
    job_script_path = os.path.join(JOBS, f"pre_anchorwave_{query}_alig_job.sh")
    with open(job_script_path, "w") as file:
        file.write("#!/bin/bash\n")
        file.write(f"#SBATCH --job-name=pre_anchorwave_{query}_job\n")  # Job name
        file.write(f"#SBATCH --output=pre_anchorwave_{query}_output\n")  # Output file name
        file.write("#SBATCH --partition=sixhour\n") # Work partition
        file.write("#SBATCH --nodes=1\n")  # Number of nodes
        file.write("#SBATCH --ntasks=1\n") # Number of tasks
        file.write("#SBATCH --cpus-per-task=1\n") # Number of parallel processes 
        file.write("#SBATCH --time=6:00:00\n")  # Time
        file.write("#SBATCH --mail-user=l338m483@ku.edu\n")  # Mail user
        file.write("#SBATCH --mail-type=FAIL\n") # Mail type
        file.write("#SBATCH --mem-per-cpu=2g\n")  # Memory limit
        file.write("\n")
        file.write("module load samtools\n")
        file.write(f"cd {ANC_DIR}\n")
        file.write("\n")
        ## Filtering the fasta files for working only with common contigs
        ref_header_path = os.path.join(INT_DIR, f"{query}_against_header_reference.txt")
        query_header_path = os.path.join(INT_DIR, f"header_{query}.txt")
        common_headers_path = os.path.join(INT_DIR, f"common_for_{query}_headers.txt")
        
        file.write(f'grep ">" {os.path.join(REF_DIR, ref_name)} > {ref_header_path}\n') # Extract the headers of the reference file
        file.write(f'grep ">" {os.path.join(QUE_DIR, query)} > {query_header_path}\n')
        file.write("\n")
        ## Finding common contigs
        file.write(r"awk 'NR==FNR{a[$0];next}($0 in a)' ")
        file.write(f"{ref_header_path} {query_header_path} > {common_headers_path}\n")

    # Make the script executable
    os.system(f"chmod a+x {job_script_path}")
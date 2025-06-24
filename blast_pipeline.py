#!/usr/bin/env python3
"""
BLASTn Remote Pipeline

This script performs the following:
1. Reads input FASTA files from the 'files-fasta' directory
2. Optionally splits large files into chunks
3. Cleans sequences (removes gaps)
4. Runs remote BLASTn against NCBI's nt database
5. Extracts metadata from BLAST output and saves to Excel/TSV
6. Fetches top hit sequences in FASTA format
7. Merges all results into unified output files

⚠️ WARNING: Avoid using this script with hundreds of sequences unless you have an API key.
NCBI may temporarily block your IP due to excessive access.

Author: João Paulo Romualdo Alarcão Bernardes
License: MIT
"""

import os
import re
import sys
import time
import subprocess
import pandas as pd
from Bio import SeqIO, Entrez
from Bio.Seq import Seq

# ======================= USER CONFIGURATION =======================

# Set Entrez email (required by NCBI)
Entrez.email = os.getenv("ENTREZ_EMAIL", "example@gmail.com")  # Use your own email here

# BLAST parameters
database = "nt"                # Remote BLAST database
evalue_cutoff = 1e-20          # E-value threshold for reporting
num_hits = 5                   # Max hits returned per sequence

# Sequence chunking to avoid NCBI IP block (This is a polite measure) - temp folder will be used to store split files
max_seqs_per_batch = 50 # Recommended to keep it low to avoid NCBI IP block (50-100 is a common value)
# Eg: if you have 100 sequences and the max_seqs_per_batch = 50, it will create 2 files with 50 sequences each inside the temp folder.

# ======================= PATH SETUP ==============================

script_dir = os.path.dirname(os.path.abspath(__file__))
fasta_folder = os.path.join(script_dir, "files-fasta")
output_folder = os.path.join(fasta_folder, "outputs")
temp_folder = os.path.join(script_dir, "temp")

for folder in [output_folder, temp_folder]:
    os.makedirs(folder, exist_ok=True)

# ======================= HELPER FUNCTIONS ========================

def format_time(seconds):
    if seconds < 60:
        return f"{seconds:.1f} seconds"
    elif seconds < 3600:
        return f"{seconds/60:.1f} minutes"
    return f"{seconds/3600:.1f} hours"

def split_fasta_file(fasta_path, max_seqs=25):
    records = list(SeqIO.parse(fasta_path, "fasta"))
    if len(records) <= max_seqs:
        return [fasta_path]
    base = os.path.splitext(os.path.basename(fasta_path))[0]
    paths = []
    for i in range(0, len(records), max_seqs):
        chunk = records[i:i+max_seqs]
        part = f"{base}_part{i//max_seqs+1}.fasta"
        path = os.path.join(temp_folder, part)
        SeqIO.write(chunk, path, "fasta")
        paths.append(path)
    return paths

# ======================= CHECK DEPENDENCIES =======================

if subprocess.run(["which", "blastn"], capture_output=True).returncode != 0:
    print("❌ Error: BLAST+ tools are not installed or not in PATH.")
    sys.exit(1)

# ======================= STARTING PIPELINE ========================

print("=" * 70)
print("🧬 BLAST ANALYSIS PIPELINE - START")
print("=" * 70)
print(f"✅ Using Entrez email: {Entrez.email}")

fasta_files = sorted(f for f in os.listdir(fasta_folder) if f.endswith((".fasta", ".fa")))
if not fasta_files:
    print("❌ No FASTA files found.")
    sys.exit(1)

# === Get first base name for naming output files ===
first_base_name = os.path.splitext(fasta_files[0])[0]

# ======================= MAIN LOOP ================================

for fasta_file in fasta_files:
    full_path = os.path.join(fasta_folder, fasta_file)
    split_files = split_fasta_file(full_path, max_seqs_per_batch)

    for split_fasta in split_files:
        base = os.path.splitext(os.path.basename(split_fasta))[0]
        cleaned_fasta = os.path.join(temp_folder, f"{base}_cleaned.fasta")
        blast_output = os.path.join(temp_folder, f"{base}_blast.tsv")
        excel_output = os.path.join(temp_folder, f"{base}_blast.xlsx")
        fasta_output = os.path.join(temp_folder, f"{base}_hits.fasta")

        print(f"\n🔹 Processing {split_fasta}")

        # STEP 1: Clean sequences
        print("📋 STEP 1: Cleaning sequences...")
        with open(cleaned_fasta, "w") as out_handle:
            count = 0
            for record in SeqIO.parse(split_fasta, "fasta"):
                record.seq = Seq(str(record.seq).replace("-", ""))
                SeqIO.write(record, out_handle, "fasta")
                count += 1
        if count == 0:
            print(f"⚠️ Skipping {split_fasta} (no valid sequences)")
            continue

        print(f"✅ Cleaned {count} sequences")

        # STEP 2: Run BLASTn remotely
        print("📋 STEP 2: Running remote BLAST...")
        print(f"🔍 Using database={database}, e-value={evalue_cutoff}, max hits={num_hits}")
        start = time.time()
        blast_cmd = [
            "blastn", "-query", cleaned_fasta,
            "-db", database, "-remote",
            "-outfmt", "6 qseqid sacc pident length evalue bitscore sseqid stitle",
            "-max_target_seqs", str(num_hits),
            "-evalue", str(evalue_cutoff),
            "-out", blast_output
        ]
        try:
            time.sleep(5)  # Politeness delay
            subprocess.run(blast_cmd, check=True)
            print(f"✅ BLAST completed in {format_time(time.time() - start)}")
        except subprocess.CalledProcessError as e:
            print(f"❌ BLAST failed: {e}")
            continue

        # STEP 3: Parse results
        print("📋 STEP 3: Parsing BLAST output...")
        try:
            cols = ["Query", "Accession", "Identity", "Length", "E-value", "Bit Score", "SeqID", "Definition"]
            df = pd.read_csv(blast_output, sep="\t", names=cols)
            df.to_excel(excel_output, index=False)
            print(f"✅ Saved results to {excel_output}")
        except Exception as e:
            print(f"❌ Error reading BLAST TSV: {e}")
            continue

        # STEP 4: Retrieve sequences
        print("📋 STEP 4: Downloading top hit sequences...")
        try:
            accessions = df["Accession"].unique().tolist()
            def chunks(lst, n): return (lst[i:i+n] for i in range(0, len(lst), n))
            with open(fasta_output, "w") as out_handle:
                for batch in chunks(accessions, 10):
                    ids = ",".join(batch)
                    try:
                        handle = Entrez.efetch(db="nucleotide", id=ids, rettype="fasta", retmode="text")
                        records = list(SeqIO.parse(handle, "fasta"))
                        handle.close()
                        SeqIO.write(records, out_handle, "fasta")
                    except Exception as err:
                        print(f"⚠️ Error retrieving batch {ids}: {err}")
            print(f"✅ FASTA sequences saved to {fasta_output}")
        except Exception as e:
            print(f"❌ Error during sequence retrieval: {e}")

# ======================= FINAL MERGE ===============================

print("\n📦 Merging results...")
merged_tsv = os.path.join(output_folder, f"{first_base_name}_combined_blast.tsv")
merged_xlsx = os.path.join(output_folder, f"{first_base_name}_combined_blast.xlsx")
merged_fasta = os.path.join(output_folder, f"{first_base_name}_combined_hits.fasta")

# Merge all .tsv
with open(merged_tsv, "w") as out:
    out.write("qseqid\tsacc\tpident\tlength\tevalue\tbitscore\tsseqid\tstitle\n")
    for f in sorted(os.listdir(temp_folder)):
        if f.endswith("_blast.tsv"):
            with open(os.path.join(temp_folder, f)) as part:
                out.writelines(part.readlines())

# Merge all Excel
sheets = []
for f in sorted(os.listdir(temp_folder)):
    if f.endswith("_blast.xlsx"):
        try:
            df = pd.read_excel(os.path.join(temp_folder, f))
            sheets.append(df)
        except:
            continue
if sheets:
    pd.concat(sheets, ignore_index=True).to_excel(merged_xlsx, index=False)

# Merge all FASTA
with open(merged_fasta, "w") as out:
    for f in sorted(os.listdir(temp_folder)):
        if f.endswith("_hits.fasta"):
            with open(os.path.join(temp_folder, f)) as part:
                out.writelines(part.readlines())

print(f"✅ Merged outputs saved in: \n- {merged_tsv}\n- {merged_xlsx}\n- {merged_fasta}")

# ======================= CLEANUP ===============================

print("\n🧹 Cleaning up temporary files...")
for f in os.listdir(temp_folder):
    try:
        os.remove(os.path.join(temp_folder, f))
    except:
        pass

print("\n✅ Pipeline complete.")
print("=" * 70)

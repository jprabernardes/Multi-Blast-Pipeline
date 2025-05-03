import re
import pandas as pd
import subprocess
import time
import sys
from Bio import SeqIO, Entrez
from Bio.Seq import Seq
try:
    from tqdm import tqdm
except ImportError:
    print("Installing tqdm library for progress bars...")
    subprocess.run([sys.executable, "-m", "pip", "install", "tqdm"], check=True)
    from tqdm import tqdm

# ============================================================
# === IMPORTANT: CONFIGURATION ==============================
# ============================================================

# Input/output files
input_fasta = "example/example.fasta" # Path to your input FASTA file
cleaned_fasta = "Sequence_Blasted_clean.fasta"
blast_output = "blast_remote.tsv" # Output file for BLAST results
excel_output = "blast_hits_remote.xlsx" # Output file for Excel results
database = "nt"  # you can change to 'refseq_rna' or another available
evalue_cutoff = 1e-20 # E-value cutoff for BLAST
num_hits = 5 # Maximum number of hits to retrieve

# === IMPORTANT: EMAIL CONFIGURATION FOR NCBI ===
# You MUST configure a valid email address below
# NCBI requires an email to use the Entrez API
# Sequence retrieval will not work without a valid email
Entrez.email = "example@example.com"  # Email configured for use with NCBI

# OPTIONAL: If you have an NCBI API key, uncomment and configure below
# An API key allows more requests per second
# Entrez.api_key = "your_api_key_here"

# File to save the sequences of the hits
blast_sequences = "blast_hit_sequences.fasta"

# Function to format estimated time
def format_time(seconds):
    if seconds < 60:
        return f"{seconds:.1f} seconds"
    elif seconds < 3600:
        return f"{seconds/60:.1f} minutes"
    else:
        return f"{seconds/3600:.1f} hours"

print("=" * 70)
print("🧬 BLAST ANALYSIS PIPELINE - START")
print("=" * 70)

# === 1. Clean dashes
print("\n📋 STEP 1: Preparing sequences...")

# Counting records for the progress bar
total_records = sum(1 for _ in SeqIO.parse(input_fasta, "fasta"))

with open(cleaned_fasta, "w") as out_handle:
    with tqdm(total=total_records, desc="Cleaning sequences", unit="seq") as pbar:
        for record in SeqIO.parse(input_fasta, "fasta"):
            # Replace the deprecated ungap() method with string replacement
            seq_str = str(record.seq).replace("-", "")
            record.seq = Seq(seq_str)
            SeqIO.write(record, out_handle, "fasta")
            pbar.update(1)

print(f"✅ Cleaned sequences saved in {cleaned_fasta} ({total_records} sequences processed)")

# === 2. Run remote BLAST against the full database
print(f"\n📋 STEP 2: Running remote BLAST against the {database} database")
print(f"⚠️  This step may take several minutes depending on the NCBI server")
print(f"🔍 Parameters: e-value={evalue_cutoff}, max hits={num_hits}")

start_time = time.time()

# Function to monitor the BLAST output file (to show progress)
def monitor_blast_progress():
    last_size = 0
    dots = 0
    while True:
        try:
            # Check if the output file exists
            current_size = 0
            try:
                current_size = os.path.getsize(blast_output)
            except:
                pass
                
            if current_size > last_size:
                # File is growing, BLAST is working
                progress = f"📊 Receiving BLAST data: {current_size/1024:.1f} KB downloaded"
                print(f"\r{progress}", end="")
                last_size = current_size
            else:
                # No update, just show it's working
                dots = (dots % 3) + 1
                print(f"\r🚀 BLAST running{'.' * dots}{' ' * 30}", end="")
            
            time.sleep(2)
        except KeyboardInterrupt:
            return

# Import threading for background monitoring
import threading
import os

# Start a thread to show progress
progress_thread = threading.Thread(target=monitor_blast_progress)
progress_thread.daemon = True
progress_thread.start()

# Run BLAST
blast_cmd = [
    "blastn",
    "-query", cleaned_fasta,
    "-db", database,
    "-remote",
    "-outfmt", "6 qseqid sacc pident length evalue bitscore sseqid stitle",
    "-max_target_seqs", str(num_hits),
    "-evalue", str(evalue_cutoff),
    "-out", blast_output
]

try:
    subprocess.run(blast_cmd, check=True)
    elapsed_time = time.time() - start_time
    print(f"\n✅ BLAST completed in {format_time(elapsed_time)}")
    
    # Check the size of the results file
    result_size = os.path.getsize(blast_output)
    print(f"📄 Results file: {blast_output} ({result_size/1024:.1f} KB)")
except subprocess.CalledProcessError as e:
    print(f"\n❌ Error running BLAST: {e}")
    sys.exit(1)
# === 3. Process output to Excel
print(f"\n📋 STEP 3: Processing BLAST results")

try:
    columns = ["Sequence-Blasted", "Accession", "Identity (%)", "Alignment Length", "E-value", "Bit Score", "SeqID", "Definition"]
    df = pd.read_csv(blast_output, sep="\t", names=columns)
    
    # Count results
    num_hits_found = len(df)
    num_sequence_blasted_with_hits = df["Sequence-Blasted"].nunique()
    
    print(f"📊 Found {num_hits_found} BLAST hits for {num_sequence_blasted_with_hits} the Sequence Blasted")
    
    # === 4. Save as Excel
    print(f"💾 Saving results to Excel...")
    df.to_excel(excel_output, index=False)
    print(f"✅ Excel saved as: {excel_output}")
except Exception as e:
    print(f"❌ Error processing BLAST results: {e}")
    sys.exit(1)

# === 5. Retrieve and save the sequences of the hits
print(f"\n📋 STEP 4: Retrieving BLAST hit sequences")

# Extract unique accession list
try:
    accessions = df["Accession"].unique().tolist()
    total_accessions = len(accessions)
    print(f"🔍 Found {total_accessions} unique accessions to retrieve")
    
    # Function to split list into chunks for batch processing
    def chunks(lst, n):
        """Split a list into chunks of size n"""
        for i in range(0, len(lst), n):
            yield lst[i:i + n]

    batch_size = 10
    sequences_retrieved = 0
    sequences_failed = 0
    
    # Create progress bar for the total process
    with tqdm(total=total_accessions, desc="Total progress", unit="seq") as main_pbar:
        # Save the sequences in FASTA format
        with open(blast_sequences, "w") as output_handle:
            # Process in batches to avoid overloading the NCBI API
            for i, batch in enumerate(chunks(accessions, batch_size)):
                batch_str = ",".join(batch)
                try:
                    # Show information about the current batch
                    tqdm.write(f"📦 Batch {i+1}/{-(-total_accessions//batch_size)}: Fetching {len(batch)} sequences...")
                    
                    # Fetch the sequences in batch
                    handle = Entrez.efetch(db="nucleotide", id=batch_str, rettype="fasta", retmode="text")
                    records = list(SeqIO.parse(handle, "fasta"))
                    handle.close()
                    
                    # Save sequences
                    for record in records:
                        SeqIO.write(record, output_handle, "fasta")
                    
                    # Update counters
                    sequences_retrieved += len(records)
                    main_pbar.update(len(batch))
                    
                    if len(records) < len(batch):
                        missing = len(batch) - len(records)
                        sequences_failed += missing
                        tqdm.write(f"⚠️  Warning: {missing} sequences were not found in this batch")
                    
                    # Respect the usage limits of the NCBI API
                    time.sleep(1)  # Wait 1 second between requests
                except Exception as e:
                    tqdm.write(f"⚠️ Error retrieving sequences: {e}")
                    sequences_failed += len(batch)
                    main_pbar.update(len(batch))
                    
                    # Retry after waiting a bit longer
                    tqdm.write("🔄 Retrying after pause...")
                    time.sleep(5)
                    try:
                        handle = Entrez.efetch(db="nucleotide", id=batch_str, rettype="fasta", retmode="text")
                        records = list(SeqIO.parse(handle, "fasta"))
                        handle.close()
                        
                        # Save sequences retrieved on the second attempt
                        for record in records:
                            SeqIO.write(record, output_handle, "fasta")
                        
                        # Update counters (correct the failure counter)
                        sequences_retrieved += len(records)
                        sequences_failed -= len(records)
                    except Exception as e:
                        tqdm.write(f"❌ Permanent failure retrieving batch: {e}")
    
    # Final report
    success_rate = (sequences_retrieved / total_accessions) * 100 if total_accessions > 0 else 0
    print(f"\n✅ BLAST hit sequences saved in {blast_sequences}")
    print(f"📊 Retrieval statistics:")
    print(f"   - Total sequences requested: {total_accessions}")
    print(f"   - Successfully retrieved sequences: {sequences_retrieved} ({success_rate:.1f}%)")
    if sequences_failed > 0:
        print(f"   - Sequences not found: {sequences_failed}")

except Exception as e:
    print(f"❌ Error retrieving sequences: {e}")

print("\n" + "=" * 70)
print("🧬 BLAST ANALYSIS PIPELINE - COMPLETED")
print("=" * 70)
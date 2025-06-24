# 🔬 Multi BLAST Pipeline

A Python script for batch BLASTn searches against NCBI's nucleotide (nt) database using remote access.  
It processes multiple FASTA files, splits large files into manageable chunks, runs BLAST remotely, parses the results, retrieves top-hit sequences, and consolidates all outputs into merged files.

## 🔧 Features

- Supports multiple input FASTA files
- Automatically splits files with many sequences to avoid NCBI IP blocking
- Cleans sequences (removes gaps)
- Runs BLASTn remotely using NCBI's servers
- Parses and exports BLAST hits to Excel and TSV
- Retrieves top-hit sequences from GenBank (FASTA)
- Merges all results into unified summary files
- Warns user when input load may trigger rate limits

## ⚠️ Important Notes

- This script uses **remote BLAST**, which is slower and limited by NCBI's servers.
- If too many sequences are submitted at once, **your IP may be temporarily blocked**. The script splits files and adds delays to reduce this risk.
- For large-scale analyses, it is **strongly recommended** to set up a local BLAST database (`nt`) and run BLAST locally.

## 📂 Folder Structure

```
project/
│
├── files-fasta/         # Input FASTA files
│   ├── example1.fasta
│   └── ...
│
├── outputs/             # Final merged results (TSV, Excel, FASTA)
├── temp/                # Temporary intermediate files
├── blast_pipeline_remote.py
└── requirements.txt
```

## 📦 Installation

Make sure you have:

- Python 3.7+
- NCBI BLAST+ installed (`blastn` available in your PATH)
- Internet connection

Then install the required Python packages:

```bash
pip install -r requirements.txt
```

## 🧬 How to Use

1. Edit the script and configure your NCBI email:
   ```python
   Entrez.email = "your_email@example.com"
   ```
2. Place your `.fasta` files inside the `files-fasta/` folder.
3. Run the pipeline:
   ```bash
   python blast_pipeline_remote.py
   ```

## 📤 Outputs

After running, you’ll get:

- `example_combined_blast.xlsx` – Clean Excel file of best BLAST hits
- `example_combined_hits.fasta` – FASTA file of matched reference sequences
- `example_combined_blast.tsv` – Raw BLAST tabular output

## 🙌 Contributions

Pull requests are welcome!  
If you find bugs or want to request features, feel free to open an issue.

## ✍️ Author

Developed by **Romualdo JP**  
📧 Contact: jpromualdo0@gmail.com

## 📄 License

This project is licensed under the MIT License – see the [LICENSE](LICENSE) file for details.

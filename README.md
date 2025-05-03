# 🔬 BLAST Sequences Pipeline

This repository contains a complete and automated Python pipeline to:

- Clean chosen sequences from gaps
- Run remote **BLASTn** searches against the **NCBI nt database**
- Parse and organize the top BLAST hits
- Export results to **Excel**
- Retrieve the matched sequences from NCBI using **Entrez**

---

## 🚀 Features

- Cleans and prepares the chosen sequences
- Performs **remote BLAST** against NCBI (`nt`, `refseq_rna`, etc.)
- Extracts best hits and generates a clean Excel table
- Retrieves FASTA sequences of matched hits using `Entrez`
- Progress bars, logging, retry logic, and failure reporting

---

## 📦 Installation

Make sure you have:

- Python 3.7+
- NCBI BLAST+ installed (`blastn` available in your PATH)
- Internet connection

Then install the required Python packages:

```bash
pip install -r requirements.txt
```

---

## 🧬 How to Use

1. Edit the script and configure your NCBI email:
```python
Entrez.email = "your_email@example.com"
```

2. Optionally, configure your NCBI API key:
```python
# Entrez.api_key = "your_api_key"
```

3. Run the pipeline:
```bash
python blast_pipeline.py
```

Make sure your input file (e.g., `example.fasta`) is located in the same folder or properly referenced.

---

## 📤 Outputs

After running, you’ll get:

- `blast_hits_remote.xlsx` – Clean Excel file of best BLAST hits
- `blast_hit_sequences.fasta` – FASTA file of matched reference sequences
- `blast_remote.tsv` – Raw BLAST tabular output
- `Sequence_Blasted_clean.fasta` – Cleaned Chosen Sequence input (no gaps)

---

## 🧪 Example Data

A sample input file is provided in the [`example/`](example/) folder:
```bash
python blast_pipeline.py
```

---

## 🙌 Contributions

Pull requests are welcome! If you find bugs or want to request features, feel free to open an issue.

---

## ✍️ Author

Developed by **Romualdo JP**  
Contact: jpromualdo0@gmail.com

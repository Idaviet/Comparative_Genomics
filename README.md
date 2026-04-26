# Comparative Genomics Toolkit

A Python toolkit for comparative genomics analysis, including sequence analysis, distance scoring, and ORF prediction.

## Features

- **Sequence Analysis**: Calculate GC content, dinucleotide frequencies, and di-amino acid frequencies
- **Distance Scoring**: Generate distance matrices and visualizations for phylogenetic analysis
- **ORF Prediction**: Identify open reading frames in DNA sequences
- **Belvu Matrix Export**: Export distance matrices in Belvu-compatible format

## Project Structure
Comparative-Genomics/
├── src/                          # Source code
│   ├── sequence_analysis.py      # DNA sequence analysis functions
│   ├── distance_scoring.py       # Distance matrix generation
│   └── orf_prediction.py         # ORF prediction algorithms
├── data/                         # Example/test data
├── results/                      # Output directory (created at runtime)
├── requirements.txt              # Python dependencies
└── README.md

## Installation

```bash
git clone https://github.com/Idaviet/Comparative-Genomics.git
cd Comparative-Genomics
pip install -r requirements.txt

python src/main.py /path/to/fasta/folder
Follow the on-screen prompts to:
  1. Select analysis type (Distance Scoring, ORF Prediction, or Both)
  2. Choose files to analyze
  3. Configure output options

# Comparative Genomics Toolkit

A Python toolkit for comparative genomics analysis, including sequence characterization, distance-based phylogenetic inference, and open reading frame (ORF) prediction.

## Features

- **Sequence Analysis**: Calculate GC content, dinucleotide frequencies, and di-amino acid frequencies across all three reading frames
- **Distance Scoring**: Generate pairwise distance matrices using Euclidean distance metrics for:
  - GC content
  - Dinucleotide composition
  - Di-amino acid composition (per reading frame)
- **ORF Prediction**: Identify open reading frames in prokaryotic genomes with customizable minimum length thresholds
- **ORF Validation**: Evaluate ORF prediction performance against reference proteomes using pairwise alignment
- **Belvu Matrix Export**: Export distance matrices in Belvu-compatible format for phylogenetic tree construction

## Installation

```bash
git clone https://github.com/Idaviet/Comparative_Genomics.git
cd Comparative_Genomics
pip install -r requirements.txt
```

## Usage
### Command Line Interface
```bash
python src/cli.py /path/to/genomes/folder
```

The CLI provides an interactive workflow to:
* Select files for analysis
* Choose analysis type (distance scoring, ORF prediction, or both)
* Configure parameters (minimum ORF length, distance metrics)
* Generate output files and visualizations

### ORF Performance Evaluation
Compare predicted ORFs against a reference proteome:

```bash
python src/orf_evaluation.py reference_proteome.fasta genome.fasta 300 90
```

#### Arguments:
* **reference_proteome.fasta**: Known protein sequences
* **genome.fasta**: Genome to analyze
* **300**: Minimum ORF length (bp)
* **90**: % identity threshold for true positive classification

## Output
### Distance Analysis
* **distance_parsing/** - Output directory
  * **sequence_analysis/** - GC/dinuc/diaa frequency tables and plots
  * **belvu_matrices/** - Belvu-format distance matrices (.fasta)

### ORF Prediction
* **Sequence_ORFs/** - Output directory  
  * **Individual ORF FASTA files per genome**
  * **ORF_summary.txt** - ORF counts per strand

## Project Structure
```plain text
Comparative-Genomics/
├── src/
│   ├── sequence_analysis.py    # DNA sequence utilities and calculations
│   ├── distance_matrix.py      # Distance scoring and matrix generation
│   ├── orf_prediction.py       # ORF finding algorithms
│   ├── orf_evaluation.py       # ORF validation against references
│   └── cli.py                  # Command-line interface
├── data/                       # Example/test data
├── results/                    # Output directory (created at runtime)
├── requirements.txt
└── README.md
```

## Algorithm Details
### Distance Metrics
All distances use Euclidean distance:

* GC distance: √((GC₁ - GC₂)²)
* Dinucleotide distance: √(Σ(freq₁ᵢ - freq₂ᵢ)²) for all 16 dinucleotides
* Di-amino acid distance: Same formula applied to 400 possible di-aa pairs

### ORF Prediction
Identifies ORFs as regions between start codons (ATG, GTG, TTG) and stop codons (TAA, TAG, TGA) on both strands.

## Author
Isaac Daviet — [GitHub](https://github.com/Idaviet)

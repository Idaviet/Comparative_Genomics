from src.sequence_analysis import parse_fasta, gc_calc, dinuc_calc, reading_frames
from src.distance_matrix import all_distances, belvu_matrix
from src.orf_prediction import extract_orfs

# Parse FASTA file
sequences = parse_fasta("genome.fasta")

# Calculate GC content
gc = gc_calc(sequences[">NC_12345"])

# Generate all distance matrices
distances = all_distances(analyzed_entries)

# Extract ORFs
orfs = extract_orfs("genome.fasta", min_length=300)
"""
Parsing and analyis functions for exercise 3 of the final project. It uses the
parse_fasta, combine_all_fasta and combine_selected_fasta functions imported 
from the 'FinalProject_group2_ex1.py' file to parse a selected set of fasta 
files from a submitted folder. Functions within thi file are then used to 
predict ORFs for both the forward and reverse strands of a given sequence.

Script contains an argparse function to call belvu_dist_parsing function from
a terminal.

INPUT FOLLOWING COMMAND ON TERMINAL AND FOLLOW ON SCREEN INSTRUCTIONS
    >python FinalProject_group2_ex3.py <folder containing fasta files to be analyzed>

FUNCTIONS:
    
    revcomp_parse_fasta: Uses the parse_fasta function from ex1 to parse all
                         sequences in a given fasta file as well as their 
                         reverse complements. Reverse complement sequences have
                         the same name as the parent sequence + '_ReverseComplement'
                         Returns a dictionary in the following format:
                         seqs{entry:(sequence as list)} (for both the forward 
                                                         and reverse sequences)
                         
    extract_orfs: Uses the revcomp_parse_fasta function to parse a given fasta
                  file and predict ORFs from the resulting dictionary. ORFs are
                  identified by their start and end positions as follows:
                  ORF_XXXXX-YYYYY.
                  Returns a ditionary of the following format:
                  all_orfs{<sequence_name>:['ORF#1', 'ORF#2',... 'ORF#n']}
                  
    write_orf_fast: Uses the dictionary returned from the previous function
                    to save the predicted ORFs into a fasta file, one per 
                    sequence parsed 
                    NOTE: this includes the reverse_complement sequences, which
                          are saved in a file separate from the parent.
                    NOTE: Names for the files are simplified and are based on
                          the assumption that the sequences are named following
                          the NCBI fasta name format
                          
    ORF_summary: Uses the dictionary returned from the extract_orfs function
                 to calculate the number of predicted ORFs on the forward and 
                 reverse strand. 
                 Returns a dictinary in the following format:
                 summary{<sequence_name>:[<# forward ORFs>, <# reverse ORFs>]}
    
    full_orf_parsing: Uses all functions previously defined to predict ORFs of
                      selected fasta files within a given folder. Returns the
                      predicted orf fasta files as well as a summary txt file
                      containing the forward, reverse and total # of ORFs for
                      all sequences.
                      THIS IS THE FUNCTION CALLED BY THE ARGPARSE MODULE

"""

import argparse
import os
import sys
file_location = r"C:\Users\isaac\Desktop\SU Semester 2023\Comparative Genomics\Final Project"
sys.path.append(file_location)
import sequence_analysis_2 as ex1
from Bio.Seq import Seq

def revcomp_parse_fasta(file):
    """
    Extracts reverse complement of sequences from dictionary resulting from 
    parse_fasta function and adds them to the dictionary resulting from said 
    function
    
    INPUT: fasta file containing all sequences to be reverse complemented
    OUTPUT: seqs{entry:(sequence as list)} (both the forward and reverse
                                           sequences)
    """
    
    seqs = ex1.parse_fasta(file)
    all_seqs = {}
    for entry in seqs:
        seq = seqs[entry]
        name = entry
        rc_name =  entry + '_ReverseComplement'
        
        rc_seq = ''.join(seq)
        rc_seq = Seq(rc_seq)
        rc_seq = rc_seq.reverse_complement()
        all_seqs[name] = seq
        all_seqs[rc_name] = list(rc_seq)
    
    return all_seqs



def extract_orfs(fasta, minimum_codon_length):
    """
    Parses sequences in fasta file and extracts potential ORFs based on the 
    'ATG' starting codon and the 'TAA', 'TAG', and 'TGA' stop codons. ORFs are 
    saved as a single string in a dictionary of dictionaries. IT assumes a 
    minimum length of 150 nucleotides.ORFs are named as follows: 
        <sequence_name>-ORF_<starting_position>-<end position>
    
    INPUT: Fasta file with sequence to analyze
    OUTPUT: all_orfs{<species_name>:['ORF#1', 'ORF#2',... 'ORF#n']}
    
    """
    fasta = revcomp_parse_fasta(fasta)
    start_codons=['ATG','atg'] ; stop_codons=['TAA','TAG','TGA','taa','tag','tga'] #reference list of start and stop codons
    
    all_orfs = {}
    
    for entry in fasta:
        seq = fasta[entry]
        name = entry
        start_read = 0
    
        orfs = {}
        for i in range(0, len(seq)-3):
            if i > start_read:
                codon1 = seq[i] + seq[i+1] + seq[i+2]
                orf_name = 'ORF_'
                orf = ''
            
                if codon1 in start_codons:
                    start_pos = i
                    start = str(i) + '-'
                    orf_name += start
                                    
                    for j in range (i+3, len(seq)-3, 3): 
                        codon2 = seq[j] + seq[j+1] + seq[j+2]
                        length = j+2 - i
                        
                        if length < minimum_codon_length and codon2 in stop_codons:
                            start_read = i+1
                            break
                        
                        elif codon2 in stop_codons:
                            end_pos = j+3
                            end = str(j+3) 
                            orf_name += end
                            start_read = j+3
                            orf = seq[start_pos:end_pos]
                            break
            
                if len(orf) != 0:
                    orf = ''.join(orf)
                    length = len(list(orf))

                    if length %3 == 0 and length >= minimum_codon_length:
                        orfs[orf_name] = orf
        print('\nORF prediction of', entry, 'is complete')
                
        all_orfs[name] = orfs
    print('\n\nORFS FROM ALL FASTA FILES HAVE BEEN PREDICTED\n\n')
    return all_orfs
 


def write_orf_fasta(all_orfs, outfolder):
    """
    Saves ORFs from all_orfs dictionary returned by extract_orfs function into 
    fasta file. Each individual sequence in the all_orfs dictionary will be saved
    in a separate fasta file. The function assumes the sequences are named 
    in the NCBI genome format and will parse the sequence names accordingly when
    naming the files.
    
    INPUT: all_orfs dictionary from extract_orfs function
           folder where the fasta files are to be saved
           
    OUTPUT: one fasta file for each sequence contained in the all_dist dictionary
            with the ORFs and their locations contained therein
    """
    
    outfolder = outfolder.replace('\\', '/') 
    for species in all_orfs:
        orfs = all_orfs[species]
        name = species.split(' ')
        name = name[1] + '_' + name[2]
        if '_ReverseComplement' in species:
            name += '_ReverseComplement'
        file = outfolder + '/ORFs_' + name + '.fasta'
        
        with open(file, 'w') as f:
            for orf in orfs:
                entry = '>' + name + '-' + orf + '\n'
                sequence = orfs[orf] + '\n\n'
                all = entry + sequence
                f.write(all)

                
def ORF_summary(all_orfs):
    """
    Produces a summary dictionary with the total number of predicted ORFs for 
    each sequence contained in the all_orfs dictionary from the extract_orfs
    function.
    
    INPUT: all_orfs dictionary
    OUTPUT: summary{<sequence_name>:[<# forward ORFs>, <# reverse ORFs>]
    """
    summary = {}
    for species1 in all_orfs:
        for species2 in all_orfs:
            if not "_ReverseComplement" in species1:
                if "_ReverseComplement" in species2: 
                    id1 = species1
                    id2 = species2.replace('_ReverseComplement', '')
                    if id1 == id2:
                        orfs1 = all_orfs[species1]
                        orfs2 = all_orfs[species2]
                        count = len(orfs1)
                        rc_count = len(orfs2)
            
                        summary[species1] = [count, rc_count]
    return summary

def full_orf_parsing(folder):
    """
    Full parsing of genomes held within a folder. Can choose to parse all files
    contained in said folder, or to choose a custom selection of files. Results
    in a new folder name 'Sequence_ORFs' containing the ORF fasta files of all
    selected files. Will also output a summary file containing the number of 
    ORFs predicted for each sequence and a compiled fasta file containing all
    the selected sequences.
    
    INPUT: folder containing fasta files to be analyzed
    OUTPUT: folder with ORF predictions for all sequences selected in separate
            fasta files
            Summary file containing number of predicted ORFs for each sequence
            Combined fasta file containing all the selected fasta files
    """
    full_analysis = input('Would you like to parse all fasta files inside this folder?\n\t[y/n]\n\n')
    folder = folder.replace('\\', '/')
 

    if full_analysis == 'n':
        ex1.combine_selected_fasta(folder)
        file = folder + '/all_selected_sequences.fasta'
    
    elif full_analysis == 'y':
        ex1.combine_all_fasta(folder)
        file = folder + '/all_sequences.fasta'
    
    outfolder = folder + '/Sequence_ORFs'
    os.mkdir(outfolder)
   
    all_orfs = extract_orfs(file)
    
    summary = ORF_summary(all_orfs)
    write_orf_fasta(all_orfs, outfolder)
    sum_file = outfolder + '/ORF_summary.txt'
    
    with open(sum_file, 'w') as f:
        for species in summary:
            forward = summary[species][0]
            reverse = summary[species][1]
            total = forward + reverse
            entry = species + ':\n\t' + str(forward) + '\t\tORFs in forward strand\n\t' + str(reverse) +'\t\tORFs in reverse strand\n\t' +str(total) +'\tORFs total\n\n' 
            f.write(entry)
            print(entry)
    print('ORFS FOR ALL SELECTED FILES HAVE BEEN COMPILED')
    return all_orfs       
 

                
if __name__ == "__main__":
    """
    Argument parser function to generate potential ORFs for all sequences within a folder
    from computer terminal.
    
    INPUT FORMAT:
    >python ex3.py <input folder containing fasta files>
    """
    parser = argparse.ArgumentParser(description="Generates potential ORFs in a\
                                     fasta file format from selected fasta files\
                                         within a given folder", \
                            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('folder', type=str, 
                        help='Please input folder containing fasta files to be analyzed')
                    
    args = parser.parse_args()
    full_orf_parsing(args.folder)          
        
            
    



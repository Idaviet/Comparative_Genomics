"""
Parsing and analyis functions for exercise 2 of the final project. It uses the
combine_all or selected_fasta functions imported from the 
'FinalProject_group2_ex1.py' file to parse a selected set of fasta files from
a submitted folder. The full_dna_analysis function from the same file to 
generate AnalyzedSeq objects of the sequence, which are then used as the input
for the belvu_dist_parsing function to generate distance scores, bar plots of 
said scores for quality control purposes and, if desired, belvu distance 
matrices using %GC, dinucleotide frequencies, di-amino acid frequencies (either)
reading frame #1, 2 or 3) or all of the above.

Script contains an argparse function to call belvu_dist_parsing function from
a terminal.

INPUT FOLLOWING COMMAND ON TERMINAL AND FOLLOW ON SCREEN INSTRUCTIONS
    >python FinalProject_group2_ex2.py <folder containing fasta files to be analyzed>

FUNCTIONS:

    gc_distances: Function to calculate distances based on %GC content from an 
                  list of AnalyzedSeq objects. 
                  Returns a dictionary in the following format:
                  gc_dist{"<sequence1>.id/<sequence2>.id":float(%GC distance)}
                  
    dinuc_distances: Function to calculate distances based on dinucleotide 
                     frequencies from an list of AnalyzedSeq objects. 
                     Returns a dictionary in the following format: 
                     dinuc_dist{"<sequence1>.id/<sequence2>.id" : float(DN distance)}
    
    diaa_distances: Function to calculate distances based on diamino acid 
                    frequencies of A SINGLE reading frame from an list of 
                    AnalyzedSeq objects. 
                    Returns a dictionary in the following format:
                    diaa_dist{"<sequence1>.id/<sequence2>.id" : float(DA distance)}
                    
    all_distances: Uses previous three functions to calculate all possible 
                   distance scores from an list of AnalyzedSeq objects including 
                   %GC, dinuc and all 3 diaa frequencies.
                   Returns a dictionary in the followig format:
                   all_distances{'%GC': gc_dist{}, 
                                 'dinuc': dinuc_dist{}, 
                                 'diaa_rf1' : diaa_dist{}, 
                                 'diaa_rf2' : diaa_dist{}
                                 'diaa_rf3' : diaa_dist{}}
                   
    plot_distances: Uses all_distances function to plot resulting distances in
                    a set of bar graphs and saves them in a given folder. 
                    It also returns the same dictionary as the all_distances
                    function, and so can be used in its stead in any program.
                    
    extract_dist_list: Extracts the information from the all_distances 
                       dictionary and reformats them into a list to be used when
                       writing the distace matrix file
                       Returns a list in the following format:
                       dist_list[[seq1, distance, seq2], [seqA, distance, seqB], ...]
                       
    belvu_matrix: uses the list produced by the previous function to write a 
                  fasta file containing the belvu distance matrix that can be 
                  used in phylogenetic tree construction.
    
    belvu_dist_parsing: Combines all functions listed above to generate a 
                        combined fasta file, a sequence analysis score txt file
                        bar plots of the distance scores and a belvu distance
                        matrix using one or all of the distance scores.
                        THIS IS THE FUNCTION CALLED BY THE ARGPARSE MODULE
"""

import argparse
import os
import sys
import math
import matplotlib.pyplot as plt
file_location = r"C:\Users\isaac\Desktop\SU Semester 2023\Comparative Genomics\Final Project"
sys.path.append(file_location)
import sequence_analysis_2 as ex1
from sequence_analysis_2 import full_dna_analysis as fda





def gc_distances(entries):
    """
    Creates a dictionary with the %GC distances between each pair of sequences 
    parsed through the full_dna_analysis function of exercise 1. Requires all
    of the sequences to be of 'AnalyzedSeq' object class. Returns a dictionary
    containing the sequence id as key and distance as value.
    
    %GC distanc = sqrt( (%GC_seq1 - %GC_seq2)^2 )
    
    INPUT: list of 'AnalyzedSeq' type sequences
    OUTPUT: gc_dist{"<sequence1>.id/<sequence2>.id" : float(%GC distance)}
    """
    gc = {}
    for entry in entries:
        gc[entry.id] = entry.gc
    
    gc_dist = {}
    for entry1, gc1 in gc.items():
        for entry2, gc2 in gc.items():
            distance = math.sqrt((gc1-gc2)**2)
            name = entry1 + '//AGAINST//' + entry2
            gc_dist[name] = distance
    return gc_dist


def dinuc_distances(entries):
    """
    Creates a dictionary with the dinucleotide distances between each pair of 
    sequences parsed through the full_dna_analysis function of exercise 1. 
    Requires all of the sequences to be of 'AnalyzedSeq' object class. 
    Returns a dictionary containing the sequence id as key and distance as 
    value.
    
    dinuc distance = 
    sqrt( (DN1_freq_seq1 - DN1_seq2)^2 + (DN2_freq_seq1 - DN2_seq2)^2 + 
         (DNn_freq_seq1 - DNn_seq2)^2)
    
    Where DNn = nth dinculeotide found in sequences
    
    INPUT: list of 'AnalyzedSeq' type sequences
    OUTPUT: dinuc_dist{"<sequence1>.id/<sequence2>.id" : float(DN distance)}
    """
    dinuc = {}
    for entry in entries:
        dinuc[entry.id] = entry.dinuc
    
    dinuc_dist = {}
    
    for entry1, dinucs1 in dinuc.items():
        for entry2, dinucs2 in dinuc.items():
            
            calculation = []
            name = entry1 + '//AGAINST//' + entry2
            
            for dinuc1, freq1 in dinucs1.items():
                for dinuc2, freq2 in dinucs2.items():
                    if dinuc1 == dinuc2:
                        distance = (freq1-freq2)**2
                        calculation.append(distance)
                        
            total_dist = math.sqrt(sum(calculation))
            dinuc_dist[name] = total_dist
    return dinuc_dist

def diaa_distances(entries, rf):
    """
    Creates a dictionary with the di-amino distances between each pair of 
    sequences parsed through the full_dna_analysis function of exercise 1. 
    Requires all of the sequences to be of 'AnalyzedSeq' object class. 
    Returns a dictionary containing the sequence id as key and distance as 
    value. Requires a secondary input specifying which reading frame to use.
    NOTE: only the specified reading frame witth be parsed.
    
    dinuc distance = 
    sqrt( (DA1_freq_seq1 - DA1_seq2)^2 + (DA2_freq_seq1 - DA2_seq2)^2 + 
         (DAn_freq_seq1 - DAn_seq2)^2)
    
    Where DNn = nth di-amino acid found in sequences
    
    INPUT: list of 'AnalyzedSeq' type sequences
           reading frame, from 1 to 3
    
    OUTPUT: diaa_dist{"<sequence1>.id/<sequence2>.id" : float(DA distance)}
    """
   
    diaa_hold = {}
    diaa_rf = {}
    rf = 'rf'+str(rf)
    
    for entry in entries:
        diaa_hold[entry.id] = entry.diaa
        
    for entry in diaa_hold:
        spec_diaas={}
        for diaa in diaa_hold[entry]:
            if diaa[0] == rf:
                spec_diaas[diaa[1]] = diaa[2]
            diaa_rf[entry] = spec_diaas
            
    diaa_dist = {}
    
    for entry1, diaas1 in diaa_rf.items():
        for entry2, diaas2 in diaa_rf.items():
            
            calculation = []
            name = entry1 + '//AGAINST//' + entry2
            
            for diaa1, freq1 in diaas1.items():
                for diaa2, freq2 in diaas2.items():
                    if diaa1 == diaa2:
                        distance = (freq1-freq2)**2
                        calculation.append(distance)
                        
            total_dist = math.sqrt(sum(calculation))
            diaa_dist[name] = total_dist
    return diaa_dist
    
def all_distances(entries):
    """
    Compiles all distances (%GC, dinuc, di-aa rf1/2/3) calculated from the
    gc/dinuc/diaa_distances functions directly into a single dictionary.
    
    INPUT: list of 'AnalyzedSeq' type sequences
    OUTPUT: all_distances{'%GC': gc_dist{}, 'dinuc': dinuc_dist{}, 
                          'diaa_rf1' : diaa_dist{}, 'diaa_rf2' : diaa_dist{}
                          'diaa_rf3' : diaa_dist{}}
    """
    all_distances = {}
    all_distances['%GC'] = gc_distances(entries)
    all_distances['dinuc'] = dinuc_distances(entries)
    all_distances['diaa_rf1']= diaa_distances(entries,1)
    all_distances['diaa_rf2'] = diaa_distances(entries,2)
    all_distances['diaa_rf3'] = diaa_distances(entries,3)
    return all_distances


def plot_distances(entries, outfolder):
    """
    Plots all distances dictionaries as bar graphs using pyplot as a means of
    quality control and analysis. Each distance type of is plotted separately.
    
    INPUT: list of 'AnalyzedSeq' type sequences
    OUTPUT: one plot for each type of distance
            returns all_distances dictionary
    """
    all_dist = all_distances(entries)
    
    for type in all_dist:
        x = []
        y = []
        for pair, dist in all_dist[type].items():
            x.append(pair)
            y.append(dist)
        
        plt.figure(figsize=(20,10))
        plt.bar(x, y, color='blue')
        plt.xticks(rotation = 90)
        plt.xlabel("PAIRS")
        plt.ylabel('DISTANCES') 
        plt.tight_layout()
        file = outfolder + "/" + type + "_dist_plot" + '.png'
        plt.savefig(file)

    return all_dist

def extract_dist_list(all_dist, type):
    """
    Extracts a distance list for all species included in the all_dist dictionary
    returned from the all_distances or plot_distances functions. Returns a distances
    list to be used in the belvu_matrix function.
    
    INPUT: all_dist dictionary returned from all_distances function
           which distances to use, input as a string (%GC, dinuc, diaa_rf1/2/3)
    
    OUTPUT: list containing all pairs and their distances
            dist_list[[seq1, distance, seq2], [seqA, distance, seqB], ...]
    """
    dist_list = []
    sample1 = []
    sample2 = []
    values = []
    for types in all_dist:
        if types == type:
            dict = all_dist[type]
            for pair in dict:
                pairs = pair.split('//AGAINST//') ###look to changing the '/' condition
                sample1 = pairs[0]
                sample2 = pairs[1]
                values = dict[pair]
                data = [sample1, values, sample2]
                dist_list.append(data)
    return(dist_list)

def belvu_matrix(dist_list, output):
    """
    Uses dist_list returned from extract_dist_list function to write a distance
    matrix in the appropriate format for phylogenetic tree construction using 
    belvu. Requires an output folder to save the resulting fasta file.
    
    INPUT: dist_list for extract_dist_list function
           folder where the resulting fasta is to be saved
           
    OUTPUT: fasta file containing distance matrix in specified folder
    """
    
    with open(output, 'w') as f:
        all_species= []
        n = 0
        for pair in dist_list:
            if pair[0] not in all_species:
                all_species.append(pair[0])
                all_species.append('\t')
                n+=1
      
        entry = ''      
        for i in range(0, len(all_species)-1):
            entry += str(all_species[i])
        f.write(entry)
        
        for i in range(0,len(dist_list)):
            pair = dist_list[i]
            value = round(pair[1], 4)
            values = str(value)
            if i == 0:
                f.write('\n')
                f.write(values)
                f.write('\t')
                
            elif i % n == 0:
                f.write('\n')
                f.write(values)
                f.write('\t')
                
            else:
                f.write(values)
                f.write('\t')
                

def belvu_dist_parsing(folder):
    """
    Full parsing of all sequences held within input folder. Results in the a
    full sequence analysis file and the bar plot graphs for all possible 
    distance scores determined from the previous functions. Includes an option 
    to extract a belvu matrix for one type of score, or all.
    
    INPUT: folder containing sequences wished to be analyzed
    OUTPUT: bar plots containing scores for all 5 possible score types
            Sequence analysis txt file containing %GC, dinucleotide and di-amino
            frequecies for all sequences analyzed
            OPTIONAL: belvu distance matrix for selected or all types.
    """

    folder = folder.replace('\\', '/')   
    
    full_analysis = input('Parse belvu distance matrix, or hold until seq analysis is completed? \n\t[y/n/hold]\n')
    type = []
    if full_analysis == 'y':
        inputs = input('Please enter which distance scores to use:\n\t%GC content [gc]\n\tDinulceotide Frequency [dn]\n\tDiamino Frequency [da]\n\tAll matrices [all]\n')  
        if inputs == 'gc':
            type.append('%GC')
        elif inputs == 'dn':
            type.append('dinuc')
        elif inputs == 'da':
            rf = input('Which reading frame? [1/2/3] \n')
            diaa = 'diaa_rf' + rf
            type.append(diaa)
        elif inputs == 'all':
            type = ['%GC','dinuc', 'diaa_rf1', 'diaa_rf2', 'diaa_rf3']
    
    simplify = input('Would you like to simplify NCBI genome names? \n\t[y/n]\n')
    ex1.combine_fasta(folder)
    
    out_folder = folder + '/belvu_distance_parsing'
    os.mkdir(out_folder)
    
    fasta_file = folder + '/all_sequences.fasta'
   
    entries = fda(fasta_file)
    
    if simplify == 'y':
        entries = ex1.ncbi_name_simplifier(entries)
    
    all_dist = plot_distances(entries, out_folder)
    
    if full_analysis == 'hold':
        q = input('Would you like to extract a belvu distance matrix?\n\t[y/n]\n')
        if q == 'y':
            inputs = input('Please enter which distance scores to use:\n\t%GC content [gc]\n\tDinulceotide Frequency [dn]\n\tDiamino Frequency [da]\n\tAll matrices [all]\n')  
            if inputs == 'gc':
                type.append('%GC')
            elif inputs == 'dn':
                type.append('dinuc')
            elif inputs == 'da':
                rf = input('Which reading frame? [1/2/3] \n')
                diaa = 'diaa_rf' + rf
                type.append(diaa)
            elif inputs == 'all':
                type = ['%GC','dinuc', 'diaa_rf1', 'diaa_rf2', 'diaa_rf3']
            
    if full_analysis == 'y' or 'hold':
        for i in type:
            output = out_folder + '/' + i + '_belvu_distance_matrix.fasta'
            dist_list = extract_dist_list(all_dist, i)
            belvu_matrix(dist_list, output)
    
     
                
if __name__ == "__main__":
    """
    Argument parser function to generate full belvu distance matrix from computer terminal.
    
    INPUT FORMAT:
    >python ex2.py <input folder containing fasta files>
    """
    parser = argparse.ArgumentParser(description="Generates bar plot graphs of\
                                     all distance types for all fasta files\
                                         within a folder, with option to generate\
                                             belvu distance matrix", \
                            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('folder', type=str, 
                        help='Please input folder containing fasta files to be analyzed')
                    
    args = parser.parse_args()
    belvu_dist_parsing(args.folder)            
        

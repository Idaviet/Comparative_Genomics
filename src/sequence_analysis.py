# -*- coding: utf-8 -*-
"""

"""

import matplotlib.pyplot as plt
import argparse
import os
import sys
import math
import matplotlib.pyplot as plt
from Bio.Seq import Seq

def ncbi_name_simplifier(entries, y_n_q):
    """
    Parses AnalyzedSeq class objects from full_dna_analysis function, 
    simplifying their NCBI names to only contain the species name.
    
    INPUT: entries returned from full_dna_analysis function
    OUTPUT: same entries, with their entry.id value changed to a only the
            include species name
    """
    for entry in entries:
        name = entry.id
        name = name.split(' ')
        name = name[1]+'_'+name[2]
        entry.id = name
    return entries

def combine_all_fasta(file_path):
    """
    Combines all fasta files in a folder into a single file named 
    'all_sequences.fasta' located in the submitted folder
    
    INPUT: filepath of folder containing files to be analyzed
    OUTPUT: all_sequences.fasta file containg all sequences in all fasta files
    
    """
    file_path = file_path.replace('\\', '/')
        
    all_genomes = file_path + '/all_sequences.fasta'
    
    with open(all_genomes, 'w') as f:
        for i in os.listdir(file_path):
            if '.fasta' in i:
                with open(os.path.join(file_path, i), 'r') as fasta:
                    text = fasta.read()
                    f.write(text)

def combine_selected_fasta(file_path):
    """
    Combines select fasta files in a folder into a single file named 
    'all_selected_sequences.fasta' located in the submitted folder. Once function
    is called with target folder, user will be asked to submit the files to be
    processed, separated by a space. 
    Eg: <sequence1.fasta sequence2.fasta... sequenceN.fasta>
    
    INPUT: filepath of folder containing files to be analyzed
           
    OUTPUT: all_sequences.fasta file containg all sequences in all fasta files
    """
    file_path = file_path.replace('\\', '/')
    
    selected_genomes = file_path + '/all_selected_sequences.fasta'
    selected_sequences = input('Please input files names to be analyzed (with fasta extension and separated by a space):\n\n')
    selected_sequences = selected_sequences.split(' ')
    seq_files = []
    
    for seq in selected_sequences:
        s = file_path + '/' + seq
        seq_files.append(s)
    
    with open(selected_genomes, 'w') as f:
        for file in seq_files:
            with open(file, 'r') as fasta:
                text = fasta.read()
                f.write(text)

class AnalyzedSeq:
    def __init__(self, id, gc, dinuc, diaa):
        self.id = id
        self.gc = gc
        self.dinuc = dinuc
        self.diaa = diaa
    
        
def parse_fasta(file):
    """
    Parses fasta files to remove any line breaks in the seqsuence and returns as
    a dictionary with the key being the sequence name and the value being the
    seqsuence as a list.
    
    INPUT: fasta file
    OUSTPUT: seqs dictionary in format seqs[name][seqsuence as a list]
    """
    
    with open(file, 'r') as f:
        lines = f.readlines()
        seqs = {}
        raw_seqs = ''
        name = 'UKNOWN SEQUENCE'
        for line in lines:
            if '>' in line:
                if raw_seqs:
                    seqs[name] = list(raw_seqs)
                    raw_seqs = ''
                name = line.strip()
            else:
                raw_seqs += line.strip()

        if raw_seqs: 
            seqs[name] = list(raw_seqs)
    return seqs
    


def gc_calc(seq):
    """
    Calculates the % GC content of a seqsuence from an input sequence and rounds
    to two decimal points. Sequence can be either a list or a string
    
    INPUTS: DNA sequence as list or string
    OUTPUT: GC percentage as a number between 0 to 100%
    """
    if type(seq) == str:
        seq = list(seq)
    length = len(seq)
    gc = 0
    for base in seq:
        if base in ['C', 'c', 'G', 'g']:
            gc +=1
    pgc = (gc / length) * 100
    pgc = round(pgc, 4)
    return pgc


#dinucleotide calculations:
def dinuc_calc(seq):
    """
    Calcuation of all dinucleotide ratios as a percentage to two decimal points
    for the 5'-3' strand of the sequences in the seqs dictinary.
    
    INPUT: seqs[sequence name][sequence as list]
    OUTPUT: dinno[dinucleotide][]
    """

    length = len(seq)
    dinno = {}
    for i in range(0, len(seq)-1):
        di = seq[i]+seq[i+1]
        if not di in dinno:
            dinno[di] = 1
        else:
            dinno[di] +=1
            
    for din in dinno:
        dinno[din] = round((dinno[din]/length) * 100, 4)
    return dinno


def dna_translate(DNA):
    """
    Establishes a DNA translation function for use in future functions.
    
    INPUT: DNA sequence as a string or list
    OUTPUT: aa_seq, holding the amino acid sequence as a string
        
    """
    aa_seq = ''
    codon_table = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'}    

    for i in range(0, len(DNA), 3):
        codon = DNA[i:i+3]
        if codon in codon_table:
            aa = codon_table[codon]
            aa_seq += aa
        else:
            aa_seq += 'X'
    if aa_seq[-1] == 'X':
        aa_seq = aa_seq[:-1]
    return aa_seq
            
            
def reading_frames(seq):
    """
    Determines the three 5'-3' reading frames of a dna sequence. Returns a 
    dictionary with all three reading frames as a string
    
    INPUT: DNA sequence as a list
    OUTPUT: rfs{'rf1': reading_frame_1, 'rf2': reading_frame_2,
                'rf3': reading_frame_3}
    """
    seq = ''.join(seq)
    rfs = {}
    rfs['rf1'] = dna_translate(seq)
    rfs['rf2'] = dna_translate(seq[1:])
    rfs['rf3'] = dna_translate(seq[2:])   
    return rfs


     

# Di-amino acid counts:
def di_amino_calc(rfs):
    """
    Calculates di-amino acid frequencies as a percentage to the 6th decimal 
    place. Input is reading frames 'rfs' dictionary returned from reading_frames
    function. Returns a dictionary of dictionaries containing percentages of
    di-amino frequencies for each reading frame.
    
    INPUT: rfs{'rf1': reading_frame_1, 'rf2': reading_frame_2,
                'rf3': reading_frame_3}
    
    OUTPUT: di_rfs{rf1: {dirf_1}, rf2:{dirf_2}, rf3:{dirf_3}}
            dirf template: dirf{di-amino : frequency of di-amino}
    """
    
    di_rfs = {}
    for frame in rfs:
        rf = list(rfs[frame])
        length = len(rf)
        dirf = {}
        for aa in range(0, len(rf)-1):
            di_aa = rf[aa] + rf[aa+1]
            if not di_aa in dirf:
                dirf[di_aa] = 1
            if di_aa in dirf:
                dirf[di_aa] += 1
        
        for entry in dirf:
            dirf[entry] = round((dirf[entry]/length)*100, 6)
        di_rfs[frame] = dirf           
    di_rfs_sorted = [(frame, di_amino, frequency) for frame in di_rfs for di_amino, frequency in di_rfs[frame].items()]
    di_rfs_sorted = sorted(di_rfs_sorted, key=lambda item: item[1], reverse=True)

    return di_rfs_sorted

def plot_frequencies(entries, file_path, simplify_q):
    """
    Plots the frequencies of the dinucleotide and diamino pairs as well as a 
    comparison graph of the %GC content of all analyzed sequences and saves as
    a bar plot.
    The %GC content data of all species is plotted in a single graph, while the
    dinucleotide and diamino frequencies for each species are plotted separately
    but saved in a single png file.
    
    INPUTS: entries resulting from full_dna_analysis
            file_path: folder containing genomes that are being analyzed
            simplify_q: answer the the question 'Would you like to simplify 
            NCBI genome names for the distance matrix? [y/n]'
            
    OUTPUT: 1 png image containing a bar plot of the %GC of all sequences
                analyzed
            1 png image containing the bar plots of all dinucleotide and diamino
                frequencies of a single sequence
            All images are saved in:
                '<file_path>/distance_parsing/sequence_anaysis' folder
    """
    names = []
    gc = []
    for entry in entries:
        gc.append(entry.gc)
        name = entry.id
        if '>' in name:
            name = name.split(' ')
            name = name[1] + '_'+ name[2]
        names.append(name)
        dinucs = []
        diaa_rf1 = []
        diaa_rf2 = []
        diaa_rf3 = []
        
        dinuc_freq = []
        diaa_freq1 = []
        diaa_freq2 = []
        diaa_freq3 = []
        
        for dinuc, freq in entry.dinuc.items():
            dinucs.append(dinuc)
            dinuc_freq.append(freq)
            
        for rf, diaas, freq in entry.diaa:
            if rf == 'rf1':
                diaa_rf1.append(diaas)
                diaa_freq1.append(freq)
            if rf == 'rf2':
                diaa_rf2.append(diaas)
                diaa_freq2.append(freq)
            if rf == 'rf3':
                diaa_rf3.append(diaas)
                diaa_freq3.append(freq)
                    
        fig, axes = plt.subplots(nrows = 4, ncols = 1, figsize = (100,20))
        axes[0].bar(dinucs, dinuc_freq, color = 'blue')
        axes[0].set_xlabel('DINUCLEOTIDES')
        axes[0].set_ylabel('FREQUENCY')
        axes[0].set_xticks(range(len(dinucs)))
        axes[0].set_xticklabels(dinucs, rotation=90)
        
        axes[1].bar(diaa_rf1, diaa_freq1, color = 'orange')
        axes[1].set_xlabel('DI-AMINO ACID')
        axes[1].set_ylabel('FREQUENCY')
        axes[1].set_xticks(range(len(diaa_freq1)))
        axes[1].set_xticklabels(diaa_rf1, rotation=90)
        
        axes[2].bar(diaa_rf2, diaa_freq2, color = 'red')
        axes[2].set_xlabel('DI-AMINO ACID')
        axes[2].set_ylabel('FREQUENCY')
        axes[2].set_xticks(range(len(diaa_freq2)))
        axes[2].set_xticklabels(diaa_rf2, rotation=90)
        
        axes[3].bar(diaa_rf2, diaa_freq2, color = 'yellow') 
        axes[3].set_xlabel('DI-AMINO ACID')
        axes[3].set_ylabel('FREQUENCY')
        axes[3].set_xticks(range(len(diaa_freq3)))
        axes[3].set_xticklabels(diaa_rf2, rotation=90)
        
        plt.tight_layout()
        out_file = file_path +  name + '_dinuc-diaa_charts.png'
        plt.savefig(out_file)
    
    plt.figure(figsize = (10,20))
    plt.bar(names, gc, color = 'green')
    plt.xticks(rotation = 90)
    plt.xlabel('SEQUENCE')
    plt.ylabel('%GC CONTENT')
    plt.tight_layout()
    out_file = file_path + '%GC_AllSeqs_chart.png'
    plt.savefig(out_file)
    
def full_dna_analysis(file, save_file_q):
    """
    Provides a full sequence analysis of an input DNA fasta file. Outputs a txt file with 
    the %GC content, di-nucleotide frequency and di-amino frequency of all 5'-3' 
    reading frames for all sequnces within the fasta. Will also return the resulting 
    text will and print on screen.
    
    INPUT: Fasta file with DNA sequence(s)
    OUTPUT: write analysis_input_file.txt containing %GC, di-nucleotide 
                and di-amino freqs
            returns text from resulting analysis file
            prints text from resulting analysis file
    """
    file = file.replace('\\', '/')
    file_path = file.split('/')

    out_name = 'SeqAnalysis_' + file_path[-1]

    file_path = file_path[0:-1]
    file_path = '/'.join(file_path)  + '/belvu_distance_parsing/'

    out_file = file_path + out_name
    seqs = parse_fasta(file)
    
    entries = []
    
    with open(out_file, 'w') as f:
        for entry in seqs:
            seq = seqs[entry]
            class_id = entry
            name = entry + '\n'

            pgc = gc_calc(seq)
            class_gc = pgc
            gc = '\n%GC Content: ' + str(pgc) + '%\n'
            
            dinu_stxt = '\nDi-Nucleotide Frequencies:'
            dn = dinuc_calc(seq)
            class_dinuc = dn
            dinuc = ''
            for din, freq in dn.items():
                dinuc += '\n\t' + str(din) + ': ' + str(freq) + '%'
            
            diaa_stxt = '\n\n\nDi-Amino Acids Frequencies: \n\nReading Frame #1\t\tReading Frame #2\t\tReading Frame #3 \n'
            
            di_aa = ''
            di_aa1 = []
            di_aa2 = []
            di_aa3 = []
            rfs = reading_frames(seq)
            di_rfs = di_amino_calc(rfs)
            class_diaa = di_rfs
            for frame, diaa, freq in di_rfs:  # Unpack the tuple
                if frame == 'rf1':
                    di_aa1.append(diaa + ': ' + str(freq) + '%\t\t\t')
                if frame == 'rf2':
                    di_aa2.append(diaa + ': ' + str(freq) + '%\t\t\t')
                if frame == 'rf3':
                    di_aa3.append(diaa + ': ' + str(freq) + '%\n')
            
            maxs = max(len(di_aa1), len(di_aa2), len(di_aa3))
            
            for i in range(0, maxs):
                di_aa += di_aa1[i] if i < len(di_aa1) else '\t\t\t\t'
                di_aa += di_aa2[i] if i < len(di_aa2) else '\t\t\t\t'
                di_aa += di_aa3[i] if i < len(di_aa3) else '\n'
                       
            final = name + gc + dinu_stxt + dinuc + diaa_stxt + di_aa + '\n\n\n\n'
            if save_file_q== 'y':
                f.write(final)
            print(final)
            entries.append(AnalyzedSeq(class_id, class_gc, class_dinuc, class_diaa))
    return entries
    

def gc_distances(entries):
    """
    Creates a dictionary with the %GC distances between each pair of sequences 
    parsed through the full_dna_analysis function of exercise 1. Requires all
    of the sequences to be of 'AnalyzedSeq' object class. Returns a dictionary
    containing the sequence id as key and distance as value.
    
    %GC distance = sqrt( (%GC_seq1 - %GC_seq2)^2 )
    
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
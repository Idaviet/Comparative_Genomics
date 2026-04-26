"""
Imports the scripts from exercises 1, 2 and 3 to create a single argparse function to perform
the sequence analysis, distance scoring/matrices and ORF prediction from a single python command.
Some functions from the previous files are re_written take in additiona inputs to  allow for a 
more seemless use of the terminal command line 

Script contains an argparse function to call belvu_dist_parsing function from
a terminal.

INPUT FOLLOWING COMMAND ON TERMINAL AND FOLLOW ON SCREEN INSTRUCTIONS
    >python FinalProject_group2.py <folder containing fasta files to be analyzed>

FUNCTIONS:
    
    the_whole_enchilada: uses all the functions from the previous folders to generate the files
                         returned from the full_dna_parsing, belvu_dist_parsing, and 
                         full_orf_parsing functions, re-written here so that the questions from 
                         the respective modules can be all answered at the start of command,
                         thus allowing for all processing to occur without interuption between
                         the modules
    
"""

import argparse
import shutil
import os
import sequence_analysis_2 as ex1
import distance_matrix as ex2
import orf_prediction as ex3

# class AnalyzedSeq:
#     def __init__(self, id, gc, dinuc, diaa):
#         self.id = id
#         self.gc = gc
#         self.dinuc = dinuc
#         self.diaa = diaa

       
def combine_selected_fasta(file_path, files):
    """
    Rewriting of combine_selected_fasta function from ex1 to take an extra 'files' variable
    
    Combines select fasta files in a folder into a single file named 
    'all_selected_sequences.fasta' located in the submitted folder. Once function
    is called with target folder, user will be asked to submit the files to be
    processed, separated by a space. 
    Eg: <sequence1.fasta sequence2.fasta... sequenceN.fasta>
    
    INPUT: filepath of folder containing files to be analyzed
           files string containing all files to be combined separated by a space
           Eg: 'sequence1.fasta sequence2.fasta... sequenceN.fasta'
           
    OUTPUT: all_sequences.fasta file containg all sequences in all fasta files
    
    """
    file_path = file_path.replace('\\', '/')
    
    selected_genomes = file_path + '/all_selected_sequences.fasta'
    selected_sequences = files
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
    
    
    
def full_dna_analysis(file, q1):
    """
    Rewriting of full_dna_analysis function from ex1 to take q1 answer variable

    Provides a full sequence analysis of an input DNA fasta file. Outputs a txt file with 
    the %GC content, di-nucleotide frequency and di-amino frequency of all 5'-3' 
    reading frames for all sequnces within the fasta. Will also return the resulting 
    text will and print on screen.
    
    INPUT: Fasta file with DNA sequence(s)
           q1: answer to fda_q in the_whole_enchilada function
           
    OUTPUT: write analysis_input_file.txt containing %GC, di-nucleotide 
                and di-amino freqs
            returns text from resulting analysis file
            prints text from resulting analysis file
    """
    file = file.replace('\\', '/')
    file_path = file.split('/')

    out_name = 'SeqAnalysis_' + file_path[-1]

    file_path = file_path[0:-1]
    file_path = '/'.join(file_path)  + '/distance_parsing/sequence_analysis/'

    out_file = file_path + out_name
    seqs = ex1.parse_fasta(file)
    
    entries = []
    
    for entry in seqs:
        seq = seqs[entry]
        class_id = entry
        name = entry + '\n'

        pgc = ex1.gc_calc(seq)
        class_gc = pgc
        gc = '\n%GC Content: ' + str(pgc) + '%\n'
        
        dinu_stxt = '\nDi-Nucleotide Frequencies:'
        dn = ex1.dinuc_calc(seq)
        class_dinuc = dn
        dinuc = ''
        for din, freq in dn.items():
            dinuc += '\n\t' + str(din) + ': ' + str(freq) + '%'
        
        diaa_stxt = '\n\n\nDi-Amino Acids Frequencies: \n\nReading Frame #1\t\tReading Frame #2\t\tReading Frame #3 \n'
        
        di_aa = ''
        di_aa1 = []
        di_aa2 = []
        di_aa3 = []
        rfs = ex1.reading_frames(seq)
        di_rfs = ex1.di_amino_calc(rfs)
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
        entries.append(ex1.AnalyzedSeq(class_id, class_gc, class_dinuc, class_diaa))
        print(class_id, ' HAS BEEN CHARACTERIZED\n')        
        if q1 == 'y':
            with open(out_file, 'a') as f:
                ex1.plot_frequencies(entries, file_path, q1)
                final = name + gc + dinu_stxt + dinuc + diaa_stxt + di_aa + '\n\n\n\n'
                f.write(final)
    print('ALL SEQUENCES HAVE BEEN CHARACTERIZED\n\n\n')
    return entries


                
def belvu_dist_parsing(file, folder, fda_q, q1, q2, q3, q4):
    """
    Rewriting of full_dna_analysis function from ex1 to take fda_q1 and belvu_q1-24 answer 
    variables

    Full parsing of all sequences held within input folder. Results in the a
    full sequence analysis file and the bar plot graphs for all possible 
    distance scores determined from the previous functions. Includes an option 
    to extract a belvu matrix for one type of score, or all.
    
    INPUT: folder containing sequences wished to be analyzed
           answers for fda_q and belvu_q1-4 from the argparse command
           
    OUTPUT: bar plots containing scores for all 5 possible score types
            Sequence analysis txt file containing %GC, dinucleotide and di-amino
            frequecies for all sequences analyzed
            OPTIONAL: belvu distance matrix for selected or all types.
    """
    print('\n\nPARSING DISTANCE SCORES...\n')
    full_analysis = q1
    type = []
    if full_analysis == 'y':
        inputs = q2 
        if inputs == 'gc':
            type.append('%GC')
        elif inputs == 'dn':
            type.append('dinuc')
        elif inputs == 'da':
            rf = q3
            diaa = 'diaa_rf' + rf
            type.append(diaa)
        elif inputs == 'all':
            type = ['%GC','dinuc', 'diaa_rf1', 'diaa_rf2', 'diaa_rf3']
    
    simplify = q4
    
    out_folder = folder + '/distance_parsing'
    seqan_folder = out_folder + '/sequence_analysis'
    matrix_folder = out_folder + '/belvu_matrices'
    os.mkdir(out_folder)
    if fda_q == 'y':
        os.mkdir(seqan_folder)
    
    
     
   
    entries = full_dna_analysis(file, fda_q)
    
    if simplify == 'y':
        entries = ex1.ncbi_name_simplifier(entries)
    
    all_dist = ex2.plot_distances(entries, out_folder)
    
    if full_analysis == 'hold':
        q = input('\n\nWould you like to extract a belvu distance matrix?\n\t[y/n]\n')
        if q == 'y':
            file_q = input('Parse all files in folder?\n\t[y/n]\n')
            
            if file_q == 'y':
                ex1.combine_all_fasta(folder)
                ds_fasta = folder + '/all_sequences.fasta'
            
            if file_q == 'n':
                files= input('Please input file names for DISTANCE SCORING (with extensions) separated by a space OR input [all] to parse all files:\n\n')                  
                if files == 'all':
                    combine_selected_fasta(folder, files )
                    fasta = folder + '/all_sequences.fasta'
                    
                if files != 'all':
                    combine_selected_fasta(folder, files)
                    fasta = folder + '/combined_dist_sequences.fasta'
                    ds_fasta = folder + '/all_selected_sequences.fasta'
                    os.rename(ds_fasta, fasta)
                        
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
            
            entries = full_dna_analysis(fasta, 'n')
            if simplify == 'y':
                entries = ex1.ncbi_name_simplifier(entries)
            
            os.mkdir(matrix_folder)    
            all_dist = ex2.plot_distances(entries, out_folder)
   
    if full_analysis == 'y':
        os.mkdir(matrix_folder)
    if full_analysis == 'y' or 'hold':
        
        for i in type:
            output = matrix_folder + '/' + i + '_belvu_distance_matrix.fasta'
            dist_list = ex2.extract_dist_list(all_dist, i)
            ex2.belvu_matrix(dist_list, output)
            
            
            
            
def full_orf_parsing(folder, file, minimum_orf_length):
    """
    Rewriting of combine_selected_fasta function from ex1 to take an extra 'file' variable
    
    Full parsing of genomes held within a folder. Can choose to parse all files
    contained in said folder, or to choose a custom selection of files. Results
    in a new folder name 'Sequence_ORFs' containing the ORF fasta files of all
    selected files. Will also output a summary file containing the number of 
    ORFs predicted for each sequence and a compiled fasta file containing all
    the selected sequences.
    
    INPUT: folder containing fasta files to be analyzed
           fasta file containing the sequence to be analyzed
           
    OUTPUT: folder with ORF predictions for all sequences selected in separate
            fasta files
            Summary file containing number of predicted ORFs for each sequence
            Combined fasta file containing all the selected fasta files
    """
    print('\n\nPREDICTING ORFs...')
    outfolder = folder + '/Sequence_ORFs'
    os.mkdir(outfolder)
   
    all_orfs = ex3.extract_orfs(file, minimum_orf_length)
    
    summary = ex3.ORF_summary(all_orfs)
    ex3.write_orf_fasta(all_orfs, outfolder)
    sum_file = outfolder + '/ORF_summary.txt'
    
    with open(sum_file, 'w') as f:
        for species in summary:
            forward = summary[species][0]
            reverse = summary[species][1]
            total = forward + reverse
            entry = species + ':\n\t' + str(forward) + '\t\tORFs in forward strand\n\t' + str(reverse) +'\t\tORFs in reverse strand\n\t' +str(total) +'\tORFs total\n\n' 
            f.write(entry)
            print(entry)
    print('ORFS FOR ALL SELECTED FILES HAVE BEEN COMPILED\n\n\n')


def the_whole_enchilada(folder):
    """
    Function combining all functions from previous modules and the rewritten functions above
    to be used for full sequence analysis of fasta file in the given folder
    
    INPUT: folder containing fasta files
    
    OUTPUTS: folder(s) containing results of belvu distance scoring and/or orf prediction
    """
    folder = folder.replace('\\', '/')   
    belvu_q1, belvu_q2, belvu_q3, belvu_q4 = 'none'
    
    full_q = input('\nSelect Analysis to be conducted:\n\tDistance Scoring/Matrix [ds]\n\tOpen Reading Frame Prediction [orf]\n\tAll [all]\n\n')

    existing_files = [os.path.join(folder, f) for f in os.listdir(folder)]
    results_folder = os.path.join(folder, 'CompGen_results')
    
    if full_q == 'orf' or full_q == 'all':
        minimum_orf_length = input('\nPlease input minimum ORF length for ORF parsing:\n(NOTE: PRESS ENTER TO USE DEFAULT 300 NUCLEOTIDE ORF LENGTH)\n\n')
        if minimum_orf_length == '':
            mol = 300
        else:
            mol = int(minimum_orf_length)
    
    q0 = input('\nParse all files in folder?:\n\t[y/n]\n\n')
    if q0== 'n':
        if full_q == 'ds' or full_q == 'all':
            q1 = input('\nPlease input file names for DISTANCE SCORING (with extensions) separated by a space OR input [all] to parse all files:\n\n')
        if full_q == 'orf' or full_q == 'all':
            q2 = input('\nPlease input file names for ORF PREDICTION (with extensions) separated by a space OR input [all] to parse all files:\n\n')
        
    if full_q == 'ds' or full_q == 'all':
        belvu_q4 = input('\nWould you like to simplify NCBI genome names for the distance matrix? \n\t[y/n]\n')

        fda_q = input('\nSave full analysis files? (contains txt files + bar plots of %GC content, dinucleotide & diamino acid freqs):\n\t[y/n]\n\n')
        
        belvu_q1 = input('\nParse belvu distance matrix, or hold until seq analysis is completed? \n\t[y/n/hold]\n')
        if belvu_q1 == 'y':
            belvu_q2 = input('\nPlease enter which distance scores to use:\n\t%GC content [gc]\n\tDinulceotide Frequency [dn]\n\tDiamino Frequency [da]\n\tAll matrices [all]\n')  
            if belvu_q2 == 'da':
                belvu_q3 = input('\nWhich reading frame? [1/2/3] \n')
            
    if q0 == 'y':
        ex1.combine_all_fasta(folder)
        ds_fasta = folder + '/all_sequences.fasta'
        orf_fasta = ds_fasta
    
    if full_q == 'ds' or full_q == 'all':
        if q0 == 'n' :
            if q1 == 'all':
                ex1.combine_all_fasta(folder)
                ds_fasta = folder + '/all_sequences.fasta'
            if q1 != 'all':
                combine_selected_fasta(folder, q1)
                ds_rename = folder + '/all_dist_sequences.fasta'
                ds_fasta = folder + '/all_selected_sequences.fasta'
                os.rename(ds_fasta, ds_rename)
                ds_fasta = ds_rename
                
    if full_q == 'orf' or full_q == 'all': 
        if q0 == 'n':
            if q2 == 'all':
                ex1.combine_all_fasta(folder)
                orf_fasta = folder + '/all_sequences.fasta'
            if q2 != 'all':
                combine_selected_fasta(folder, q2)
                orf_fasta = folder + '/all_selected_sequences.fasta'
                orf_rename = folder + '/all_orf_sequences.fasta'
                orf_fasta = folder + '/all_selected_sequences.fasta'
                os.rename(orf_fasta, orf_rename)
                orf_fasta = orf_rename
        
    
    if full_q == 'orf' or full_q == 'all':
        full_orf_parsing(folder, orf_fasta, mol)
        
    if full_q == 'ds' or full_q == 'all':
        belvu_dist_parsing(ds_fasta, folder, fda_q, belvu_q1, belvu_q2, belvu_q3, belvu_q4)
        

    for f in os.listdir(folder):
        file = os.path.join(folder, f)
        if file not in existing_files:
            os.makedirs(results_folder, exist_ok = True)
            shutil.move(file, results_folder)

    print('\n\n\n------------RUN COMPLETE-----------')
    

if __name__ == "__main__":
    """
    Argument parser function to generate Distance matrices and potential ORFs for all sequences 
    within a folder from computer terminal.
    
    INPUT FORMAT:
    >python FinalProject_group2.py <input folder containing fasta files>
    """
    parser = argparse.ArgumentParser(description="Generates Distance Scores/matrices and potential ORFs in a\
                                     fasta file format from selected fasta files within a given folder", \
                            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('folder', type=str, 
                        help='Please input folder containing fasta files to be analyzed')
                    
    args = parser.parse_args()
    the_whole_enchilada(args.folder) 
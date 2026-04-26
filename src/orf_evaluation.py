import sys
file_location = r"C:\Users\isaac\Desktop\SU Semester 2023\Comparative Genomics\Final Project\FinalProject_group2"
sys.path.append(file_location)
import orf_prediction as ex3
import sequence_analysis_2 as ex1
from Bio.Seq import Seq
from Bio import pairwise2
import argparse

#orf_file = r'C:\Users\isaac\Desktop\SU Semester 2023\Comparative Genomics\Final Project\FinalProject_group2\genomes\Leuconostoc.fasta'
#rp_file = r'C:\Users\isaac\Desktop\SU Semester 2023\Comparative Genomics\Final Project\FinalProject_group2\13_ref_prot'

def parse_reference(rp_file, orf_file):
    print('\n\n')
    ref_prot = ex1.parse_fasta(rp_file)
    orf_gen = ex1.parse_fasta(orf_file)
    for_rf = []
    rev_rf = []
    rp_orfs = {}

    for entry, seq in orf_gen.items():
        sequence = Seq(''.join(seq))
        rev_comp = sequence.reverse_complement()
        for i in range(3):
            rf = str(sequence[i:].translate())
            rcrf = str(rev_comp[i:].translate())
            for_rf.append(rf)
            rev_rf.append(rcrf)



    for orf, seq in ref_prot.items():
        sequence = ''.join(seq)
        if '_rev' in orf:
            for rf in rev_rf:
                orf_name = orf+ 'ORF_'
                start = rf.find(sequence) * 3
                end = start + len(sequence) * 3 - 1
                orf_name += str(start) + '-' + str(end)
                if start != -3:
                    rp_orfs[orf_name] = sequence
        else: 
            for rf in for_rf:
                orf_name = orf+ 'ORF_'
                start = rf.find(sequence) * 3
                end = start + len(sequence) * 3 - 1
                orf_name += str(start) + '-' + str(end)
                if start != -3:
                    rp_orfs[orf_name] = sequence



    print('\n\n# of Reference ORFs: ',len(rp_orfs))
    return rp_orfs



def orf_performance_evaluation(rp_file, orf_file, minimum_orf_length, true_positive_threshold):
   
    mol = int(minimum_orf_length)
    tpt = int(true_positive_threshold)
    
    predicted_orfs = ex3.extract_orfs(orf_file, mol)
    ref_orfs = parse_reference(rp_file, orf_file)
    p_orfs = {}
    
    for entry, orfs in predicted_orfs.items():
        for orf, sequence in orfs.items():
            dna_seq = Seq(sequence)
            protein_seq = str(dna_seq.translate())
            p_orfs[orf]=protein_seq
        
       
    print('\n\nCALCULATING TRUE POSITIVES, PRECISION AND RECALL VALUES...\n\n')
    pcids=[]
    tracker  = []
    print('# of ORFS predicted: ', len(p_orfs), '\n\n\n')

    i = 0
    range1 = 10
    range2 = 11
    for p_orf, seq1 in p_orfs.items():
        i +=1
        pc_done = ((i)/len(p_orfs)) *100
        if range1< pc_done < range2:
            print(range1, '% Complete')
            range1 += 10
            range2 += 10
            

        p_pos = p_orf.split('ORF_')
        if len(p_pos) != 2:
            continue 
        p_pos = p_pos[1]
        p_pos = p_pos.split('-')
        if len(p_pos) != 2:
            continue 
        p_start = p_pos[0]
        p_end = p_pos[1]
        if not p_start.isdigit() or not p_end.isdigit():
            continue
        p_start = int(p_start)
        p_end = int(p_end)


        for rp_orf, seq2 in ref_orfs.items():
            if rp_orf not in tracker:                        

                rp_pos = rp_orf.split('ORF_')
                if len(rp_pos) != 2:
                    continue
                rp_pos = rp_pos[1]
                rp_pos = rp_pos.split('-')
                if len(rp_pos) != 2:
                    continue
                rp_start = rp_pos[0]
                rp_end = rp_pos[1]
                
                if not rp_start.isdigit() or not rp_end.isdigit():
                    continue
                
                rp_start = int(rp_start)
                rp_end = int(rp_end)
                
                if p_start in range(rp_start, rp_end) or p_end in range(rp_start, rp_end) or rp_end in range(p_start, p_end) or rp_start in range(p_start, p_end):

                    alignments = pairwise2.align.globalxx(seq1, seq2)
                    
                    best_alignment = alignments[0]
                    aligned_seq1 = best_alignment[0]
                    aligned_seq2 = best_alignment[1]
                    
                    alignment_length = len(aligned_seq1)
                    identity_count = sum(aa1 == aa2 for aa1, aa2 in zip(aligned_seq1, aligned_seq2))
                    percent_identity = (identity_count / alignment_length) * 100
                    
                    if percent_identity >= tpt:
                        tracker.append(rp_orf)
                        pcids.append(percent_identity)
                        break
    
    average = round(sum(pcids)/len(pcids), 2)
    
    tps = len(pcids)
    
    fns = len(ref_orfs) - tps
    
    fps = len(p_orfs) - tps
    
    precision = (tps / (tps + fns)) *100
    
    recall = (tps / (tps + fns + fps)) *100
    
    print('\n\n# OF PREDICTED ORFS: ', len(p_orfs), '\n# OF ORFS IN REFERENCES: ', len(ref_orfs))
    print('\n\n# OF TRUE POSITIVES: ',tps)
    print('# OF FALSE POSITIVES: ',fps)
    print('# OF FALSE NEGATIVES: ',fns)
    print('\nAVERAGE %IDENTITY OF TRUE POSITIVES: ',average, '%')
    
    print('\n\nPRECISION: ',round(precision, 2), '%')
    print('RECALL: ',round(recall, 2), '%')


if __name__ == "__main__":
    """
    Argument parser function to generate analysis file from computer terminal.
    
    INPUT FORMAT:
    >python ex1.py <input_fasta_file>
    """
    parser = argparse.ArgumentParser(description="Generate DNA analysis file for \
                            fasta input file. Includes %GC, Di-nucleotide \
                                and Di-amino acid frequencies for all 5'-3' frames", \
                            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('Reference_Proteome', type=str, 
                        help='Please input Reference Proteome used to determine positives')  
    parser.add_argument('Genome', type=str, 
                        help='Please input fasta file containing the genome of the reference proteome')  
    
    parser.add_argument('mol', type=str, 
                        help='Please input minimum orf length in base pairs')  
    
    parser.add_argument('tpt', type=str, 
                        help='Please input % identity threshold for an alignment to be considered a true positive')  
    
    args = parser.parse_args()
    orf_performance_evaluation(args.Reference_Proteome, args.Genome, args.mol, args.tpt)
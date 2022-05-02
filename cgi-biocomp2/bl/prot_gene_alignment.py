# -*- coding: utf-8 -*-
#!/usr/bin/python3
"""
Created on Fri Apr 29 14:52:10 2022
This is the script to produce a protein sequence: gene sequence alignment
@author: angli
"""
def prot_gene_alignment(protein_seq, gene_sequence):
    """
    This function takes in a protein sequence and genomic sequence and aligns them.
    If there is mismatch in length, the alignment does not occur.

    Parameters
    ----------
    protein_seq : str
        DESCRIPTION.
    gene_sequence : str
        DESCRIPTION.

    Returns
    -------
    List of tuples of format (amino acid, codon)
    """
    alignment_list =[]
    if len(protein_seq) != (len(gene_sequence)/3):
        alignment_list.append("ERROR: The sequences cannot be aligned")
    else:
        for i in range(0,len(protein_seq)):
            alignment_list.append((protein_seq[i], gene_sequence[i*3:(i*3)+3].upper()))
            
    master_codon_tuples = [
        ("ATT","I"), ("ATC", "I"), ("ATA", "I"),
        ("CTT", "L"), ("CTC", "L"), ("CTA", "L"), ("CTG", "L"), ("TTA", "L"), ("TTG", "L"),
        ("GTT", "V"), ("GTC", "V"), ("GTA", "V"), ("GTG", "V"),
        ("TTT", "F"), ("TTC", "F"),
        ("ATG", "M"),
        ("TGT", "C"), ("TGC", "C"),
        ("GCT", "A"), ("GCC", "A"), ("GCA", "A"), ("GCG", "A"),
        ("GGT", "G"), ("GGC", "G"), ("GGA", "G"), ("GGG", "G"),
        ("CCT", "P"), ("CCC", "P"), ("CCA", "P"), ("CCG", "P"),
        ("ACT", "T"), ("ACC", "T"), ("ACA", "T"), ("ACG", "T"),
        ("TCT", "S"), ("TCC", "S"), ("TCA", "S"), ("TCG", "S"), ("AGT", "S"), ("AGC", "S"),
        ("TAT", "Y"), ("TAC", "Y"),
        ("TGG", "W"),
        ("CAA", "Q"), ("CAG", "Q"),
        ("AAT", "N"), ("AAC", "N"),
        ("CAT", "H"), ("CAC", "H"),
        ("GAA", "E"), ("GAG", "E"),
        ("GAT", "D"), ("GAC", "D"),
        ("AAA", "K"), ("AAG", "K"),
        ("CGT", "R"), ("CGC", "R"), ("CGA", "R"), ("CGG", "R"), ("AGA", "R"), ("AGG", "R"),
        ("TAA", "Stop"), ("TAG", "Stop"), ("TGA", "Stop")
        ]
    codon_error_check_counter = 0
    for (aa, codon) in alignment_list:
        if (codon, aa) in master_codon_tuples:
            codon_error_check_counter +=1
    if codon_error_check_counter == len(protein_seq):
        return(alignment_list)
    else:
        return("Error checking failed, alignment not sucessful")
    
#-----------------------
#Test code
#ps = "LLDGKEIKRLNVQWLRAHLGIVSQEPILFDCSIAENIAYGDNSRVVSQEEIVRAAKEANIHAFIESLPN"
#gs  = "ctgcttgatggcaaagaaataaagcgactgaatgttcagtggctccgagcacacctgggcatcgtgtcccaggagcccatcctgtttgactgcagcattgctgagaacattgcctatggagacaacagccgggtggtgtcacaggaagagatcgtgagggcagcaaaggaggccaacatacatgccttcatcgagtcactgcctaat"
#test = prot_gene_alignment(ps, gs)

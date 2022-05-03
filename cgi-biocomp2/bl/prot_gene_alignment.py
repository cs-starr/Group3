#!/usr/bin/python3.7.9
# -*- coding: utf-8 -*-

"""
Created on Sat Apr 30 12:37:54 2022
This function
@author: angli
"""

def prot_gene_alignment(protein_seq, gene_sequence):
    """
    This function takes in a protein sequence and genomic sequence and aligns them.
    If there is mismatch, the function will attempt to match up to any codon discrepancy.

    Parameters
    ----------
    protein_seq : str
        The protein sequence to be aligned
    gene_sequence : str
        DThe DNA sequence to be aligned

    Returns
    -------
    List of tuples of format (amino acid, codon)
    """
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
    alignment_list =[]
    if len(protein_seq) != (len(gene_sequence)/3):
        for i in range(0,len(protein_seq)):
            if (gene_sequence[i*3:(i*3)+3].upper(),protein_seq[i]) in master_codon_tuples:
                alignment_list.append((protein_seq[i], gene_sequence[i*3:(i*3)+3].upper()))
            else:
                break
        #The for loop allows for potential counting of which position the break has occured.
    else:
         for i in range(0,len(protein_seq)):
             if (gene_sequence[i*3:(i*3)+3].upper(),protein_seq[i]) in master_codon_tuples:
                 alignment_list.append((protein_seq[i], gene_sequence[i*3:(i*3)+3].upper()))
    
    return(alignment_list)

"""
#Test code
ps = "LLDGKEIKRLNVQWLRAHLGIVSQEPILFDCSIAENIAYGDNSRVVSQEEIVRAAKEANIHAFIESLPN"
gs  = "ctgcttgatggcaaagaaataaagcgactgaatgttcagtggctccgagcacacctgggcatcgtgtcccaggagcccatcctgtttgactgcagcattgctgagaacattgcctatggagacaacagccgggtggtgtcacaggaagagatcgtgagggcagcaaaggaggccaacatacatgccttcatcgagtcactgcctaat"
test = prot_gene_alignment(ps, gs)
print(test)
print(len(ps))
print(len(gs)/3)
"""

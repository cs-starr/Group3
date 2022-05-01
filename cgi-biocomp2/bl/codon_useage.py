# -*- coding: utf-8 -*-
#!/usr/bin/python3
"""
Created on Fri Apr 29 14:52:10 2022
This is the script for calculation of codon useage within a sequence
@author: angli
"""
def codon_useage(coding_dna_sequence):
    from collections import Counter
    
    raw_seq = coding_dna_sequence.upper()
    
    len_rs = len(raw_seq)
    exp_iter = int(len_rs/3)
    current_iter = 0
    seq_slice = slice(0,3)
    
    codon_list = []
    
    while current_iter != exp_iter:
        codon = raw_seq[seq_slice]
        codon_list.append(codon)
        if codon == "TAA" or codon == "TAG" or codon == "TGA":
            break
        raw_seq = raw_seq[3:]
        current_iter += 1
    
    seq_total_codons = len(codon_list)
    
    count_of = Counter(codon_list)
        
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
    
    codon_use_data={}
    for (m_codon, AA) in master_codon_tuples:
        codon_use_data.update({m_codon:(count_of[m_codon], AA, (count_of[m_codon]/seq_total_codons)*100)})
    
    return(codon_use_data)

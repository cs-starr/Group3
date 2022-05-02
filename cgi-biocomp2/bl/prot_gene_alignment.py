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
    return(alignment_list)
    
#-----------------------
#Test code
#ps = "LLDGKEIKRLNVQWLRAHLGIVSQEPILFDCSIAENIAYGDNSRVVSQEEIVRAAKEANIHAFIESLPN"
#gs  = "ctgcttgatggcaaagaaataaagcgactgaatgttcagtggctccgagcacacctgggcatcgtgtcccaggagcccatcctgtttgactgcagcattgctgagaacattgcctatggagacaacagccgggtggtgtcacaggaagagatcgtgagggcagcaaaggaggccaacatacatgccttcatcgagtcactgcctaat"
#test = prot_gene_alignment(ps, gs)

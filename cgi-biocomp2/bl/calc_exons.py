# -*- coding: utf-8 -*-
#!/usr/bin/python3

"""
Created on Fri Apr 29 14:50:05 2022

This is the script for calculation of the coding DNA sequence given input information on 
the overall DNA sequence, exon locations, forward or complimentary strand location of the gene and 
whether the uploaded sequence is offset.

@author: angli
"""

def calc_exons(dna_seq, exon_locations, comp_strand, codon_start):
    """
    This function calculates the exons when given information of:
    dna_seq - str - the overall genomic sequence
    exon_location - list of tuples of format ("string","string")
    comp_string - boolean - 1 indicates that the coding sequence is present on the complementary strand
    codon_start - int - the codon start offset, as indicated by the genbank file
    --------
    Returns:
    String containing the coding DNA sequence
    """
    
    temp_exon_store = []
    for (start,stop) in exon_locations:
        if start[0] == "<" or start[0] == ">":
            start = int(start[1:])
        else:
            start = int(start)
        if stop[0] == "<" or stop[0] == ">":
            stop = int(stop[1:])
        else:
            stop = int(stop)
        temp_exon_store.append(dna_seq[(start-1):(stop)])
    working_exon_string = "".join(temp_exon_store)
    
    #Addressing whether the sequence begins after the given co-ordinates. 
    #(Codon start is which nucleotide it start from with 1 being the co-ordinate given within te genbank file)
    trim = (codon_start - 1) 
    final_exon_str = ""
    
    if comp_strand == 1:
        rev_ex_str = working_exon_string[::-1]
        rev_list = []
        final_exon_str = ""
        for i in rev_ex_str:
            if i == "c":
                s = "g"
            if i == "g":
                s = "c"
            if i == "a":
                s = "t"
            if i == "t":
                s = "a"
            rev_list.append(s)
        final_exon_str = "".join(rev_list)
    else:
            final_exon_str = working_exon_string
            
    final_exon_str = final_exon_str[trim:]
    
    return(final_exon_str)

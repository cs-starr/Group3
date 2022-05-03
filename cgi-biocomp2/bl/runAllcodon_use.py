#!/usr/bin/python3.7.9
# -*- coding: utf-8 -*-

"""
Created on Fri Apr 29 14:51:00 2022
This function calculates all coding regions within the database

@author: angli

"""

def runAllcodon_use():
    """
    This function process all entries from the database, calculating the codon useage for each gene.
    This data is collated and written to a text file using a for loop, presenting information in format (and order):
        codon, amino acid code, raw counts and percentage use
    " " is given as a separator between data items, "\n" separates each codon entry.
    
    Returns
    -------
    Writes data into a text file, stored in the same directory
    Returns a confirmatory message :"Calculation complete"
    """
    import bl_api
    import db_api
    import calc_exons as ce

    #Set up codon counter
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
    master_codon_counter=[]
    for (codon_code, aa_code) in master_codon_tuples:
        master_codon_counter.append((codon_code, aa_code, 0))
    
    #Get all entries from DB later
    entry_list = bl_api.getgetAllEntries()
    
    #iterating on accession number
    #db_API output: return[('accession_number', 'gene_id', 'protein name','chromosomal location')]
    master_coding_string_list = []
    
    for acc_number in entry_list[0]:
        temp_data = db_api.sequence(acc_number)
        dna_seq = temp_data[w] #expecting string
        exon_locations = temp_data[x] #expecting list of tuples
        comp_strand = temp_data[y] #expecting comp strand location
        codon_start = temp_data[z] #expecting codon start loc
        exon_str = ce.calc_exons(dna_seq, exon_locations, comp_strand, codon_start)
        master_coding_string_list.append(exon_str)
    
    #break down each coding_string into individual codons
    for coding_string in master_coding_string_list :
        raw_seq = coding_string().upper()
        len_rs = len(raw_seq)
        exp_iter = len_rs//3
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
        
        k = Counter(codon_list)
        for (codon_code, aa_code, count) in master_codon_counter:
            update_count = k[codon_code]
            count = count + update_count
    
    total_count = 0
    for (codon_info) in master_codon_counter:
        total_count = total_count + codon_info[2] 
    
    with open("total_codon_use.txt", 'w') as tcu:
        for (codon_info) in master_codon_counter:
            tcu.write("\n" + str(codon_info[0])+" ")
            tcu.write(str(codon_info[1] + " "))
            tcu.write((str(codon_info[2]) + " "))
            tcu.write(str((codon_info[2]/total_count)*100) + " ")

    return("Calculation complete")

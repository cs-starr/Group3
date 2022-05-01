# -*- coding: utf-8 -*-
#!/usr/bin/python3
"""
Created on Fri Apr 29 14:49:28 2022
This is the script to run the calculation for all codon useage statistics within the database
@author: angli
"""
def runAllcodon_use():
    """
    This function is a one time run event.  To calculate all coding regions within the database
    and to calculate

    Returns
    -------
    Writes data into a text file. (name TBD)
    (Need to include - date/time of calculation)
    """
    import getgetAllEntries
    import bl_api
    import db_API
    import calc_exons
    
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
    master_exon_list = []
    
    for acc_number in entry_list[0]:
        temp_data = db_API.sequence(acc_number)
        dna_seq = temp_data["a"] #expecting string
        exon_locations = temp_data["x"] #expecting list of tuples
        comp_strand = temp_data["y"] #expecting comp strand location
        codon_start = temp_data["z"] #expecting codon start loc
        exon_str = calc_exons(dna_seq, exon_locations, comp_strand, codon_start
        master_exon_list.append(exon_str)
    
    #break down each exon into individual codons
    for exon in master_exon_list :
        raw_seq = exon.upper()
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
        
        k = Counter(codon_list)
        for (codon_code, aa_code, count) in master_codon_counter:
            update_count = k[codon_code]
            count = count + update_count
    
    total_count = 0
    for (codon_info) in master_codon_counter:
        total_count = total_count + codon_info[2] 
    
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

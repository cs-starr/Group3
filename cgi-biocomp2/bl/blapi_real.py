#!/usr/bin/python3
"""
...Comment header goes here...
This is the business logic API
"""

# Add the bl sub-directory to the module path (for testing this routine)
# and the directory above to import the config file
import sys
sys.path.insert(0, "../db/")
sys.path.insert(0, "../")

import dbapi   # Import the database api

def getAllEntries():
    """
    This function is designed to be called by the frontend when a list of gene
    identifiers alone is required.
    
    Returns
    -------
    
    """    
    return(db_api.all_genebank())

def calc_exons(dna_seq, exon_locations, comp_strand, codon_start):
    """
    This function calculates the exons
    dna_seq - string
    exon_location - list of tuples
    comp_string - boolean
    codon_start - int
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
        temp_exon_store.append(ex_dna_seq[start:stop])
    working_exon_string = "".join(temp_exon_store)
    
    trim = (codon_start - 1)*3
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

def rest_enzyme_activity(dna_seq, rest_enzyme):
    import re
    """
    This function searches a genetic sequence for restriction enzyme activity
    Inputs
    genetic sequence in string format
    restriction enzyme - EcoR1, BamH1, BsuM1

    Returns
    -------
    Restriction enzyme activity, location
    """
    
    """
    EcoR1- forward G/AATTC, complementary CTTAA/G
    BAMH1 - forward G/GATCC, complementary CCTAG/G
    BsuM - forward CTCGAG, complementary CTCGAG
    """
    if rest_enzyme == "ecor1":
        print("ECOR1") #placeholder
    elif rest_enzyme == "bamh1":
        print("BAMH1") #placeholder
    elif rest_enzyme == "bsum":
        print("BSUM") #placeholder
    
    ##The following code identifies any RE site
    matched_loc = re.finditer("insert entry here variable", dna_seq) #Need to updat variable
    rest_enzyme_locations=[]
    for loc in matched_loc:
        rest_enzyme_locations.append(loc.span()) #tuples
    
    ## If no match found.
    if len(rest_enzyme_locations) == 0:
        print("There is no restriction enzyme match for", rest_enzyme_locations)
        
    return("NILL") ##To update

def runAllcodon_use():
    """
    This function is a one time run event.  To calculate all coding regions within the database
    and to calculate

    Returns
    -------
    Writes data into a text file. (name TBD)
    (Need to include - date/time of calculation)
    """
    return("into a text file")

def getAllcodon_use():
    """
    This function allows the web frontend to access the pre-calculated - all codon use file

    Returns
    -------
    None.
    """
    return("textfile.txt")
    
def frontend_input(request):
    """
    This function takes a x item tuple or list of format:
        (searchfield, data_type, codon_useage)
        searchfield - string object
        data_type - string object
        codon_useage - boolean (True/False)
    Returns
    -------
    
    """
    searchfield = request[0]
    data_type = request[1]
    codon_flag = request[2]
    
    
    temp_data=[]
    if data_type == "gene_id":
       temp_data = ["db_gene_id", "db_accession_code", "db_product", 
                    "db_location", "db_translation", "db_dna_seq"]
       """"temp_data = db_api.getgeneetries(searchfield)"""
    if data_type == "gen_acc":
        temp_data = "accessionn success"
        """"temp_data = db_api.sequence(searchfield)"""
    if data_type == "prot_prod":
        temp_data = "protein success"
           
    if data_type == "chro_loc":
        temp_data = ("tuple1", "tuple2")
          
    dna_seq = temp_data[5]
    exon_locations = ()
    comp_strand_flag =""
    coding_regions = ()
          
    if codon_flag == 1:
        codon_useage("dummy input")
    
    """
    Storage
    """
    return("variables to be returned")

def codon_useage(seq_input):
    from collections import Counter
    import re
    
    input_sequence =str(seq_input) 
    
    start_split = re.split("ATG", input_sequence, maxsplit=1)
    print(start_split[1])
    
    """
    May need additional work to work out exons/intron locations
    before this can be passed into raw_seq
    
    AJKGFJHKSADJHKJKHTJKHLATGAAABBBAAACCCAABCCCTGAADSFGKLJASDFJLKDSAFLKJ
    """
    
    raw_seq = (start_split[1])
    len_rs = len(raw_seq)
    exp_iter = int(len_rs/3)
    current_iter = 0
    seq_slice = slice(0,3)
    codon_list = ["ATG",]
    
    while current_iter != exp_iter:
        codon = raw_seq[seq_slice]
        print(codon)
        print(current_iter, exp_iter)
        codon_list.append(codon)
        if codon == "TAA" or codon == "TAG" or codon == "TGA":
            break
        raw_seq = raw_seq[3:]
        current_iter += 1
    
    print(codon_list)
    seq_total_codons = len(codon_list)
    
    k = Counter(codon_list)
        
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
    for (i,j) in master_codon_tuples:
        codon_use_data.update({i:(k[i],j, (k[i]/seq_total_codons)*100)})
    
    return(codon_use_data)

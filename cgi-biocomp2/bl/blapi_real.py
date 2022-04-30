# -*- coding: utf-8 -*-
#!/usr/bin/python3

"""
Created on Thu Apr 28 11:05:49 2022
This is the business logic API
@author: angli
"""

# Add the bl sub-directory to the module path (for testing this routine)
# and the directory above to import the config file
import sys
sys.path.insert(0, "../db/")
sys.path.insert(0, "../")

import dbapi   # Import the database api


def frontend_input(request):
    """
    This function takes a x item tuple or list of format:
        (searchfield, data_type, codon_useage)
        searchfield - string object
        data_type - string object
        rest_enzyme_flag - boolean (True/False)
        rest_enzyme_name 
    Returns
    -------
    #TBD
    """
    ##This is an assumption that the request is parsed as an iterable
    ##To check with Coco and update, could be passed as 3 variables instead
    searchfield = request[0]
    data_type = request[1]
    rest_enzyme_flag = request[2]
    if rest_enzyme_flag ==1:
        rest_enzyme_name = request[3]
    
    import db_API
    
    temp_data=[]
    
    #All data type variable names to confirm with Coco
    
    # Search by gene ID
    if data_type == "gene_id":
        temp_data = db_API.getgeneentries(searchfield)
       #["db_gene_id", "db_accession_code", "db_product", "db_location", "db_translation", "db_dna_seq"]  
       #also need exon_locations, comp_strand, codon_start
       dbgeneid =
       dbaccession =
       dbproduct =
       dblocation =
       dbtranslation = 
       db_dna_seq =
       exon_locations =
       comp_strand =
       codon_start =
    # Search by accession number
    if data_type == "gen_acc":
        temp_data = db_API.sequence(searchfield)
        ["gene_id","accession","protein id","location","translation","dna seq"]
        dbgeneid =
        dbaccession =
        dbproduct =
        dblocation =
        dbtranslation = 
        db_dna_seq =
        exon_locations =
        comp_strand =
        codon_start =
    # Search by protein product
    if data_type == "prot_prod":
        temp_data = db_API."variablename"(searchfield) #to complete once db_api updated
        dbgeneid =
        dbaccession =
        dbproduct =
        dblocation =
        dbtranslation = 
        db_dna_seq =
        exon_locations =
        comp_strand =
        codon_start =
    # Search by chromosome location
    if data_type == "chro_loc":
        temp_data = db_API."variablename"(searchfield) #to complete once db_api updated
        dbgeneid =
        dbaccession =
        dbproduct =
        dblocation =
        dbtranslation = 
        db_dna_seq =
        exon_locations =
        comp_strand =
        codon_start =   
    
    import calc_exons
    coding_string = calc_exons(db_dna_seq, exon_locations, comp_strand, codon_start)
    
    codon_data = {}
    import codon_useage
    codon_data = codon_useage(coding_dna_sequence)
    
    if renzyme == 1:
        import codon_useage
        rest_enzyme_activity(db_dna_seq, rest_enzyme_name, exon_locations)
    
    """
    Storage
    """
    #Putting everything togeather for output
    output_dictionary = {}
    ##Need to finalise format
    
    return("variables to be returned")

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def calc_exons(dna_seq, exon_locations, comp_strand, codon_start):
    """
    This function calculates the exons
    dna_seq - string
    exon_location - list of tuples
    comp_string - boolean
    codon_start - int
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
        temp_exon_store.append(dna_seq[start:(stop+1)])
    working_exon_string = "".join(temp_exon_store)
    
    #Addressing whether the sequence begins after the given co-ordinates. 
    #(Codon start is which nucleotide it start from with 1 being the co-ordinated given)
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

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def rest_enzyme_activity(dna_seq, rest_enzyme, input_exon_locations):
    import re
    """
    This function searches a genetic sequence for restriction enzyme activity
    Inputs
    genetic sequence in string format
    restriction enzyme - EcoR1, BamH1, BsuM1
    exon locations

    Returns
    -------
    Restriction enzyme activity in dictionary format {Enzyme name:[[(Start, Stop), Exon interference]]}
    Enzyme name - string
    Start - int
    Stop - int
    Exon_interference - bool - 0 is no mid coding
    """
    
    """
    EcoR1- forward G/AATTC, complementary CTTAA/G
    BAMH1 - forward G/GATCC, complementary CCTAG/G
    BsuM - forward CTCGAG, complementary CTCGAG
    All 3 are palindromes, therefore no need to specify/search on complementary strand.
    """
    enzyme_target_seq = []
    if rest_enzyme == "ecor1":
        enzyme_target_seq.append(("gaattc", "ecor1"))
    elif rest_enzyme == "bamh1":
        enzyme_target_seq.append(("ggatcc", "bamh1"))
    elif rest_enzyme == "bsum":
        enzyme_target_seq.append(("ctcgag", "bsum"))
    elif rest_enzyme == "all":
        enzyme_target_seq.append(("gaattc", "ecor1"))
        enzyme_target_seq.append(("ggatcc", "bamh1"))
        enzyme_target_seq.append(("ctcgag", "bsum"))
    
    enzyme_activity_dict = {} #storing enzyme name (key) against activity locations (list of tuples)
    for (e_seq, e_name) in enzyme_target_seq:
    ##The following code identifies any RE site
        matched_loc = re.finditer(e_seq, dna_seq)
        rest_enzyme_locations=[]
        for loc in matched_loc:
            rest_enzyme_locations.append(loc.span()) #tuples
    
        if len(rest_enzyme_locations) != 0:
            exon_locations = []
            for (start,stop) in input_exon_locations:
                if start[0] == "<" or start[0] == ">":
                    start = int(start[1:])
                else:
                    start = int(start)
                if stop[0] == "<" or stop[0] == ">":
                    stop = int(stop[1:])
                else:
                    stop = int(stop)
                exon_locations.append((start,stop))
                
            exon_start_boundary = min(exon_locations[0])
            exon_stop_boundary = max(exon_locations[1])         
                
        else:
            print("No match for:", e_name) #expand on if no match found
        enz_store = []    
        for (enzyme_start,enzyme_stop) in rest_enzyme_locations:
            if (enzyme_start < exon_start_boundary and enzyme_stop < exon_start_boundary):
                enz_store.append([(enzyme_start, enzyme_stop), 0])
            elif (enzyme_start > exon_stop_boundary and enzyme_stop > exon_start_boundary):
                enz_store.append([(enzyme_start, enzyme_stop), 0])
            else:
                enz_store.append([(enzyme_start, enzyme_stop), 1])
            enzyme_activity_dict[e_name] = enz_store
            ##Needs additional code to calculate whether there is any interference.  Planned
            ##as 
    return(enzyme_activity_dict) 
##To update - dictionary {restriction enzyme name (key),[activity sites, interference}
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def codon_useage(coding_dna_sequence):
    from collections import Counter
    import re
    
    input_sequence=coding_dna_sequence.upper()
    start_split = re.split("ATG", input_sequence, maxsplit=1) #this is probably redundant - review needed
    print(start_split[1])
    
    
    raw_seq = (start_split[1])
    len_rs = len(raw_seq)
    exp_iter = int(len_rs/3)
    current_iter = 0
    seq_slice = slice(0,3)
    
    #add ATG to the start - previous split method lopped off ATG, likely redundant if exon calculation done first
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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def getAllEntries():
    """
    This function is designed to be called by the frontend when a list of gene
    identifiers alone is required.
    
    Returns
    -------
    
    """    
    return(db_api.all_genebank())
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
    Text file
    """
    return("textfile.txt")

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
    This function takes a 3 or 4 item iterable object (tuple or list)
    with format: (searchfield, data_type, rest_enzyme_flag, rest_enzyme_name)
        searchfield - string object
        data_type - string object of options "gene_id", "gen_acc", "prot_prod" or "chro_loc"
        rest_enzyme_flag - boolean (True/False)
        rest_enzyme_name - string object
        
    Returns
    -------
    Dictionary with the following structure for 1 key:value pair:
        {Gene name:([Accession, Protein product name, Chromosomal location],[raw DNA sequence, Exon locations], alignment, codon_data, renzyme_output)}
        Gene name - string - dictionary key
        [Accession, Protein product name, Chromosomal location] - 3 item list of string objects
        [raw DNA sequence, Exon locations] - 2 item list, of string oject and list of tuples (of length 1->x)
        alignment - list of tuples with of (Amino acid, codon) pairing.  An error message is returned if alignment fails
        codon_data - dictionary of codon useage in format {codon:(number_of_occurances, amino_acid_letter, codon_use_perc)}
            codon - str
            number_of_occurances - int
            amino_acid_letter - str
            codon_use_perc - float
        renzyme_output - returns a tuple of 2 dictionaries
            renzyme_output[0] - {enzyme_name :[(start, stop), midcut]}
                enzyme_name - str
                (start,stop) - int (if activity in overall sequence), otherwise returns str ("No","activity")
                midcut - int (0 - no mid exon cut, 1 - mid exon cut, 2 - no activity on DNA strand) - specific to the (start,stop) site
            renzyme_output[1] - {enzyme_name:overall_midcut}
                overall_midcut - int (0 - activity present on DNA strand but suitable to isolate exons, 1 - enzyme not suitable to isolate exons)
        
    """
    ##Request is passed from frontend as 3 or 4 item iterable
    searchfield = request[0]
    data_type = request[1]
    rest_enzyme_flag = request[2]
    if rest_enzyme_flag ==1:
        rest_enzyme_name = request[3]
    
    import db_API
    
    temp_data=[]
    
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
        #["gene_id","accession","protein id","location","translation","dna seq"]
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
        temp_data = db_API."functionname"(searchfield) #to complete once db_api updated
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
        temp_data = db_API."functionname"(searchfield) #to complete once db_api updated
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
    codon_data = codon_useage(coding_string)
    
    import prot_gene_alignment
    alignment = prot_gene_alignment(dbtranslation, coding_string)
    
    if rest_enzyme_flag == 1:
        import codon_useage
        renzyme_output = rest_enzyme_activity(db_dna_seq, rest_enzyme_name, exon_locations)
    else:
        renzyme_output = None
    
    #Combining output
    output_dictionary = {dbgeneid:([dbaccession, dbproduct, dblocation],[db_dna_seq, exon_locations], alignment, codon_data, renzyme_output)}
    
    return(output_dictionary)

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

    """
    This function searches a genetic sequence for restriction enzyme activity and identifieds enzymes
    with activity outside of the exon start/stop region.
    
    Inputs
    -------
    genetic sequence in string format
    restriction enzyme - str - EcoR1, BamH1, BsuM1, all
    input_exon_locations - list of tuples

    Returns
    -------
    Returns a tuple of 2 dictionaries
       enzyme_activity_dict - {enzyme_name :[(start, stop), midcut]}
           enzyme_name - str
           (start,stop) - int (if activity in overall sequence), otherwise returns str ("No","activity")
           midcut - int (0 - no mid exon cut, 1 - mid exon cut, 2 - no activity on DNA strand) - specific to the (start,stop) site
       interference_output - {enzyme_name:overall_midcut}
           overall_midcut - int (0 - activity present on DNA strand but suitable to isolate exons, 1 - enzyme not suitable to isolate exons)
    """
    
    """
    Additional notes
    EcoR1- forward GAATTC, complementary CTTAA/G
    BAMH1 - forward GGATCC, complementary CCTAG/G
    BsuM - forward CTCGAG, complementary CTCGAG
    Most restriction enzymes are palindromes, therefore no need to specify/search on complementary strand based on entry requirements.
    Modifications to code may be required for non palindromic REs
    """
    import re
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
        else:
            enz_store=[(("No","activity"),2)]
            enzyme_activity_dict[e_name] = enz_store #expand on if no match found
    
    interference_output = []
    for (e_seq, e_name) in enzyme_target_seq:
        Exon_interefence = 0
        
        for i in range(0,len(enzyme_activity_dict.get(e_name))):
            Exon_interefence += int(enzyme_activity_dict.get(e_name)[i][1])
        if Exon_interefence == 0:
            interference_output.append((e_name, 0))
        else:
            interference_output.append((e_name, 1))
    return(enzyme_activity_dict, interference_output) 

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
    import getgetAllEntries
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
    entry_list = getgetAllEntries()
    
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

def getAllcodon_use():
    """
    This function allows the web frontend to access the pre-calculated - all codon use file

    Returns
    -------
    Total codon use textfile
    Format Codon, amino_acid, count, percentage
    Separated by " "
    """
    return(total_codon_use.txt)

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def prot_gene_alignment(protein_seq, gene_sequence):
    """
    This function takes in a protein sequence and genomic sequence and aligns them.
    If there is mismatch in length, the alignment does not occur.

    Parameters
    ----------
    protein_seq : str
    gene_sequence : str

    Returns
    -------
    List of tuples of format (amino acid, codon)
    """
    alignment_list =[]
    if len(protein_seq) != (len(gene_sequence)/3):
        alignment_list.append("ERROR: The sequences cannot be aligned")
    else:
        for i in range(0,len(protein_seq)):
            alignment_list.append((protein_seq[i], gene_sequence[i*3:(i*3)+3]))
    return(alignment_list)

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def rest_enzyme_list:
    """
    This function passes a list of available restriction enzymes to the frontend.
    """
    enzyme_list = ["ecor1", "bamh1", "bsum", "all"]
return(enzyme_list)

#!/usr/bin/python3.7.9
# -*- coding: utf-8 -*-

"""
Created on Thu Apr 28 11:05:49 2022
This is the business logic API, containing 4 of the 5 business layer functions that may be required by the front end.
@author: angli
"""

# Add the bl sub-directory to the module path (for testing this routine)
# and the directory above to import the config file
import sys
sys.path.insert(0, "../db/")
sys.path.insert(0, "../")
sys.path.insert(0, "../bl/")
import dbapi   # Import the database api

def frontend_input(request):
    """
    This function takes a 3 or 4 item iterable object (tuple or list), using the information to perform a search and carry out necessary processing to present data back to the front end.
    
    The input should be of format: (searchfield, data_type, rest_enzyme_flag, rest_enzyme_name)
        searchfield - string object
        data_type - string object of options "gene_id", "gen_acc", "prot_prod" or "chro_loc" - corresponding to the type of search data
        rest_enzyme_flag - boolean (True/False)
        rest_enzyme_name - string object (Optional)

    Returns
    -------
    Dictionary with the following structure for key:value pair:
        {index:([Gene name, Accession, Protein product name, Chromosomal location],[raw DNA sequence, Exon locations], alignment, codon_data, renzyme_output)} 
        
        index 0->x (number of entries returned, this is also the dictionary key) 
        Gene name - string - dictionary key
        [Accession, Protein product name, Chromosomal location] - 3 item list of string objects
        [raw DNA sequence, Exon locations] - 2 item list, of string object and list of tuples ()
        alignment - list of tuples with of (Amino acid, codon) pairing.  An error message is returned if alignment fails
        codon_data - dictionary of codon useage in format {codon:(number_of_occurances, amino_acid_letter, codon_use_perc)}
            codon - str - three letter codon
            number_of_occurances - int
            amino_acid_letter - str - single letter amino acid code
            codon_use_perc - float - percentage use of the codon within the sequence
        renzyme_output - returns a tuple of 2 dictionaries
            renzyme_output[0] - {enzyme_name :[(start, stop), midcut]}
                enzyme_name - str
                (start,stop) - int (if activity in overall sequence), formatted in format of the range function (i.e. the true end site is stop-1)
                    If no activity found for the restriction enzyme, the module returns the tuple ("No","activity")
                midcut - int (0 - no mid exon cut, 1 - mid exon cut, 2 - no activity on DNA strand) - specific to the (start,stop) site
            renzyme_output[1] - {enzyme_name:overall_midcut}
                overall_midcut - int (0 - activity present on DNA strand but suitable to isolate exons, 1 - enzyme not suitable to isolate exons)
    """

    # Request is passed from frontend as 3 or 4 item iterable
    searchfield = request[0]
    data_type = request[1]
    rest_enzyme_flag = request[2]
    if rest_enzyme_flag == 1:
        rest_enzyme_name = request[3]

    import db_api

    temp_data = []

    # Search by gene ID
    if data_type == "gene_id":

        temp_data = db_API.getgeneentries(searchfield)
       # From db_api_dummy:
       # db_api.getgeneentries
       # ["db_gene_id", "db_accession_code", "db_product", "db_location", "db_translation", "db_dna_seq"]
       # Assumption that this will return a list of lists
                       
       """
       This section of code should go on to include information stored from the cds sections associated with each gene entry
       With appropriate appending to each list entry within "temp_data".
       An idealised solution will be to call an additional db_api function to :
           1. Recieve exon_locations as a list of tuples of format [("string", "string"), ("string","string")...]
           2. Comp_strand as a boolean object
           3. Codon_start as integer
       """
    # Search by accession number
    if data_type == "gen_acc":
        temp_data = db_api.sequence(searchfield)
        #["db_gene_id", "db_accession_code", "db_product", "db_location", "db_translation", "db_dna_seq"]
        """
        As per description for the search by gene_id, this would have duplicated code to ensure that the appropriate db_api functions were called
        and appropriately joined.
        """

    # Search by protein product
    if data_type == "prot_prod":
        # to complete once db_api updated
        temp_data = db_api."functionname"(searchfield)
        """
        As per description for the search by gene_id, this would have duplicated code to ensure that the appropriate db_api functions were called
        and appropriately joined.
        """
    # Search by chromosome location
    if data_type == "chro_loc":
        # to complete once db_api updated
        temp_data = db_api."functionname"(searchfield)
        """
        As per description for the search by gene_id, this would have duplicated code to ensure that the appropriate db_api functions were called
        and appropriately joined.
        """
        
    #Iteration through each search resule
    for each_entry in temp_data:
        dbgeneid = each_entry[0]
        dbaccession = each_entry[1]
        dbproduct = each_entry[2]
        dblocation = each_entry[3]
        dbtranslation = each_entry[4]
        db_dna_seq = each_entry[5]
        exon_locations = each_entry[6]
        comp_strand = each_entry[7]
        codon_start = each_entry[8]
        
        #Calculate the coding sequence
        import calc_exons as ce
        coding_string = ce.calc_exons(db_dna_seq, exon_locations, comp_strand, codon_start)

        #Calculation of codon useage
        codon_data = {}
        import codon_useage as cu
        codon_data = cu.codon_useage(coding_string)

        import prot_gene_alignment as pga
        alignment = pga.prot_gene_alignment(dbtranslation, coding_string)

        if rest_enzyme_flag == 1:
            import rest_enzyme_activity as r_e_a
            renzyme_output = r_e_a.rest_enzyme_activity(
                db_dna_seq, rest_enzyme_name, exon_locations)
        else:
            renzyme_output = None

    # Combining output
        output_dictionary = {temp_data.index(each_entry): ([dbgeneid, dbaccession, dbproduct, dblocation], [
            db_dna_seq, exon_locations], alignment, codon_data, renzyme_output)}

    return(output_dictionary)

# ------------------------------------------

def getAllEntries():
    """
    This function is designed to be called by the frontend when a list of gene
    identifiers alone is required.

    Returns
    -------
    From db_API: [('accession', 'gene', 'protein_id', 'sequence', 'aa_seq', 'chromosomal location')]
    """
    import db_api
    return(db_api.all_genebank())

# ------------------------------------------

def getAllcodon_use():
    """
    This function allows the web frontend to access the pre-calculated - all codon use file

    Returns
    -------
    Total codon use textfile
    Format Codon, amino_acid, count, percentage
    Separated by " "
    """
    file_location = total_codon_use.txt
    return(file_location)

# -----------------------------------------

def rest_enzyme_list():
    """
    This function passes a list of available restriction enzymes to the frontend.
    """
    enzyme_list = ["ecor1", "bamh1", "bsum", "all"]

return(enzyme_list)

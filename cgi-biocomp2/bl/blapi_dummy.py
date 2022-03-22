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
import config  # Import configuration information (if needed)

def getAllEntries():
    """
    ...Function comment header goes here...
    This is a very simple function that just calls the database API to do the SQL to 
    obtain the full list of entries. It doesn't need to do anything else.
    """
    return(dbapi.getAllEntries()) 
    """This directly accesses the DB API, therefore need to discuss with Janani as to formatting for Coco."""
    """From the feedback session - ?cached vs denovo search from clicking from the list - if caching needsxz to happen, may need to store on DB"""
    
def frontend_input():
    """
    This function processes inputs from the HTML layer, performs error checking against standard formats
    Input:
        - search_field - the user inputted search parameter [string]
        - search_data_type - the type of data of the seearch parameter - i.e. Gene ID, ACCESSION number etc.. [string]
                -(gene_id, prot_prod, gen_acc, chro_loc) [string]
    Processing:
        This uses the search field 
    Returns:
    
    Stored variables
    """
    return()

def db_input():
    """
    This function uses stored variables from frontend_input() and presents the request to the DB layer.
    Input
    
    Returns:
    """
    return("database API search function")
    
def codon_useage():
    """
    This function calculates the codon useage of a given gene
    Inputs:
        - DNA sequence from database [string]
        - Coding region markers[list]
    Returns:
        - Codon_use_data [dictionary{codon(key), count, amino acid, percentage use}]
    """
    return(codon_use_data)
    
def renzyme_acitivity():
    """
    This function takes 
    Will be called in renzyme activity searches
        - EcoRI, BamHI, BsuMI - via dropdown menu
    """
    
    
Output to Weblayer:
List: ('gene_id'[string],protein product name [string], accession [string], chromosome location [tuple], complete DNA seq [string], coding region markers [list of tuples (start:stop), amino acid sequence[string], codon useage [dictionary{codon(key), count, amino acid, percentage use}])

Restriction enzyme:
List (restriction enzyume name [string], cut5'[bool], cut3'[bool], middle[bool], cutsite (list of tuples (start:stop) - may be NULL/NONE)

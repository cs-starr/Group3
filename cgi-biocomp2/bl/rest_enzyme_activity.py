# -*- coding: utf-8 -*-
#!/usr/bin/python3

"""
Created on Sun 01 May 2022 14:36
This is the script for calculation of restriction enzyme activity
@author: angli
"""

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

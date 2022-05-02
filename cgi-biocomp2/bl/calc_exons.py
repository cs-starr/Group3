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
            stop = int(stop)-1
        temp_exon_store.append(dna_seq[(start-1):(stop)])
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

"""
#Test code:
ds = "aggaaggcatgcggtgggccctactggtgcttctagctttcctgtctcctggtgagtacgctgcctacagagaggctcacaggttgggttttgttttgttttcttcttgaaaggggtgccatacaaaggaatacctcattgtattttgtgttgttcccattgcagccagtcagaaatcttccaacttggaagggagaacgaagtcagtcaccaggcagactgggtcatctgctgaaatcacttgcgatcttactgtaacaaataccttctacatccactggtacctacaccaggaggggaaggccccacagcgtcttctgtactatgacgtctccaccgcaagggatgtgttggaatcaggactcagtccaggaaagtattatactcatacacccaggaggtggagctggatattgagactgcaaaatctaattgaaaatgattctggggtctattactgtgccacctggcggacgaattattataagaaactctttggcagtggaacaacacttgttgtcacaggtaagtatcggaagaatacaacatttccaaggtaatagagggaaggcaggaaatgattaaactggaataatgt"
el =  [("9","51"),("166",">525")]
cs = 0
csod = 1

ds2 = "ctgcttgatggcaaagaaataaagcgactgaatgttcagtggctccgagcacacctgggcatcgtgtcccaggagcccatcctgtttgactgcagcattgctgagaacattgcctatggagacaacagccgggtggtgtcacaggaagagatcgtgagggcagcaaaggaggccaacatacatgccttcatcgagtcactgcctaat"
el2 =  [("<1",">207")]
cs2 = 0
csod2 = 1
print(calc_exons(ds,el,cs,csod))
print(len(calc_exons(ds,el,cs,csod)))
print("space-------")
print(calc_exons(ds2,el2,cs2,csod2))
print(len(calc_exons(ds2,el2,cs2,csod2)))
"""

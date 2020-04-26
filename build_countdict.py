# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 01:24:15 2020

@author: Hooooaaanng
Building a dictionary storing the count of bacteria in bacteria dictionary for the specified kmer length
"""

import time
from Bio import SeqIO
import pickle

bacteria_dict = {'Acinetobacter baumannii':'abauman','Bacteroides fragilis':'bfragil','Escherichia coli':'ecolik12','Enterobacter cloacae':'entcloac', 'E. faecium':'entfaec', 'Klebsiella pneumoniae':'klebpneu','Proteus mirabilis':'pmirabil','Pseudomonas aeruginosia':'pseudaer','Staph. epidermidis':'stapepid','Staph. aureus':'staphaur','Staphylococcus saprophyticus':'stapsapr','Streptococcus agalactiae':'strepaga','Streptococcus pneumoniae':'strepneu'}

# Import all bacteria file and set up global variables with bacteria files as name
for bacteria in bacteria_dict.values():
    record = 0
    globals()[bacteria] = []
    for seq_record in SeqIO.parse(bacteria + ".fna", "fasta"):
        #print(seq_record)
        globals()[bacteria + '_' + str(record)] = seq_record
        globals()[bacteria].append(globals()[bacteria + '_' + str(record)])
        #print(globals()[bacteria])
        record += 1

# Construct the count matrix 
def kmer_count(k, bacteria_dict):
    start_time = time.time()
    bacteria_index = 0
    count_dict = {}
    for bacteria in bacteria_dict.values():
        count_dict[str(bacteria)] = {}
        sequence = ''
        print("Starting for ", str(bacteria))
        # concatenate all the records found for each bacterial sequence
        for records in globals()[bacteria]:
            sequence = sequence + str(records.seq)

        for i in range(len(sequence) - k + 1):

            kmer =  sequence[i:i+k]
            if kmer in count_dict[str(bacteria)]:
                count_dict[str(bacteria)][kmer] += 1
            else:
                count_dict[str(bacteria)][kmer] = 1       
        
        print("Time taken for ", str(bacteria), ' is ', time.time() - start_time)
        bacteria_index += 1
    return(count_dict)
    
# main function to build and dump the dictionary into a pickle file
if __name__ == '__main__':
    start_time = time.time()
    print("Building kmer set time takes ", time.time()-start_time)
    count_matrix = kmer_count(10,bacteria_dict)
    pickle.dump(count_matrix, open('10mer_all_count_dict','wb'))

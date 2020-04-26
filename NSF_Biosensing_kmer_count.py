from Bio.Seq import Seq
import time
from Bio.Alphabet import IUPAC
from Bio import SeqIO
import itertools
import random
import numpy as np
import numpy.linalg as linalg
import scipy.linalg as la
import pickle
import pandas as pd
import matplotlib.pyplot as plt
from kmer_classes_functions import *
from copy import copy
from tqdm import tqdm_notebook
import multiprocessing

bacteria_dict = {'Acinetobacter baumannii':'abauman','Bacteroides fragilis':'bfragil','Escherichia coli':'ecolik12','Enterobacter cloacae':'entcloac', 'E. faecium':'entfaec', 'Klebsiella pneumoniae':'klebpneu','Proteus mirabilis':'pmirabil','Pseudomonas aeruginosia':'pseudaer','Staph. epidermidis':'stapepid','Staph. aureus':'staphaur','Staphylococcus saprophyticus':'stapsapr','Streptococcus agalactiae':'strepaga','Streptococcus pneumoniae':'strepneu'}
kmer_num = 3

def clear_all():
    """Clears all the variables from the workspace of the spyder application."""
    gl = globals().copy()
    for var in gl:
        if var[0] == '_': continue
        if 'func' in str(globals()[var]): continue
        if 'module' in str(globals()[var]): continue

        del globals()[var]
#Enumerate all possible permutations of 8-mer
with open('10mer_all', 'rb') as filehandle:
    # read the data as binary data stream
    kmer_enumerate = pickle.load(filehandle)


    
def multiprocess_ransearch(seed_list, ransearch_fitness_best, ransearch_fitness_average, ransearch_probe_set, ransearch_count, ransearch_L2pairwise, ransearch_L2column, ransearch_seed, seed_index, alpha, offset, *args):
    print("Starting a random search process")
    start_time = time.time()
    seed = seed_list[seed_index]
    [rand_cycle_track, probe_set, count, pop_fitness_best, pop_fitness_average, L2_pairwise, L2_column]  = random_search(100, kmer_enumerate, kmer_num, bacteria_dict, seed, 30, alpha, offset)
    ransearch_fitness_average.append(pop_fitness_average)
    ransearch_fitness_best.append(pop_fitness_best)
    ransearch_probe_set.append(probe_set)
    ransearch_count.append(count)
    ransearch_seed.append(seed)
    ransearch_L2pairwise.append(L2_pairwise)
    ransearch_L2column.append(L2_column)
    print("Ending a random search process. Time taken for this process is ", time.time() - start_time)

def multiprocess_ga(seed_list, ga_fitness, ga_ave_pop_fitness, ga_probe_set, ga_count, ga_L2pairwise, ga_L2column, ga_seed, seed_index, alpha, offset, *args):
    print("Starting a GA process")
    start_time = time.time()
    seed = seed_list[seed_index]
    ga_seed.append(seed)
    cycle_track, fitness_score, pop_ave_fitness, probe_set, count, exe_time, L2pairwise, L2column = gene_alg(100, 30, kmer_num, kmer_enumerate, bacteria_dict, seed, alpha, offset)
    ga_fitness.append(fitness_score)
    ga_ave_pop_fitness.append(pop_ave_fitness)
    ga_probe_set.append(probe_set)
    ga_count.append(count)
    ga_L2pairwise.append(L2pairwise)
    ga_L2column.append(L2column)
    print("Ending a GA process. Time taken for this process is ", time.time() - start_time)
    
#if __name__ == '__main__':
#    # Dictionary with bacteria full name as keys and variable name/fna file name as values
#    bacteria_dict = {'Acinetobacter baumannii':'abauman','Bacteroides fragilis':'bfragil','Escherichia coli':'ecolik12','Enterobacter cloacae':'entcloac', 'E. faecium':'entfaec', 'Klebsiella pneumoniae':'klebpneu','Proteus mirabilis':'pmirabil','Pseudomonas aeruginosia':'pseudaer','Staph. epidermidis':'stapepid','Staph. aureus':'staphaur','Staphylococcus saprophyticus':'stapsapr','Streptococcus agalactiae':'strepaga','Streptococcus pneumoniae':'strepneu'}
#
#    # look at setting a seed for an rng 
#    start_time = time.time()
#    # Initialize variables
#    kmer_num = 3
#        
#    #Enumerate all possible permutations of 8-mer
#    with open('10mer_all', 'rb') as filehandle:
#        # read the data as binary data stream
#        kmer_enumerate = pickle.load(filehandle)
#
#
#    seed_list = pickle.load(open('ransearch_seed','rb'))
#    start_time = time.time()
#    seed_index = [0,1,2,3,4,5,6,7,8,9]
#    manager = multiprocessing.Manager()
#    
#    ransearch_fitness_best = manager.list([]*len(seed_index))
#    ransearch_fitness_average = manager.list([]*len(seed_index))
#    ransearch_probe_set = manager.list([]*len(seed_index))
#    ransearch_count = manager.list([]*len(seed_index))
#    ransearch_seed = manager.list([]*len(seed_index))
#    ransearch_L2pairwise = manager.list([]*len(seed_index))
#    ransearch_L2column = manager.list([]*len(seed_index))
#    alpha = 1000
#    offset = -0.0000001 
#    p = []
#    
#    for index in range(len(seed_index)):
#        p.append(multiprocessing.Process(target=multiprocess_ransearch, args = (seed_list, ransearch_fitness_best, ransearch_fitness_average, ransearch_probe_set, ransearch_count, ransearch_L2pairwise, ransearch_L2column, ransearch_seed, seed_index[index], alpha, offset)))
#        p[index].start()
#    
#    for index in range(len(seed_index)):
#        p[index].join()
#    
#
#    print("random search, 10 seeds, 100 cycles", time.time() - start_time, " sec")
# 
#    pickle.dump(list(ransearch_seed), open('10mer_ransearch_seed','wb'))
#    pickle.dump(list(ransearch_fitness_best), open('10mer_ransearch_fitness_best','wb'))  
#    pickle.dump(list(ransearch_fitness_average), open('10mer_ransearch_fitness_average','wb'))
#    pickle.dump(list(ransearch_probe_set), open('10mer_ransearch_probe_set','wb'))
#    pickle.dump(list(ransearch_L2pairwise), open('10mer_ransearch_L2pairwise','wb'))
#    pickle.dump(list(ransearch_L2column), open('10mer_ransearch_L2column','wb'))
#    pickle.dump(list(ransearch_count), open('10mer_ransearch_count','wb'))
#   
#    print("Total time taken, random search: ", time.time() - start_time)
#    
#    
#    clear_all()
#    # Dictionary with bacteria full name as keys and variable name/fna file name as values
#    bacteria_dict = {'Acinetobacter baumannii':'abauman','Bacteroides fragilis':'bfragil','Escherichia coli':'ecolik12','Enterobacter cloacae':'entcloac', 'E. faecium':'entfaec', 'Klebsiella pneumoniae':'klebpneu','Proteus mirabilis':'pmirabil','Pseudomonas aeruginosia':'pseudaer','Staph. epidermidis':'stapepid','Staph. aureus':'staphaur','Staphylococcus saprophyticus':'stapsapr','Streptococcus agalactiae':'strepaga','Streptococcus pneumoniae':'strepneu'}
#    
#    # look at setting a seed for an rng 
#    start_time = time.time()
#    # Initialize variables
#    kmer_num = 3
#            
#     #Enumerate all possible permutations of 8-mer
#    with open('10mer_all', 'rb') as filehandle:
#        # read the data as binary data stream
#        kmer_enumerate = pickle.load(filehandle)
#
#
#    seed_list = pickle.load(open('ransearch_seed','rb'))
#    start_time = time.time()
#    seed_index = [0,1,2,3,4,5,6,7,8,9]
#    manager = multiprocessing.Manager()
#    ga_fitness = manager.list([]*len(seed_index))
#    ga_ave_pop_fitness = manager.list([]*len(seed_index))
#    ga_probe_set = manager.list([]*len(seed_index))
#    ga_count = manager.list([]*len(seed_index))
#    ga_L2pairwise = manager.list([]*len(seed_index))
#    ga_L2column = manager.list([]*len(seed_index))
#    ga_seed = manager.list([]*len(seed_index))
#    alpha = 1000
#    offset = -0.0000001 
#    p = []
#    
#    for index in range(len(seed_index)):
#        p.append(multiprocessing.Process(target=multiprocess_ga, args = (seed_list, ga_fitness, ga_ave_pop_fitness, ga_probe_set, ga_count, ga_L2pairwise, ga_L2column, ga_seed, seed_index[index], alpha, offset)))
#        p[index].start()
#    
#    for index in range(len(seed_index)):
#        p[index].join()
#    
#    print("GA search, 10 seeds, 100 cycles", time.time() - start_time, " sec")
#    
#    pickle.dump(list(ga_seed), open('10mer_ga_seed','wb')) 
#    pickle.dump(list(ga_fitness), open('10mer_ga_fitness','wb')) 
#    pickle.dump(list(ga_probe_set), open('10mer_ga_probe_set','wb'))   
#    pickle.dump(list(ga_count), open('10mer_ga_count','wb'))
#    pickle.dump(list(ga_ave_pop_fitness), open('10mer_ga_ave_pop_fitness','wb'))
#    pickle.dump(list(ga_L2pairwise), open('10mer_ga_L2pairwise','wb'))
#    pickle.dump(list(ga_L2column), open('10mer_ga_L2column_alpha1000','wb'))
#    
#    print("Total time taken, GA: ", time.time() - start_time)

# read the data as binary data stream
#rand_count = pickle.load(open('8mer_ransearch_count6','rb'))
rand_fitness_pairwise = np.array(pickle.load(open('12mer_ransearch_L2pairwise6','rb')) + pickle.load(open('12mer_ransearch_L2pairwise6_1','rb')) + pickle.load(open('12mer_ransearch_L2pairwise6_2','rb')))
rand_fitness_column = np.array(pickle.load(open('12mer_ransearch_L2column6','rb')) + pickle.load(open('12mer_ransearch_L2pairwise6_1','rb')) + pickle.load(open('12mer_ransearch_L2pairwise6_2','rb')))
rand_fitness_best = np.array(pickle.load(open('12mer_ransearch_fitness_best6','rb')) + pickle.load(open('12mer_ransearch_L2pairwise6_1','rb')) + pickle.load(open('12mer_ransearch_L2pairwise6_2','rb')))
rand_fitness_average = np.array(pickle.load(open('12mer_ransearch_fitness_average6','rb')) + pickle.load(open('12mer_ransearch_L2pairwise6_1','rb')) + pickle.load(open('12mer_ransearch_L2pairwise6_2','rb')))
#rand_count = pickle.load(open('12mer_ransearch_count6','rb')) + pickle.load(open('12mer_ransearch_L2pairwise6_1','rb')) + pickle.load(open('12mer_ransearch_L2pairwise6_2','rb'))
rand_probe = pickle.load(open('12mer_ransearch_probe_set6','rb')) + pickle.load(open('12mer_ransearch_L2pairwise6_1','rb')) +    pickle.load(open('12mer_ransearch_L2pairwise6_2','rb'))

ga_count = pickle.load(open('12mer_ga_count6','rb')) + pickle.load(open('12mer_ga_count6_2','rb')) +  pickle.load(open('12mer_ga_count6_3','rb')) + pickle.load(open('12mer_ga_count6_4','rb'))
ga_fitness = np.array(pickle.load(open('12mer_ga_fitness6','rb')) + pickle.load(open('12mer_ga_fitness6_2','rb')) + pickle.load(open('12mer_ga_fitness6_3','rb')) + pickle.load(open('12mer_ga_fitness6_4','rb')))
ga_probe = pickle.load(open('12mer_ga_probe_set6','rb')) + pickle.load(open('12mer_ga_probe_set6_2','rb')) + pickle.load(open('12mer_ga_probe_set6_3','rb')) + pickle.load(open('12mer_ga_probe_set6_4','rb'))
ga_ave_pop_fitness = np.array(pickle.load(open('12mer_ga_ave_pop_fitness6','rb')) + pickle.load(open('12mer_ga_ave_pop_fitness6_2','rb')) + pickle.load(open('12mer_ga_ave_pop_fitness6_3','rb')) + pickle.load(open('12mer_ga_ave_pop_fitness6_4','rb')) )
ga_fitness_pairwise = np.array(pickle.load(open('12mer_ga_L2pairwise6','rb')) + pickle.load(open('12mer_ga_L2pairwise6_2','rb')) + pickle.load(open('12mer_ga_L2pairwise6_3','rb')) + pickle.load(open('12mer_ga_L2pairwise6_4','rb')))
ga_fitness_column = np.array(pickle.load(open('12mer_ga_L2colum6','rb')) + pickle.load(open('12mer_ga_L2colum6_2','rb')) + pickle.load(open('12mer_ga_L2colum6_3','rb')) + pickle.load(open('12mer_ga_L2colum6_4','rb')))

cycle = []
for i in range(100):
    cycle.append(i+1)
#
rand_fitness_best_mean = np.mean(rand_fitness_best, 0)
rand_fitness_best_stdev = np.std(rand_fitness_best, 0)
rand_fitness_ave_mean = np.mean(rand_fitness_average, 0)
rand_fitness_ave_stdev = np.std(rand_fitness_average, 0)
rand_fitness_pairwise_mean = np.mean(rand_fitness_pairwise, 0)
rand_fitness_pairwise_stdev = np.std(rand_fitness_pairwise, 0)
rand_fitness_column_mean = np.mean(rand_fitness_column, 0)
rand_fitness_column_stdev = np.std(rand_fitness_column, 0)

ga_fitness_mean = np.mean(ga_fitness,0)
ga_fitness_stdev = np.std(ga_fitness,0)
ga_ave_pop_fitness_mean = np.mean(ga_ave_pop_fitness,0)
ga_ave_pop_fitness_stdev = np.std(ga_ave_pop_fitness,0)
ga_fitness_pairwise_mean = np.mean(ga_fitness_pairwise, 0)
ga_fitness_pairwise_stdev = np.std(ga_fitness_pairwise, 0)
ga_fitness_column_mean = np.mean(ga_fitness_column, 0)
ga_fitness_column_stdev = np.std(ga_fitness_column, 0)

plt.figure(figsize = (25,10))
plt.errorbar(cycle,list(rand_fitness_best_mean),list(rand_fitness_best_stdev), linestyle = 'None', marker = 'o',color = 'blue',label = 'Random Search Best Solution')
plt.errorbar(cycle,list(ga_fitness_mean),list(ga_fitness_stdev), linestyle = 'None', marker = 'o',color = 'red', label = "Genetic Algorithm Best Solution")

plt.xlabel('cycle', fontsize = 20)
plt.ylabel('fitness', fontsize = 20)
plt.legend(fontsize = 20, loc='upper left')

plt.figure(figsize = (25,10))
plt.errorbar(cycle,list(rand_fitness_ave_mean),list(rand_fitness_ave_stdev), linestyle = 'None', marker = 'o',color = 'grey',label = 'Random Search Population Average Fitness')
plt.errorbar(cycle,list(ga_ave_pop_fitness_mean),list(ga_ave_pop_fitness_stdev), linestyle = 'None', marker = 'o',color = 'red',label = 'Genetic Algorithm Population Average Fitness')
plt.xlabel('cycle', fontsize = 20)
plt.ylabel('fitness', fontsize = 20)
plt.legend(fontsize = 20, loc='upper left')


plt.figure(figsize = (25,10))
plt.errorbar(cycle,list(rand_fitness_pairwise_mean),list(rand_fitness_pairwise_stdev), linestyle = 'None', marker = 'o',color = 'blue',label = 'Random Search pairwise component')
plt.errorbar(cycle,list(ga_fitness_pairwise_mean),list(ga_fitness_pairwise_stdev), linestyle = 'None', marker = 'o',color = 'red', label = "Genetic Algorithm pairwise component")

plt.xlabel('cycle', fontsize = 20)
plt.ylabel('fitness', fontsize = 20)
plt.legend(fontsize = 20, loc='upper left')

plt.figure(figsize = (25,10))
plt.errorbar(cycle,list(rand_fitness_column_mean),list(rand_fitness_column_stdev), linestyle = 'None', marker = 'o',color = 'blue',label = 'Random Search column component')
plt.errorbar(cycle,list(ga_fitness_column_mean),list(ga_fitness_column_stdev), linestyle = 'None', marker = 'o',color = 'red', label = "Genetic Algorithm column component")

plt.xlabel('cycle', fontsize = 20)
plt.ylabel('fitness', fontsize = 20)
plt.legend(fontsize = 20, loc='upper left')





        


    
    


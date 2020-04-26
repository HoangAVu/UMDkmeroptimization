# -*- coding: utf-8 -*-
"""
Created on Sat Apr  4 16:48:26 2020

@author: Hooooaaanng
"""

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
from copy import copy
from tqdm import tqdm_notebook

# Dictionary with bacteria full name as keys and variable name/fna file name as values
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

# Function to generate a set containing all possible permutation of a kmer given k, as DNA format
def kmer_perm(k):
    kmer_set = set([''.join(i) for i in itertools.product('ATCG', repeat = k)])
    return kmer_set

# Random search function for maximizing fitness score
def random_search(cycles, kmer_enumerate, kmer_num, bacteria_dict, seed, pop_size, alpha, offset):
    random.seed(seed)
    
    bacteria_num = len(bacteria_dict)
    cycle_track = []
    min_L2norm_track = []
    fitness_count = []
    total_records_length = np.zeros((1,bacteria_num))
    pop_best_fitness = []
    pop_average_fitness = []
    L2min_pairwise_diff_track = []
    L2min_column_track = []
    

    for cycle_num in range(cycles):
        cycle_track.append(cycle_num)
        
        pop = []
        pop_count = []
        pop_fitness_cycle = []
        pairwise_cycle = []
        min_column_cycle = []
        
        for i in range(pop_size):
            # Choose kmer_num random 8-mer from the set of all possible 8-mer
            kmer_set = random.sample(tuple(kmer_enumerate),kmer_num)
            pop.append(kmer_set)
            # Build count matrix for the given set of probe
            Count_Matrix = kmer_count(kmer_num, bacteria_num, bacteria_dict, kmer_set)
            pop_count.append(Count_Matrix)
            # Calculate L2 norms for pairwise differences among bacteria's column vector
                # Note: bacteria's column vector is the length of the number of probes
          
            fitness = UMD_fitness_score(Count_Matrix, bacteria_dict)
            pop_fitness_cycle.append(fitness.L2norm_min_pairwise_diff() + alpha*np.log(fitness.L2norm_min_column() - offset))
            min_L2norm_track.append(fitness.L2norm_min_column())
            pairwise_cycle.append(fitness.L2norm_min_pairwise_diff())
            min_column_cycle.append(fitness.L2norm_min_column())
        
        
            if max(min_L2norm_track) == min_L2norm_track[-1]:
                best_probe_set = kmer_set
                if cycle_num == 0:
                    fitness_count.append(fitness.L2norm_min_pairwise_diff() + alpha*np.log(fitness.L2norm_min_column() - offset))
                else:
                    fitness_count.append(fitness_count[-1])
                best_count = Count_Matrix
            else:
                fitness_count.append(fitness.L2norm_min_pairwise_diff() + alpha*np.log(fitness.L2norm_min_column() - offset))
            
        pop_best_fitness.append(max(pop_fitness_cycle))
        pop_average_fitness.append(sum(pop_fitness_cycle)/len(pop_fitness_cycle))
        L2min_pairwise_diff_track.append(max(pairwise_cycle))
        L2min_column_track.append(max(min_column_cycle))
        
        
    return(cycle_track, best_probe_set, best_count, pop_best_fitness, pop_average_fitness, L2min_pairwise_diff_track, L2min_column_track)

# genetic algorithm for maximizing fitness score
# a single kmer is representative of a genes, a sets of kmer constitute a chromosome, a set of chromosomes represent a population
# initial population size is randomly chosen, initial chromosome size is the number of kmers
def gene_alg(cycles, chrom_per_pop, genes_per_chrom, kmer_enumerate, bacteria_dict, seed, alpha, offset):
    random.seed(seed)
    
    start_time = time.time()
    bacteria_num = len(bacteria_dict)
    gen_alg = genetic_alg(cycles, chrom_per_pop, genes_per_chrom, kmer_enumerate, bacteria_num, bacteria_dict, alpha, offset)
    gen_alg.pop_ini()
    
    cycle_track = []
    fitness_score = []
    pop_fitness_average = []
    exe_time = []
    L2_pairwise = []
    L2_column = []
    
    for generation in range(cycles):
        cycle_track.append(generation)
        
        gen_alg.pop_fitness_cal()
        gen_alg.select_mating_pool()
        gen_alg.crossover_mutation()
        fitness_score.append(np.max(gen_alg.pop_fitness))
        L2_pairwise.append(np.max(gen_alg.pop_L2pairwise))
        L2_column.append(np.max(gen_alg.pop_L2column))
        
        if generation == 0:
            best_probe_set = gen_alg.best_probe_set_cycle
            best_count = gen_alg.best_count_cycle
            
        if (np.max(gen_alg.pop_fitness) > np.max(fitness_score)):
            best_probe_set = gen_alg.best_probe_set_cycle
            best_count = gen_alg.best_count_cycle
            
        pop_fitness_average.append(gen_alg.pop_fitness_average)
#        print(gen_alg.pop_fitness)
        exe_time.append(time.time() - start_time)
        
    return(cycle_track, fitness_score, pop_fitness_average, best_probe_set, best_count, exe_time, L2_pairwise, L2_column)
   
    

    


# Construct the count matrix 
def kmer_count(kmer_num,bacteria_num, bacteria_dict, kmer_set):
    start_time = time.time()
    total_records_length = np.zeros((1,bacteria_num))
    bacteria_index = 0
    Count_Matrix = np.zeros((kmer_num,bacteria_num))
    
    for bacteria in bacteria_dict.values():
        for records in globals()[bacteria]:
            sequence = records.seq
            probe_index = 0
            total_records_length[0,bacteria_index] += len(sequence)
            for probe in kmer_set:
                
                true_probe_index = list(kmer_set).index(probe)
                Count_Matrix[true_probe_index,bacteria_index] += sequence.count(probe)
                probe_index += 1
                if (probe_index%1000) == 0:
                    print(time.time() - start_time)
        
        bacteria_index += 1
    
    return(Count_Matrix)

# Obtain the top ranked N probes in a set of kmer

    
# a class containing different ways the fitness score can be calculated
class UMD_fitness_score:
    def __init__(self, Count_Matrix, bacteria_dict):
        self.Count_Matrix = Count_Matrix
        self.bacteria_num = len(Count_Matrix[1,:])
        self.bacteria_dict = bacteria_dict
        self.L2norm_pairwise = {}
        for i in range(self.bacteria_num):
            for j in range(self.bacteria_num): 
                if i > j:
                   self. L2norm_pairwise[str(list(self.bacteria_dict.keys())[i]) + '-' + str(list(self.bacteria_dict.keys())[j])] = linalg.norm(self.Count_Matrix[:,i] - self.Count_Matrix[:,j])         
    
    def __str__(self):
        return('Fitness score class to calculate different fitness score for a kmer count matrixfor universal microbial diagnostic platform')
     
    def L2norm_pairwise(self):
        return(self.L2norm_pairwise)
        
    def L2norm_max(self):
        return(max(list(self.L2norm_pairwise.values())))
        
    def L2norm_ave(self):
        return(sum(list(self.L2norm_pairwise.values()))/len(list(self.L2norm_pairwise.values())))
        
    def L2norm_min_pairwise_diff(self):
        return(min(list(self.L2norm_pairwise.values())))
        
    def L2norm_min_column(self):
        self.column_norm = []
        for i in range(len(self.Count_Matrix[0])):
            self.column_norm.append(linalg.norm(self.Count_Matrix[:,i]))
            
        self.min_column = min(self.column_norm)
        return(self.min_column)
    
    
# class of genetic algorithm containing the steps for genetic algorithms
# a single kmer is representative of a genes, a sets of kmer constitute a chromosome, a set of chromosomes represent a population
# initial population size is randomly chosen, initial chromosome size is the number of kmers
        
class genetic_alg:
    def __init__(self, cycles, chrom_per_pop, genes_per_chrom, kmer_enumerate, bacteria_num, bacteria_dict, alpha, offset):
        self.kmer_enumerate = kmer_enumerate
        self.bacteria_num = bacteria_num
        self.bacteria_dict = bacteria_dict
        self.cycles = cycles
        self.chrom_per_pop = chrom_per_pop
        self.genes_per_chrom = genes_per_chrom
        self.pop_size = (chrom_per_pop, self.genes_per_chrom)
        self.pop = []
        self.pop_fitness = []
        self.pop_L2pairwise = []
        self.pop_L2column = []
        self.num_parent = 2
        self.parents = []
        self.offsprings = []
        self.pop_fitness_average = 0
        self.best_probe_set_cycle = 0
        self.best_count_cycle = 0 
        self.alpha = alpha
        self.offset = offset
        
    def __str__(self):
        return('Class of genetic algorithm containing steps customized for optimizing kmer for Universal Microbial Diagnostic platform')
    
    # initialize the initial population
    def pop_ini(self):
        for index in range(self.chrom_per_pop):
#            for kmerindex in range(self.genes_per_chrom):
#                chrom = random.sample(tuple(self.kmer_enumerate), self.genes_per_chrom)
#                self.pop[index][:] = chrom[kmerindex]
#            self.pop[index][:] = list(random.sample(tuple(self.kmer_enumerate), self.genes_per_chrom))
            self.pop.append(random.sample(tuple(self.kmer_enumerate), self.genes_per_chrom))
        return self.pop
    
    # calculate the fitness of each sets of kmer in the population
    def pop_fitness_cal(self):
        self.pop_fitness = []
        for i in range(self.chrom_per_pop):
            Count_Matrix = kmer_count(self.genes_per_chrom, self.bacteria_num, self.bacteria_dict, self.pop[i])
            fitness = UMD_fitness_score(Count_Matrix, bacteria_dict)
            self.pop_fitness.append(fitness.L2norm_min_pairwise_diff() + self.alpha*np.log(fitness.L2norm_min_column() - self.offset))
            self.pop_L2pairwise.append(fitness.L2norm_min_pairwise_diff())
            self.pop_L2column.append(fitness.L2norm_min_column())
            if (fitness.L2norm_min_pairwise_diff() + self.alpha*np.log(fitness.L2norm_min_column() - self.offset)) >= np.max(self.pop_fitness):
                self.best_probe_set_cycle = self.pop[i]
                self.best_count_cycle = Count_Matrix
            
        self.pop_fitness_average = sum(self.pop_fitness)/len(self.pop_fitness)
        return(self.pop_fitness)
            
    def select_mating_pool(self):
        fitness = copy(self.pop_fitness)
        self.parents = []
        for parent_num in range(self.num_parent):
            max_fitness_idx = np.argmax(fitness)
            self.parents.append(self.pop[max_fitness_idx])
            fitness[max_fitness_idx] = -999999
        
        return self.parents
    
    def crossover_mutation(self):
        # point at which kmer set is switched between 2 sets, here I take the middle
        crossover_point = int(self.genes_per_chrom/2)
        offspring_size = self.chrom_per_pop - len(self.parents)
        self.offsprings = []
        self.offsprings.append(self.parents[0][:crossover_point] + self.parents[1][crossover_point:])
        self.offsprings.append(self.parents[1][:crossover_point] + self.parents[0][crossover_point:])
#        print('parent size',len(self.parents))
#        print('offsprings size', len(self.offsprings))
        
        for k in range(offspring_size - 2):
            if int(k/2) == 0:
                mutated_offspring = copy(self.offsprings[0])
            else:
                mutated_offspring = copy(self.offsprings[1])
            
            mutated_offspring[random.randint(0,len(mutated_offspring)-1)] = random.sample(self.kmer_enumerate,1)[0]
            #print(self.offsprings)
            self.offsprings.append(mutated_offspring)
            
        self.pop = []
        
        for j in range(len(self.parents)):
            self.pop.append(self.parents[j])
        for k in range(len(self.offsprings)):
            self.pop.append(self.offsprings[k])
        
        return self.pop

            
        
        

    
        
        
        
    
    
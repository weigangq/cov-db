"""NK_model.py
Referenced paper:
    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9018209/

Original file is located at
    https://github.com/song88180/fitness-landscape-error
"""
import numpy as np
import numpy.random as nrand
import matplotlib.pyplot as plt
import pandas as pd
import itertools
import copy
import pickle
import multiprocessing
import argparse

# Functions to generate interaction matrices ------------------------------

def imatrix_rand(N,K):
    '''
    This function takes the number of N elements and K interdependencies
    and creates a random interaction matrix.
    '''
    Int_matrix_rand = np.zeros((N, N))
    for aa1 in np.arange(N):
        Indexes_1 = list(range(N))
        Indexes_1.remove(aa1)  # remove self
        np.random.shuffle(Indexes_1)
        Indexes_1.append(aa1)
        Chosen_ones = Indexes_1[-(K+1):]  # this takes the last K+1 indexes
        for aa2 in Chosen_ones:
            Int_matrix_rand[aa1, aa2] = 1  # we turn on the interactions with K other variables
    return(Int_matrix_rand)

def calc_fit(NK_land_, inter_m, Current_position, Power_key_):
    '''
    Takes the landscape and a given combination and returns a vector of fitness
    values for the vector of the N decision variables.
    '''
    Fit_vector = np.zeros(N)
    for ad1 in np.arange(N):
        Fit_vector[ad1] = NK_land_[np.sum(Current_position * inter_m[ad1]
                                          * Power_key_), ad1]
    return(Fit_vector)

def comb_and_values(NK_land_, Power_key_, inter_m):
    '''
    Calculates values for all combinations on the landscape. The resulting
    array contains:
    - the first columns indexed from 0 to N-1 are for each of the combinations
    - the column indexed N is for the total fit (average of the entire vector)
    '''
    Comb_and_value = np.zeros((2**N, N+1))  # to capture the results
    c1 = 0  # starting counter for location
    for c2 in itertools.product(range(2), repeat=N):
        # this takes time so be carefull with landscapes of bigger size
        Combination1 = np.array(c2)  # taking each combination
        fit_1 = calc_fit(NK_land_, inter_m, Combination1, Power_key_)
        Comb_and_value[c1, :N] = Combination1  # combination and values
        Comb_and_value[c1, N] = np.mean(fit_1)
        c1 = c1 + 1
    return(Comb_and_value)

def normalize(array):
    '''
    Normalize an array of value to the scale of 0 to 1
    '''
    MAX = np.max(array)
    MIN = np.min(array)
    return (array - MIN)/(MAX - MIN)

# NK landscape ---------------------------------------------------------------

# Change parameters to get fitness landscape of different variable site.
parser = argparse.ArgumentParser(description = 'Add number of variable sites')
parser.add_argument('N', type=int, help='Number of variable sites')
args = parser.parse_args()
N = args.N

# Create NK fitness landscape
NK_landscape_list = {i:[] for i in range(1, 51)}
Power_key = np.power(2, np.arange(N - 1, -1, -1))
K_loop = itertools.cycle(range(1,N)) # K is evenly sampled from 1 to N
for i in range(1,2):
    for j in range(1):
        #print(i,j)
        Landscape_data = []
        K = next(K_loop)
        Int_matrix = imatrix_rand(N,K).astype(int)
        
        # Create random fitness effect table for fitness calculation
        NK_land = np.random.rand(2**N, N)
        
        # Calculate fitness and combine that to FL
        NK_landscape_list[i].append(comb_and_values(NK_land, Power_key, Int_matrix))
        
        # Normalize fitness in FL
        NK_landscape_list[i][j][:,N] = normalize(NK_landscape_list[i][j][:,N])

# Take only unique haplotypes
NK_items = []
unique_haps = []
for k in range(NK_landscape_list[1][0].shape[0]):
    haps = [str(int(n)) for n in NK_landscape_list[1][0][k][0:10]]
    haps = "".join(haps)
    if haps in unique_haps:
        continue
    else:
        unique_haps.append(haps)
        NK_items.append((haps, NK_landscape_list[1][0][k][-1]))
NK_items.sort(key=lambda x: x[-1], reverse=True)

# Take only 200 haplotypes with unique fitness values
hap_num = 200
pos = len(NK_items)/hap_num
pos = int(pos)
NK_list = []
for i in range(0,hap_num):
    NK_list.append(NK_items[pos*i])

# Save output
with open('nk-landscape.tsv', 'w') as f:
    f.write("ranked_id\thap\tfit\tmodel\n")
    # Rank haps
    h_id = 0
    for n in NK_list:
        f.write(f"H{h_id:03d}\t{n[0]}\t{n[1]}\tNK\n")
        h_id += 1
"""NK model and RMF model
Referenced paper:
    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9018209/

Original file is located at
    https://github.com/song88180/fitness-landscape-error
"""
from asyncore import ExitNow
import numpy as np
import numpy.random as nrand
import matplotlib.pyplot as plt
import pandas as pd
import itertools
import copy
import pickle
import multiprocessing
import argparse
from scipy import special
from sklearn.metrics import jaccard_score

# Parameters -------------------------------------------------------------------
parser = argparse.ArgumentParser(description='Add number of variable sites')
parser.add_argument('-N', type=int, help='Number of variable sites (i.e., SNVs, single-nucleotide variants) (Default 10)', default = 10)
parser.add_argument('-hap_num', type=int, help='Number of haplotypes in output (default 100)', default = 100)
parser.add_argument('-model', '--model', choices = ['nk', 'rmf', 'poly', 'cov'], help='Fitness landscape model (default: nk)', default='nk')
parser.add_argument('-ref', '--ref', choices = ['beta', 'e484k', 'n501y', 'wuhan', 'delta'], help='Reference strand (default Wuhan)', default='wuhan')
parser.add_argument('-in', '--infile')

args = parser.parse_args()
N = args.N or 10
hap_num = args.hap_num or 2**N


# Landscapes --------------------------------------------------------------------
# NK landscape
def nk():
    # Functions to generate interaction matrices
    def imatrix_rand(N, K):
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
            Chosen_ones = Indexes_1[-(K + 1):]  # this takes the last K+1 indexes
            for aa2 in Chosen_ones:
                Int_matrix_rand[aa1, aa2] = 1  # we turn on the interactions with K other variables
        return (Int_matrix_rand)

    def calc_fit(NK_land_, inter_m, Current_position, Power_key_):
        '''
        Takes the landscape and a given combination and returns a vector of fitness
        values for the vector of the N decision variables.
        '''
        Fit_vector = np.zeros(N)
        for ad1 in np.arange(N):
            Fit_vector[ad1] = NK_land_[np.sum(Current_position * inter_m[ad1] * Power_key_), ad1]
        return (Fit_vector)

    def comb_and_values(NK_land_, Power_key_, inter_m):
        '''
        Calculates values for all combinations on the landscape. The resulting
        array contains:
        - the first columns indexed from 0 to N-1 are for each of the combinations
        - the column indexed N is for the total fit (average of the entire vector)
        '''
        Comb_and_value = np.zeros((2 ** N, N + 1))  # to capture the results
        c1 = 0  # starting counter for location
        for c2 in itertools.product(range(2), repeat=N):
            # this takes time so be carefull with landscapes of bigger size
            Combination1 = np.array(c2)  # taking each combination
            fit_1 = calc_fit(NK_land_, inter_m, Combination1, Power_key_)
            Comb_and_value[c1, :N] = Combination1  # combination and values
            Comb_and_value[c1, N] = np.mean(fit_1)
            c1 = c1 + 1
        return (Comb_and_value)

    def normalize(array):
        '''
        Normalize an array of value to the scale of 0 to 1
        '''
        MAX = np.max(array)
        MIN = np.min(array)
        return (array - MIN) / (MAX - MIN)

    # Create NK fitness landscape
    NK_landscape_list = {i: [] for i in range(1, 51)}
    Power_key = np.power(2, np.arange(N - 1, -1, -1))
    K_loop = itertools.cycle(range(1, N))  # K is evenly sampled from 1 to N
    for i in range(1, 2):
        for j in range(1):
            # print(i,j)
            Landscape_data = []
            K = next(K_loop)
            Int_matrix = imatrix_rand(N, K).astype(int)

            # Create random fitness effect table for fitness calculation
            NK_land = np.random.rand(2 ** N, N)

            # Calculate fitness and combine that to FL
            NK_landscape_list[i].append(comb_and_values(NK_land, Power_key, Int_matrix))

            # Normalize fitness in FL
            NK_landscape_list[i][j][:, N] = normalize(NK_landscape_list[i][j][:, N])

    # Combine haplotype and its corresponding fitness
    NK_items = []
    for k in range(NK_landscape_list[1][0].shape[0]):
        haps = [str(int(n)) for n in NK_landscape_list[1][0][k][0:-1]]
        haps = "".join(haps)
        NK_items.append((haps, NK_landscape_list[1][0][k][-1]))
    NK_items.sort(key=lambda x: x[-1], reverse=True)

    # Take only certain number of haplotypes with unique fitness values as output
    pos = len(NK_items) / hap_num
    pos = int(pos)
    NK_list = []
    for i in range(0, hap_num):
        NK_list.append(NK_items[pos * i])
    return NK_list, NK_landscape_list[1][0]

# RMF landscape
def rmf():
    # Fitness function
    # F(gt) = -cD(wt,gt)+N(std)

    # Fixed distance effect
    c = 1

    # Initialize genotype 0-1 space
    gt_lst = np.array(list(map(list, itertools.product([0, 1], repeat=N))))

    # Create Polynomial fitness landscape
    RMF_landscape_list = {i:[] for i in range(1,2)}
    for idx in range(1,2,2):
        for std in [10/n for n in range(1,2)]: # std is sampled from 0.5 to 10
            wt = nrand.randint(2,size=N) # Set wildtype genotype

            # Calculate fitness and combine that to FL
            fitness_lst = -c*np.sum(wt != gt_lst,axis=1)+nrand.normal(scale=std,size = np.power(2,N))

            # Normalize fitness
            MIN = np.min(fitness_lst)
            MAX = np.max(fitness_lst)
            fitness_lst = (fitness_lst - MIN) / (MAX - MIN)

            # Combine fitness to genotype 0-1 space
            fitness_landscape = np.concatenate((gt_lst, fitness_lst.reshape([-1, 1])), axis=1)

            # if there are 10 landscapes in the current RMF_landscape_list[idx],
            # go to the next idx
            if len(RMF_landscape_list[idx]) < 10:
                RMF_landscape_list[idx].append(fitness_landscape)
            else:
                RMF_landscape_list[idx + 1].append(fitness_landscape)

    # Combine haplotype and its corresponding fitness
    RMF_items = []
    for k in range(RMF_landscape_list[1][0].shape[0]):
        haps = [str(int(n)) for n in RMF_landscape_list[1][0][k][0:-1]]
        haps = "".join(haps)
        RMF_items.append((haps, RMF_landscape_list[1][0][k][-1]))
    RMF_items.sort(key=lambda x: x[-1], reverse=True)

    # Take only certain number of haplotypes with unique fitness values as output
    pos = len(RMF_items) / hap_num
    pos = int(pos)
    RMF_list = []
    for i in range(0, hap_num):
        RMF_list.append(RMF_items[pos * i])
    return RMF_list, RMF_landscape_list[1][0]

# Polynomial landscape
def poly():

    def get_fit(gt,b1,b2,b3):
        '''
        Takes a genoytpe and b1, b2, b3 parameters to calculate the fitness
        '''

        # Calculate additive fitnesse contribution
        fit_1 = np.sum(gt*b1)

        # Calculate 2-way epistatic fitness contribution
        gt_2 = np.array([])
        for i,row in enumerate(gt*gt.reshape(-1,1)):
            gt_2 = np.concatenate([gt_2,row[i+1:]])
        fit_2 = np.sum(gt_2*b2)

        # Calculate 3-way epistatic fitness contribution
        gt_3 = np.array([])
        gt_33 = gt.reshape(-1,1,1)*gt.reshape(1,-1,1)*gt.reshape(1,1,-1)
        for pos_1 in range(N-2):
            for pos_2 in range(pos_1+1,N-1):
                gt_3 = np.concatenate([gt_3,gt_33[pos_1,pos_2,pos_2+1:]])
        fit_3 = np.sum(gt_3*b3)

        # Return overall fitness
        return fit_1+fit_2+fit_3

    def get_fitness(gt_lst,b1,b2,b3):
        '''
        This function takes a list of genotypes and calculate their fitness
        '''
        # get fitness for all genotypes
        fitness_lst = np.array([get_fit(gt,b1,b2,b3) for gt in gt_lst])

        # normalize fitness
        MIN = np.min(fitness_lst)
        MAX = np.max(fitness_lst)
        fitness_lst = (fitness_lst - MIN) / (MAX - MIN)
        return fitness_lst

    def normalize(array):
        '''
        Normalize an array of value to the scale of 0 to 1
        '''
        MAX = np.max(array)
        MIN = np.min(array)
        return (array - MIN)/(MAX - MIN)

    # NK landscape parameters -----------------------------------------
    # Change parameters to get fitness landscape of different variable site.

    N = 15  # number of variable site

    # Initialize genotype 0-1 space
    gt_lst = np.array(list(map(list, itertools.product([0, 1], repeat=N))))

    # Create Polynomial fitness landscape
    POLY_landscape_list = {i:[] for i in range(1,101)}
    v1_loop = itertools.cycle(np.linspace(0.95,0.05,20)) # v1 is evenly sampled from 0.05 to 0.95
    for i in range(1,2):
        for j in range(1):
            Landscape_data = []
            v1 = next(v1_loop)
            v2 = np.sqrt((1-v1**2)/2); v3 = np.sqrt((1-v1**2)/2) # calculate v2 and v3

            # calculate b1 and b2 and b3
            b1 = nrand.normal(scale=v1,size=N)
            b2 = nrand.normal(scale=v2,size=special.comb(N,2).astype(int))
            b3 = nrand.normal(scale=v3,size=special.comb(N,3).astype(int))

            # calculate fitness for all genotypes
            fitness_lst = get_fitness(gt_lst,b1,b2,b3)
            fitness_landscape = np.concatenate((gt_lst,fitness_lst.reshape([-1,1])),axis=1)
            fitness_landscape = normalize(fitness_landscape)
            POLY_landscape_list[i].append(fitness_landscape)

    # Combine haplotype and its corresponding fitness
    POLY_items = []
    for k in range(POLY_landscape_list[1][0].shape[0]):
        haps = [str(int(n)) for n in POLY_landscape_list[1][0][k][0:-1]]
        haps = "".join(haps)
        POLY_items.append((haps, POLY_landscape_list[1][0][k][-1]))
    POLY_items.sort(key=lambda x: x[-1], reverse=True)

    # Take only certain number of haplotypes with unique fitness values as output
    pos = len(POLY_items) / hap_num
    pos = int(pos)
    POLY_list = []
    for i in range(0, hap_num):
        POLY_list.append(POLY_items[pos * i])
    return POLY_list, POLY_landscape_list[1][0]

    # with open(f'../FL_data_100X10/Polynomial_{N}_landscape_list_100X10.pkl','wb') as f:
    #     pickle.dump(Polynomial_landscape_list,f)

# COVID landscape
def cov(infile):

    # Get the reference strand
    if args.ref == 'beta':
        model_name = 'cov-beta'
        ref = 'Beta'

    if args.ref == 'e484k':
        model_name = 'cov-e484k'
        ref = 'E484K'

    if args.ref == 'n501y':
        model_name = 'cov-n501y'
        ref = 'N501Y'

    if args.ref == 'wuhan':
        model_name = 'cov-wuhan'
        ref = 'Wuhan-Hu-1'

    if args.ref == 'delta':
        model_name = 'cov-delta'
        ref = 'Delta'

    # Define fitness for each position as average of mutants fitness 
    def get_fit(input_file):
        '''
        Produces a list of all fitness values
        '''
        df = pd.read_csv(input_file)
        reference = df.groupby(['target']).get_group(ref)
        fit_list = reference.query('delta_bind > 0 or delta_bind <0').groupby('position').mean().reset_index().delta_bind.tolist()
        return fit_list
        
    def get_fitness(input_file):
        '''
        Converts the fitness list into an array
        '''
        fit_list = get_fit(input_file)
        fit_array = np.array(fit_list)
        return fit_array

    # COV landscape parameters -----------------------------------------

    N = 15  # number of variable site
    
    # calculate fitness for all genotypes
    fit_list = get_fit(infile)

    # Generate COV_list
    def COV_list():
        haps = []
        i = 1
        for item in fit_list:
            pasti = i
            j = 1
            k = fit_list.index(item) + 1
            a = f'{j:0{k}b}'
            single_hap = [a]
            norm = '0'
            while (i < len(fit_list)):
                single_hap.append(norm)
                i +=1
            single_hap = "".join(single_hap)
            haps.append(single_hap)
            i = pasti
            i += 1

        COV_list = []
        for item in haps:
            COV_item = (haps[haps.index(item)], fit_list[haps.index(item)])
            COV_list.append(COV_item)
        return COV_list

    # Generate COV_landscape_list
    def COV_landscape_list():
        haps = []
        i = 1
        for item in fit_list:
            pasti = i
            a = 1
            j = fit_list.index(item) + 1
            norm = 0
            one = 1
            single_hap = []
            while a < j:
                single_hap.append(norm)
                a +=1
            single_hap.append(one)
            while (i < len(fit_list)):
                single_hap.append(norm)
                i +=1
            haps.append(single_hap)
            i = pasti
            i += 1
        #print(haps)

        fitness = get_fitness(infile)

        COV_landscape = {i:[] for i in range(1,2)}
        COV_item_list = []
        #print(COV_landscape)
        for i in range(1,2):
            for item in haps:
                COV_item = haps[haps.index(item)]
                COV_item.append(fitness[haps.index(item)])
                COV_item_list.append(COV_item)
            COV_landscape[i].append(np.array(COV_item_list))
        return COV_landscape[1][0]

    COV_list = COV_list()
    COV_landscape_list = COV_landscape_list()

    return COV_list, COV_landscape_list, model_name

# Ruggedness -------------------------------------------------------
# Number of maxima (N max)
def get_N_max(landscape):
    N = landscape.shape[1] - 1
    N_max = 0
    if args.model == 'cov':
        hap = landscape[0]
        fit = hap[-1]
        flag = True
        for item in landscape:
            fit_ = item[-1]
            if fit < fit_:
                flag = False
            if flag == True:
                N_max += 1
            flag = True
            fit = fit_
    else:
        for gt in landscape:
            seq = gt[0:N]
            fit = gt[N]
            flag = True
            for i,_ in enumerate(seq):
                seq_ = copy.deepcopy(seq)
                seq_[i] = 1 - seq_[i]
                tmp = ''.join(seq_.astype(int).astype(str))
                idx = int(tmp, 2)
                fit_ = landscape[idx,N]
                if fit < fit_:
                    flag = False
                    break
            if flag == True:
                    N_max += 1
    return N_max

# Roughness to slope ratio (r/s)
from sklearn.linear_model import LinearRegression
def cal_r_s(landscape):
    N = landscape.shape[1] - 1
    X = landscape[:,:N]
    y = landscape[:,-1]
    reg = LinearRegression().fit(X, y) # fit_intercept default=True
    y_predict = reg.predict(landscape[:,:N])
    roughness = np.sqrt(np.mean(np.square(y - y_predict)))
    slope = np.mean(np.abs(reg.coef_))
    return roughness/slope

# Save output -------------------------------------------------------------------
# out_land is for landscape, out_rug if for ruggedness
if args.model == 'nk':
    model_name = 'nk'
    out_land = nk()[0]
    out_rug = nk()[1]

if args.model == 'rmf':
    model_name = 'rmf'
    out_land = rmf()[0]
    out_rug = rmf()[1]

if args.model == 'poly':
    model_name = 'poly'
    out_land = poly()[0]
    out_rug = poly()[1]

if args.model == 'cov':
    if args.infile == None:
        print("Need a COVID-binding file")
        exit
    infile = args.infile
    model_name = cov(infile)[2]
    out_land = cov(infile)[0]
    out_rug = cov(infile)[1]

N_max = get_N_max(out_rug)
r_s = cal_r_s(out_rug)
outFile_land = model_name + "-landscape.tsv"
with open(outFile_land, 'w') as f:
    f.write("ranked_id\thap\tfit\tmodel\tN_max\tr/s\n")
    # Rank haps
    h_id = 0
    for n in out_land:
        f.write(f"H{h_id:03d}\t{n[0]}\t{n[1]}\t{model_name}\t{N_max}\t{r_s}\n")
        h_id += 1

'''
outFile_rug = model_name + "-ruggedness.tsv"
with open(outFile_rug, 'w') as f:
    f.write(f"N_max\t{get_N_max(out_rug)}\n")
    f.write(f"r/s\t{cal_r_s(out_rug)}\n")
    f.write("hap\tfit\tmodel\n")
    for n in out_rug:
        f.write(f"{n[0:-2]}\t{n[-1]}\t{model_name}\n")
        h_id += 1
'''

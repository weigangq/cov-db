#!/usr/bin/env python
# Coding: utf-8
# obtaining evolutionary stats from CoV simulator

import logging
import sys
import argparse
import itertools
import pandas as pd
from progress.bar import Bar
import datetime
import re
from numpy.random import default_rng
import csv
import math

rng = default_rng() # Random number generator

# Initialize parser
parser = argparse.ArgumentParser(
    description="Calculate evolutionary statistics based on output files from cov-sim.py")

################################
# Arguments
##############################

parser.add_argument("-s", '--sites',
                    help='sites.tsv file')

parser.add_argument("-v", '--vars',
                    help='vars.tsv file')

parser.add_argument("-g", '--genes',
                    help='genes.tsv file')

parser.add_argument('-f', '--freq_cut', type=int, default = 500,
                    help='Minimal SNP frequency (default = 500, 0.5 percent of simulated samples)')

parser.add_argument('-p', '--pi', action = 'store_true',
                    help = 'print average pairwise nucleotide diffs per generation. Requires: -s sites.tsv')

parser.add_argument('-e', '--entropy', action = 'store_true',
                    help = 'print entropy for syn, mis, and total. Requires: -s sites.tsv')

'''
parser.add_argument('-q', '--pi_gene', action = 'store_true',
                    help = 'print average pairwise nucleotide diffs per gene. Requires: -s sites.tsv')
'''
parser.add_argument('-t', '--temporal', type = int,
                    help = 'print nucleotide diffs to a reference sample (default: a randome sequence in generation 1). Required: sites.tsv')

parser.add_argument('-a', '--trace_var', action = 'store_true',
                    help='Trace freq of most frequent SNPs, set by -f; Requires: -s sites.tsv -v vars.tsv files')

parser.add_argument('-b', '--trace_gene', action = 'store_true',
                    help='Trace var by genes; Requires: -s sites.tsv -v vars.tsv files')
#parser.add_argument('-g', '--gene', type = int, default = 100,
#                    help='Trace var freq with at least 100 (1%) presence. Required: vars.tsv')
'''
parser.add_argument('-c', '--coal',
                    help='Output: a coalescence tree. Required: lineages.tsv (to be implemented)')

parser.add_argument('-r', '--rate4site',
                    help='Map per site rates to SNPs. Required: vars.tsv (to be implemented)')
'''

######################################
# Program settings
#########################################

args = parser.parse_args()
logging.basicConfig(level = logging.DEBUG)
freqCut = args.freq_cut
################################
# Define all functions
##################################

def parse_gene_file(fname):
    genes = {}
    with open(args.genes, "r") as fh:
        lines = fh.readlines()
        for line in lines:
            if line.strip(): # line not empty
                data = line.split()
                genes[data[1]] = {'geneID': data[2],
                                'gene_length': data[2],
                                'gene_name': data[3]
                                    }
    return(genes)

def parse_var_file(fname):
    vars = {}
    with open(args.vars, "r") as fh:
        lines = fh.readlines()
        for line in lines:
            if line.strip(): # line not empty
                data = line.split()
                vars[data[2]] = {
                'tag': data[0],
                'pos': data[1],
                'varID': data[2],
                'locus': data[3],
                'refNT': data[4],
                'altNT': data[5],
                'refAA': data[6],
                'altAA': data[7],
                'refCodon': data[8],
                'altCodon': data[9],
                'conseq': data[10],
                'fit': data[11],
                'cts': data[12],
                'mult': data[13]
                }
    return(vars)

def parse_site_file(fname):
    ### load the site file, parse it
    df = pd.read_table(fname, header=None, sep="\t")
    df = df.fillna("")
    df.columns = ['runTag', 'gen', 'seqID', 'sites']
    df_stat = df.groupby(['gen']).size() #get numbers of generation
    samples = []
    tag = df.iloc[0,0]
    for gen in df_stat.index: # multiple samples in each gen
        subdata = df[df.gen==gen][['seqID', 'sites']].reset_index(drop=True)
        sites = [ subdata['sites'][i].split("|") for i in range(len(subdata)) ] # a list of sample sites
        syn = []
        mis = []
        all = []
        for sample_sites in sites: # for each sample
            sample_syn = []
            sample_mis = []
            sample_all = []
            for site in sample_sites:
                if re.search("_syn", site):
                    sample_syn.append(site)
                    sample_all.append(site)
                elif re.search("_mis", site):
                    sample_mis.append(site)
                    sample_all.append(site)
                elif re.search("_NA", site): # igs
                    sample_all.append(site)
                else: # nonsense, skip
                    continue
            #if len(sample_all) > 0:
            # append even if empty (diff = 0)
            all.append(sample_all)
            #if len(sample_mis) > 0:
            mis.append(sample_mis)
            #if len(sample_syn) > 0:
            syn.append(sample_syn)
        #print(syn)
        samples.append({
            'tag': tag,
            'generation':gen,
            'all_sites': all, # a list of 20 sample_all, one for each sample
            'syn_sites': syn,
            'mis_sites': mis
        })
    return(samples)

def avg_pair_diff(sites):
    '''
        Author: Lily (3/29/2021)
        A function that prints polymorphism change over generation to check mutation-drift equilibrium (at the expected level pi = 2Nu)

        Input:
            sites.tsv

        Returns:
            pairwise nt diff per generation
    '''
    diff = []
    #print(sites)
    #if len(sites) < 1: # no sample
    #    return 0
    #elif len(sites) < 2: # one sample
    #    return 1
#    print(sites)
    for k1, k2 in itertools.combinations(sites, 2):
        # empty sets works
        diff_pair = len(set(k1)|set(k2)) - len(set(k1)&set(k2))
        #Wrong: diff_pair = len(set(k1) - set(k2))
#        print(k1, "-", k2, "\t", diff_pair)
        diff.append(diff_pair)
    pi = sum(diff)/len(diff)
    return(pi)

def entropy(file):
    # model = file.split("/")[-1].split("-")[0]
    with open(file, "r") as test_sites:
        sim_data = list(csv.reader(test_sites, delimiter="\t"))

    # Find size of generation automatically
    gen_size = 0
    for indiv in sim_data:
        if int(indiv[1]) == 1:
            gen_size += 1
        else:
            break

    # Find number of generations automatically
    num_generations = int(sim_data[-1][1])
    model = sim_data[-1][0]
    # Parse sites into dictionary
    # Format: {gen 1: {site_num: {mut: #, mut, #}, site_num: {mut: #, mut:#}, ... } 
    #              gen: {site_num: {mut: #, ...}, site_num: {...} }
    #              ...
    #              gen 500: {...}
    #             }

    # e.g., {1: {13090: {'T_mis': 1}, 24637: {'A_mis': 1, 'C_syn': 1}}, 2: {13090: {'T_mis': 2}, 24637: {'A_mis': 2, 'C_syn': 1}}}

    sites = {}
    for gen in range(1, num_generations+1):
        sites[gen] = {}

    # indiv: ['neut', '1', 'Seq_1_1', '1202_C_mis|10989_G_syn']
    for indiv in sim_data:
        # Skip unmutated sequences
        if indiv[-1] == '':
            continue
        gen = int(indiv[1])
        mutations = indiv[-1].split("|") # ['1202_C_mis', '10989_G_syn']
        
        # For each indiv, if a new mutation occurs on a previously mutated site, overwrite the previous mutation
        filtered_mutations = {}
        for mutation in mutations: 
            # mutation: '1202_C_mis'
            # Skip NA mutations 
            if mutation[-2:] == "NA":
                continue
            m = mutation.split("_") # ['1202', 'C', 'mis']
            site = int(m[0]) # 1202
            mut = m[1] + "_" + m[2] # 'C_mis'
            filtered_mutations[site] = mut

        # Add the count to the sites dict
        for site in filtered_mutations:
            if site not in sites[gen]:
                sites[gen][site] = {}
                
            mut = filtered_mutations[site]
            if mut in sites[gen][site]:
                sites[gen][site][mut] += 1
            else:
                sites[gen][site][mut] = 1

    # Entropy calculation: for site i, H = summation( p(base) * log p(base)  for base = A, T, C, G)

    # gen_entropy: total entropy per generation
    # Format: {gen 1: 10 bits, 2: 5 bits, ...}
    gen_entropy = {}
    # Separate entropy into synonymous vs non-synonymous sites
    syn_entropy = {}
    mis_entropy = {}

    for gen in range(1, num_generations+1):
        total_entropy = 0 # Sum of entropy for all sites in a generation
        site_entropy = {} # Entropy for each mutated site    
        syn_total_entropy = 0 # Sum entropy for syn mutations
        syn_site_entropy = {} # For each syn mutations  
        mis_total_entropy = 0 # Sum entropy for non-syn mutations
        mis_site_entropy = {} # For non-syn mutations
        
        for site in sites[gen]: # Mutated sites 13090, 9403, ...
            site_entropy[site] = 0
            syn_site_entropy[site] = 0
            mis_site_entropy[site] = 0
            count_mutated = 0 # Total number of individuals with a mutation at that site
            syn_count = 0 # Count number of syn mutations at the site
            mis_count = 0 # Count number of mis mutations.

            for mut in sites[gen][site]: # C_syn, A_mis, ...
                freq = sites[gen][site][mut]
                count_mutated += freq
                if mut[-3:] == "mis":
                    mis_count += freq
                elif mut[-3:] == "syn":
                    syn_count += freq
            ancestral_count = gen_size - count_mutated
             # Calculate syn site entropy with population size only of ancestral + syn
            syn_pop = ancestral_count + syn_count
            mis_pop = ancestral_count + mis_count 
            
            for mut in sites[gen][site]: # C_syn, A_mis, ...
                freq = sites[gen][site][mut]
                count_mutated += freq
                
                # Calculate total entropy for the site
                site_entropy[site] += -1 * freq/gen_size * math.log(freq/gen_size, 2)
                
                # Calculate entropy with population size only of ancestral + mis or ancestral + syn
                if mut[-3:] == "mis":           
                    mis_site_entropy[site] += -1 * freq/mis_pop * math.log(freq/mis_pop, 2)
                elif mut[-3:] == "syn":
                    syn_site_entropy[site] += -1 * freq/syn_pop * math.log(freq/syn_pop, 2)

            if ancestral_count != 0:
                site_entropy[site] += -1 * ancestral_count/gen_size * math.log(ancestral_count/gen_size, 2)
                if mut[-3:] == "mis":           
                    mis_site_entropy[site] += -1 * ancestral_count/mis_pop * math.log(ancestral_count/mis_pop, 2)
                elif mut[-3:] == "syn":
                    syn_site_entropy[site] += -1 * ancestral_count/syn_pop * math.log(ancestral_count/syn_pop, 2)

        for s in site_entropy:
            total_entropy += site_entropy[s]
            syn_total_entropy += syn_site_entropy[s]
            mis_total_entropy += mis_site_entropy[s]
        
        gen_entropy[gen] = total_entropy
        syn_entropy[gen] = syn_total_entropy
        mis_entropy[gen] = mis_total_entropy
    
    for n in range(1, num_generations+1):
        print("{0}\t{1:.4f}\t{2}\t{3}".format(n, gen_entropy[n], "Total", model))
        print("{0}\t{1:.4f}\t{2}\t{3}".format(n, syn_entropy[n], "Synonymous", model))
        print("{0}\t{1:.4f}\t{2}\t{3}".format(n, mis_entropy[n], "Missense", model))
        

def diff_to_a_sample():
    '''
        A function that prints diffs to a chosen reference sample. This is for temporal analysis to identify viral origin

        Parameters:
            sites.tsv
            generation

        Returns:
            diffs to a reference
    '''

def coal():
    '''
        A function that reconstruct the coalescent tree of a simulated evolution. This helps to visualize effects of selection (background selection, adaptive, or neutral) to tree shape

        Parameters:
            lineages.tsv

        Returns:
            a tree in NEWICK format
    '''

def rate4site():
    '''
        This function maps rate4site to simulated SNPs, as controls for identifying adaptive mutations

        Parameters:
            rate4site
            vars.tsv

        Returns:
            mapped rates
    '''

logging.info("Start timestamp: %s", datetime.datetime.now())

sample_sites = []
if args.pi is True or args.entropy is True or args.trace_var is not None or args.temporal is not None or args.trace_gene is not None:
    if args.sites is None:
        logging.info("--pi, --trace_var, --temporal all need a -s site.tsv file as input")
        sys.exit()
    else:
        sample_sites = parse_site_file(args.sites)
        logging.info("sites file parsed")

if args.temporal is not None:
    genRef = args.temporal
    genSample = [ s for s in sample_sites if s['generation'] == genRef ][0]
    # print(genSample)
    countSite = len(genSample['all_sites'])
    if countSite < 1:
        logging.info("no mutations in generation %s. Pick another generation", genRef)
        sys.exit()
    pickRef = rng.choice(genSample['all_sites'])

    #print(pickRef)
    # pick one ind from each generation
    for sample in sample_sites:
        gen = sample['generation']
        #pickInd = rng.choice(sample['all_sites']) # pick one ind
        #print(pickRef, "\t", pickInd)
        for pickInd in sample['all_sites']: # all samples
            diff_pair = len(set(pickRef) | set(pickInd)) - len(set(pickRef) & set(pickInd))
            print(genRef, "\t", gen, "\t", diff_pair)

if args.trace_var is True:
    tag = ''
    if args.vars is None:
        logging.info("--trace_var needs a -v vars.tsv file as input")
        sys.exit()

    varDict = parse_var_file(args.vars)

    print("\t".join(["model", "snpID", "geneID", "half_life" ,"fit", "cts", "multi"]), "\t", "\t".join(map(str, range(1,501))))
    for varId in varDict:
        if int(varDict[varId]['cts']) < freqCut or varDict[varId]['conseq'] != 'missense':
            continue
        snpCts = list()
        tag = varDict[varId]['tag']
        begin = 0
        end = 0
        seen = 0
        for gen_sample in sample_sites:
            snpCt = 0
            gen = gen_sample['generation']
            misSites = gen_sample['mis_sites'] # all samples
            for site_sample in misSites: # one sample
                for site in site_sample: # one site
                    snpId = re.sub(r"_.{2,3}$", r"", site)
                    if snpId == varId:
                        snpCt += 1

            if snpCt > 0 and seen == 0: # identify first emergence
                begin = gen
                seen = 1

            if seen > 0 and (snpCt == 0 or gen == 500): # identify extinction or end
                end = gen

            snpCts.append(snpCt)
        snpCts = map(str, snpCts)
        print("\t".join([tag, varId, varDict[varId]['locus'],  str(end-begin+1), varDict[varId]['fit'], varDict[varId]['cts'], varDict[varId]['mult']]), "\t", "\t".join(snpCts))

if args.trace_gene is True:
    # print syn and mis by gene & generation
    tag = ''
    if args.vars is None:
        logging.info("--trace_gene needs a -v vars.tsv file as input")
        sys.exit()
    varDict = parse_var_file(args.vars)

    if args.genes is None:
        logging.info("--trace_gene needs a -g genes.tsv file as input")
        sys.exit()

    geneInfo = parse_gene_file(args.genes)

    print("\t".join(["model", "geneID", "gene_leng", "gen", "syn", "mis", "sum"]))
    geneDict = {}
    inLocus = {}
    for varId in varDict:
        tag = varDict[varId]['tag']
        gene = varDict[varId]['locus']
        inLocus[varDict[varId]['varID']] = gene
        if gene in geneDict:
            geneDict[gene].append(varDict[varId])
        else:
            geneDict[gene] = [ varDict[varId] ]

    for gene in geneDict: # for each gene
        if gene not in geneInfo:
            continue
        gene_len = geneInfo[gene]['gene_length']
        for gen_sample in sample_sites: # for each gen
            conseqCts = {'syn':0, 'mis':0, 'sum':0}
            gen = gen_sample['generation']

            for sample in gen_sample['mis_sites']:
                for site in sample:
                    snpId = re.sub(r"_.{2,3}$", r"", site)
                    if gene == inLocus[snpId]:
                        conseqCts['mis'] += 1

            for sample in gen_sample['syn_sites']:
                for site in sample:
                    snpId = re.sub(r"_.{2,3}$", r"", site)
                    if gene == inLocus[snpId]:
                        conseqCts['syn'] += 1

            for sample in gen_sample['all_sites']:
                for site in sample:
                    snpId = re.sub(r"_.{2,3}$", r"", site)
                    if gene == inLocus[snpId]:
                        conseqCts['sum'] += 1

            print("\t".join([tag, gene, gene_len, str(gen), str(conseqCts['syn']), str(conseqCts['mis']), str(conseqCts['sum'])]))

if args.pi is True:
    for sample in sample_sites: # for each generation
        tag = sample['tag']
        gen = sample['generation']
        all_sites = sample['all_sites']
        pi_all = avg_pair_diff(all_sites)
        syn_sites = sample['syn_sites']
        pi_syn = avg_pair_diff(syn_sites)
        mis_sites = sample['mis_sites']
        pi_mis = avg_pair_diff(mis_sites)
        #print(all_sites)
        print(tag, "\t", gen, "\t",
            "%.4f" % pi_all, "\t",
            "%.4f" % pi_syn, "\t",
            "%.4f" % pi_mis)

if args.temporal:
    diff_to_a_sample()

if args.entropy:
    entropy(args.sites)
'''
if args.coal:
    coal()

if args.rate4site:
    rate4site()
'''

logging.info("End timestamp: %s", datetime.datetime.now())
sys.exit()

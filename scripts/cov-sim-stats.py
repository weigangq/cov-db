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

rng = default_rng() # Random number generator

# Initialize parser
parser = argparse.ArgumentParser(
    description="Calculate evolutionary statistics based on output files from cov-sim.py")

################################
# Arguments
##############################

parser.add_argument('filename',
                    help='filename')

parser.add_argument('-p', '--pi', action = 'store_true',
                    help = 'print average pairwise nucleotide diffs per generation. Required file: sites.tsv')

parser.add_argument('-t', '--temporal', type = int,
                    help = 'print nucleotide diffs to a reference sample (default: a randome sequence in generation 1). Required: sites.tsv')

parser.add_argument('-c', '--coal',
                    help='Output: a coalescence tree. Required: lineages.tsv')

parser.add_argument('-r', '--rate4site',
                    help='Map per site rates to SNPs. Required: vars.tsv')

parser.add_argument('-a', '--trace',
                    help='Trace freq of a given SNP Required: sites.tsv')

######################################
# Program settings
#########################################

args = parser.parse_args()
logging.basicConfig(level = logging.DEBUG)

if args.filename is None:
    logging.info("Meed a sim out file as input")
    sys.exit()

fileName = args.filename

################################
# Define all functions
##################################

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
            if len(sample_all) > 0:
                all.append(sample_all)
            if len(sample_mis) > 0:
                mis.append(sample_mis)
            if len(sample_syn) > 0:
                syn.append(sample_syn)
        #print(syn)
        samples.append({
            'tag': tag,
            'generation':gen,
            'all_sites': all,
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
    if len(sites) < 1: # no sample
        return 0
    elif len(sites) < 2: # one sample
        return 1
#    print(sites)
    for k1, k2 in itertools.combinations(sites, 2):
        diff_pair = len(set(k1)|set(k2)) - len(set(k1)&set(k2))
        #Wrong: diff_pair = len(set(k1) - set(k2))
#        print(k1, "-", k2, "\t", diff_pair)
        diff.append(diff_pair)
    pi = sum(diff)/len(diff)
    return(pi)

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
if args.pi is True or args.trace is not None or args.temporal is not None:
    if re.search("sites", args.filename):
        sample_sites = parse_site_file(fileName)
        logging.info("sites file parsed")
    else:
        logging.info("--pi, --trace, --temporal all need a site.tsv file as input")

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

if args.trace is not None:
    # read site info in var.tsv format
    snps = []
    with open(args.trace, "r") as fh:
        lines = fh.readlines()
        for line in lines:
            if line.strip(): # line not empty
                data = line.split()
                snps.append(data[2])

    for snp in snps:
        for gen_sample in sample_sites:
            snpCt = 0
            gen = gen_sample['generation']
            misSites = gen_sample['mis_sites'] # all samples
            for site_sample in misSites: # one sample
                for site in site_sample: # one site
                    snpId = re.sub(r"_.{2,3}$", r"", site)
                    if snpId == snp:
                        snpCt += 1
            print(snp, "\t", gen, "\t", snpCt)

if args.pi is True:
    for sample in sample_sites:
        tag = sample['tag']
        gen = sample['generation']
        all_sites = sample['all_sites']
        pi_all = avg_pair_diff(all_sites)
        syn_sites = sample['syn_sites']
        pi_syn = avg_pair_diff(syn_sites)
        mis_sites = sample['mis_sites']
        pi_mis = avg_pair_diff(mis_sites)
        print(tag, "\t", gen, "\t",
            "%.4f" % pi_all, "\t",
            "%.4f" % pi_syn, "\t",
            "%.4f" % pi_mis)

if args.temporal:
    diff_to_a_sample()

if args.coal:
    coal()

if args.rate4site:
    rate4site()

logging.info("End timestamp: %s", datetime.datetime.now())
sys.exit()

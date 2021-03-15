#!/usr/bin/env python
# Coding: utf-8
# obtaining evolutionary stats from CoV simulator

import logging
import sys
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from progress.bar import Bar
import numpy as np
from numpy.random import default_rng
import datetime

rng = default_rng() # Random number generator

# Initialize parser
parser = argparse.ArgumentParser(
    description="Calculate evolutionary statistics based on output files from cov-sim.py")

################################
# Arguments
##############################

parser.add_argument('-p', '--pi_stat',
                    help = 'print average pairwise nucleotide diffs per generation. Required: sites.tsv')

parser.add_argument('-t', '--temporal',
                    help = 'print nucleotide diffs to a reference sample (default: a randome sequence in generation 1). Required: sites.tsv')

parser.add_argument('-g', '--gen', type = int, default = 1,
                    help = 'pick a sample as reference at a generation (default 1)')

parser.add_argument('-c', '--coal',
                    help='Output: a coalescence tree. Required: lineages.tsv')

parser.add_argument('-r', '--rate4site',
                    help='Map per site rates to SNPs. Required: vars.tsv')

######################################
# Program settings
#########################################

args = parser.parse_args()
logging.basicConfig(level = logging.DEBUG)

################################
# Define all functions
##################################
def avg_pair_diff():
    '''
        A function that prints polymorphism change over generation to check mutation-drift equilibrium (at the expected level pi = 2Nu)

        Input:
            sites.tsv

        Returns:
            pairwise nt diff per generation
    '''

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
if args.pi_stat:
    avg_pair_diff()

if args.temporal:
    diff_to_a_sample()

if args.coal:
    coal()

if args.rate4site:
    rate4site()

logging.info("End timestamp: %s", datetime.datetime.now())
sys.exit()

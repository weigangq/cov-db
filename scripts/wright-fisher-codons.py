#!/usr/bin/env python
# Coding: utf-8
# Simulate CoV genome evolution under Wright-Fisher model
# Main simulator
# To do:
# Done: 1. fix recurrence (for samples only; identify reversal)
# Done: 2. customize mutation counts
# Exponential pop growth?
# make mutation bias a constant

import logging
import sys
import argparse
import pprint
#from BCBio import GFF
#from BCBio.GFF import GFFExaminer
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from progress.bar import Bar
import copy
import numpy as np
from numpy.random import default_rng
import datetime

rng = default_rng() # Random number generator

# Initialize parser
parser = argparse.ArgumentParser(
    description="Wright-Fisher Simulator for SARS-CoV-2 genome evolution (version: 1.00). Required input: a Genbank file. Outputs: six files: tag-run.log (log of simulation parameters), tag-samples.tsv (mutations per sample, for plotting num mutations over generation), tag-genes.tsv (mutations per gene), tag-vars.tsv (snp info, VCF style), tag-sites.tsv (mutated sites per sample, for e.g., pi stats for use by sim-stats.py) and tag-lineages.tsv (lineage per sample, for reconstruction of coalescent tree). Three selection schemes are implemented: default is neutral (-u 0 -v 0). For background seleciton, use -v 0 (and -u e.g. 0.5). For adaptive selection, use -u 0 (and -v e.g., 0.5). Selection coefficient could be changed with -x (negative coefficeint, default 0.95) and -y (positive coefficient, default 1.05).")

################################
# Arguments
##############################

# Arguments for input & outputs
parser.add_argument('-b', '--genbank',
                    help='GenBank file (from NCBI). Required.')

# parser.add_argument('-f', '--gff',
#                    help='GFF file (from NCBI)')

parser.add_argument('-t', '--tag', default = 'test',
                    help='prefix for output files (default "out")')
parser.add_argument('-s', '--sample', type=int, default=10,
                    help='Sample size (default = 10)')

parser.add_argument('-f', '--fasta', action = 'store_true',
                    help='output sampled genome seqs (caution: large file size)')

# Arguments for basic evolution parameters
parser.add_argument('-g', '--generations', type=int, default=200,
                    help='Number of generations (default = 100)')

parser.add_argument('-p', '--population', type=int, default=100,
                    help='Population size (default = 100)')

#parser.add_argument('-n', '--gametes', type=int, default=20,
#                    help='Number of gametes produced by each genome (default = 20)')

parser.add_argument('-m', '--mutation', type=float, default=0.1,
                    help='Mutation rate per genome per generation (default = 0.1)')

parser.add_argument('-r', '--recombination', type=float, default = 0,
                    help='Recombination rate in the form of template switching(default=0.1)')

# Arguments for selection schemes

#parser.add_argument('-w', '--fitness', type = float, default = 1,                    help='fitness cost of a missense mutation (with respect to ref genome) under negative/purifying selection (default 0.95)')


parser.add_argument('-x', '--fit_neg', type = float, default = 0.95,
                    help='fitness cost of a missense mutation (with respect to ref genome) under negative/purifying selection (default 0.95)')

parser.add_argument('-y', '--fit_pos', type = float, default = 1.05,
                    help='fitness gain of a missense mutation (with respect to ref genome) under positive/adaptive selection (default 1.05)')

parser.add_argument('-u', '--prob_neg', type = float, default = 0,
                    help = 'fraction of missense mutations (with respect to ref genome) under negtive selection (default 0; set e.g., u=0.5 for background selection')

parser.add_argument('-v', '--prob_pos', type = float, default = 0,
                    help = 'fraction of missense mutations (with respect to ref genome) under positive selection (default 0; set e.g., v=0.05 for adaptive evolution)')

######################################
# Program settings
#########################################

args = parser.parse_args()
tagRun = args.tag
logging.basicConfig(filename = "%s-run.log" % tagRun,
                    filemode = "w",
                    level = logging.DEBUG)

if args.genbank is None:
    logging.info("Need a Genbank file as --genbank NC_045512.2.gb")
    sys.exit()

ref_gb = SeqIO.read(args.genbank, 'genbank')
genomeSeq = list(str(ref_gb.seq)) # make a list

# fit = args.fitness
fitNeg = args.fit_neg
fitPos = args.fit_pos
probNeg = args.prob_neg
probPos = args.prob_pos

if probNeg + probPos > 1:
    logging.info("-u and -v add up to more than one. Quit")
    sys.exit()

'''
def examine_gff():
    in_file = args.gff
    examiner = GFFExaminer()
    in_handle = open(in_file)
    pprint.pprint(examiner.parent_child_map(in_handle))
    pprint.pprint(examiner.available_limits(in_handle))
    in_handle.close()

def parse_gff():
    in_seq_file = args.fasta
    in_seq_handle = open(in_seq_file)
    seq_dict = SeqIO.to_dict(SeqIO.parse(in_seq_handle, "fasta"))
    in_seq_handle.close()

    in_handle = open(args.gff)
#    limit_info = dict(gff_id=["chr1"], gff_source=["Coding_transcript"])
    for rec in GFF.parse(in_handle, base_dict = seq_dict):
#        print(dir(rec))
        for feat in rec.features:
            print(feat)
    in_handle.close()

#examine_gff()
#parse_gff()
#sys.exit()
'''
################################
# Define all functions
##################################
def gene_locations(genbank):
    '''
        A function that reads a GenBank file and returns gene location range
        and translated protein.

        Parameters:
            genbank: GenBank file.

        Returns:
            gene_locations (list): A list of dictionaries containing
                                   information about translated protein
                                   and its coding gene location
                                   range in the genome.
    '''
    proteins = {}
    seen_product = {}
# exclude ORF1ab & nsp11
    exclude_product = { 'YP_009724389.1': 1, # ORF1ab
                        'YP_009725295.1': 1, # ORF1ab
                        'YP_009725312.1': 1 # nsp11, only 13aa, enclosed within nsp12
    }

    for feature in genbank.features:
        if feature.type == 'CDS' or feature.type == 'mat_peptide':
            # remove duplicated nsp's:
            if feature.qualifiers['product'][0] in seen_product:
                continue
            else:
                seen_product[feature.qualifiers['product'][0]] = 1

            # skip orf1ab & nsp11
            # but ORF7a and ORF7b overlap were kept
            for id in feature.qualifiers['protein_id']:
                if id in exclude_product:
                    continue
                else:
                    cds = feature.location.extract(genbank.seq) # properly handles compound location
                    pep = cds.translate()
                    proteins[id] = { 'id': id,
                                     'gene': feature.qualifiers['gene'][0],
                                     'location': feature.location, # compound for nsp12
                                     'product': feature.qualifiers['product'][0],
                                     'cds_seq': cds,
                                     'pep_seq': pep,
                                     'sample_synon': 0,
                                     'sample_missense': 0
                                   }
    return proteins
#    for pep in proteins:
#        print(pep)
#                    print(feature)
#            protein = feature.qualifiers['translation'][0] + '*' # Add 'stop'
#            gene_locations.append( {'location': feature.location, 'protein': protein} )



def get_codons(cds_location):
    '''
        A function that returns a list of coordinates for each codon.

        Parameters:
            gene_location (dict): A dictionary containing information about
                                  translated protein and its coding gene
                                  location range in the genome.

        Returns:
            codons (list): A list of tuples representing codons. Each tuple
                           value is a nucleotide position in a codon
                           (N1, N2, N3).
    '''
    codons = []
    positions = list(cds_location['location'])
    for index, position in enumerate(positions):
        if index % 3 == 0:
            second_codon_position = positions[index + 1]
            third_codon_position = positions[index + 2]
            codons.append(
                (position, second_codon_position, third_codon_position))
    return codons

def position_info(genbank):
    '''
        A function that takes in a reference GenBank file as argument and
        returns information about each coding nucleotide position.

        This information includes the nucleotide located at that position,
        the codon it is part of, codon coordinates, and translated amino acid.

        Parameters:
            genbank: GenBank file.

        Returns:
            position_data (dict): A dictionary containing information about
                                  a specific coding position in a genome.

                           position (int): {
                               'nucleotide': nucleotide at that position (str),
                               'codons': [
                                   'codon': codon found at that position (str),
                                   'amino_acid': translated amino acid (str),
                                   'location':  codon position
                                               (tuple of integer values)
                               ]
                           }
    '''
    # Get translated protein and coding range
#    genes = gene_locations(genbank)
    position_data = {}
#    position_data['seq'] = genbank.seq
    seen_position = {}

    for geneId in cdsObj: # for each gene
        # Get codon locations
        # print(gene['location'])
        gene = cdsObj[geneId]
        codons = get_codons(gene)
        # geneId = gene['id']
        geneProd = gene['product']
        geneName = gene['gene']
        # print(codons)
        # continue
        for index, codon_location in enumerate(codons): # for each codon; the location object handles slippage site correctly: site 13467 in two codons "AAC" & "CGG"
            # Translated amino acid
            amino_acid = gene['pep_seq'][index]
            codon = str(genbank[codon_location[0]:codon_location[2] + 1].seq)
            # Base position in a codon
            for position in codon_location: # for each codon position
#                validate slippage position:
#                if position > 13465 and position < 13470:
#                    print(position, end = "\t")
#                    print(codon)
                if position in seen_position: # assuming one position one frame (ribosome slippage site [13468, or 13467 index] actually serves in two reading frames; 2nd reading frame is ignored); also avoid overlapping genes (skipped in the get_locations function)
                    continue
                seen_position[position] = 1
                # Nucleotide at that position
                nucleotide = genbank[position]
                # Create a codon dictionary
#                codon_object = {
#                    'codon': codon, 'amino_acid': amino_acid, 'location': codon_location}
                position_data[position] = {
                    'nt': nucleotide,
                    'aa': amino_acid,
                    'location': codon_location,
                    'codon': codon,
                    'geneId': geneId,
                    'geneName':  geneName,
                    'geneProduct': geneProd
                }
#               position_data[position]['codons'].append(codon_object)
    return position_data

#                else:
                    # Only add unique codons for overlapping positions
#                    duplicates = 0
#                    for codon_obj in position_data[position]['codons']:
#                        if codon_obj == codon_object:
#                            duplicates += 1
#                    if not duplicates:


#reference_seq = SeqIO.read('ref.fas', 'fasta').seq.tomutable()
'''
for index in (13466, 13467, 13468):
    print(posInfo[index])
'''

def fitness(individual, mutation_site, new_base):
    '''
        Neutral + Negative selection model
        All relative to the ancestral genome, not with each other
        A function that checks if a mutation is synonymous or nonsynonymous
        and tracks nonsynonymous mutation sites.
        1. Nonsyn: fitness decreased (--fitness)
        2. Syn: fitness = 1.
        3. Stop to sense: fitness = 0
        4. Sense to stop: fitness = 0

        Parameters:
            individual (dict): A dictionary representing an individual in a
                               population.
            mutation_site (int): A position where mutation occurs in
                                 an individual's sequence.
            new_base (str): New nucleotide after mutation.
            positions (dict): A dictionary containing information
                                  about each individual coding position
                                  in a genome.

        Returns:
            individual (dict): Individual with a changed fitness value
                               if there was a nonsynonymous mutation.
    '''
    ind_seq = individual['seq']
    site = {
        'mut_site': mutation_site,
        'nt_ref': genomeSeq[mutation_site],
        'nt_pre': ind_seq[mutation_site],  # before mutation
        'nt_post': new_base,
        'gene_id': 'NA',
        'codon_ref': 'NA',
        'codon_pre': 'NA',
        'codon_post': 'NA',
        'aa_ref': 'NA',
        'aa_post': 'NA',
        'aa_pre': 'NA',
        'conseq': 'NA',
        'fit': 1.0
    }

    # Return individual as it is if a mutation occurs at a non-coding site
    if mutation_site not in posInfo:
        individual['igs'] += 1
    else: # in cds
        # Get information for selected position (codon and translated amino acid)
        position_info = posInfo[mutation_site]
        #    print(position_info)
        #    print(str(mutation_site))
        #    for codon in position_info['location']:
        gene = cdsObj[position_info['geneId']]
        index = position_info['location'].index(mutation_site) # What part of codon
        site['gene_id'] = position_info['geneId']
        site['nt_ref'] = position_info['nt']
        site['aa_ref'] = position_info['aa']
        site['codon_ref'] = position_info['codon']

        pre_codon = ''
        new_codon = ''
        if index == 0: # First base of codon
            base2 = mutation_site + 1
            base3 = mutation_site + 2
            pre_codon = ''.join(ind_seq[mutation_site : mutation_site + 3])
            new_codon = new_base + ind_seq[base2] + ind_seq[base3]
#           bug: new_codon = new_base + genomeSeq[base2] + genomeSeq[base3]
        elif index == 1: # Second base of codon
            base1 = mutation_site - 1
            base3 = mutation_site + 1
            pre_codon = ''.join(ind_seq[mutation_site-1 : mutation_site + 2])
            new_codon = ind_seq[base1] + new_base + ind_seq[base3]
#           bug: new_codon = genomeSeq[base1] + new_base + genomeSeq[base3]
        else: # Third base of codon
            base1 = mutation_site - 2
            base2 = mutation_site - 1
            pre_codon = ''.join(ind_seq[mutation_site-2 : mutation_site + 1])
            new_codon = ind_seq[base1] + ind_seq[base2] + new_base
#           bug: new_codon = genomeSeq[base1] + genomeSeq[base2] + new_base
        # BioPython Seq object
        pre_codon = Seq(pre_codon)
        new_codon = Seq(new_codon)
        site['codon_post'] = str(new_codon)
        site['codon_pre'] = str(pre_codon)
        pre_aa = str(pre_codon.translate())
        new_amino_acid = str(new_codon.translate())
        site['aa_pre'] = pre_aa
        site['aa_post'] = new_amino_acid

        # obtain fitness with respect to ref genome
        if new_amino_acid == position_info['aa']: # synonymous (including stop)
            individual['synon'] += 1
            site['fit'] = 1.0
            site['conseq'] = 'synonymous'
            individual['fitness'] *= 1.0
        elif new_amino_acid == '*' or position_info['aa'] == '*': # sense <=> stop
            individual['fitness'] = 0 # remove from gametes
            site['fit'] = 0
            site['conseq'] = 'nonsense'
        else: # nonsyn; ind fitness is cumulative site/codon fitness
            fit = rng.choice([fitNeg, 1.0, fitPos], p=[probNeg, 1 - probNeg - probPos, probPos]) # 89% missense neutral, 10% negative, 1% positive
            individual['missense'] += 1
            individual['fitness'] *= fit
            site['fit'] = fit
            site['conseq'] = 'missense'
#    print(site)
    individual['sites'].append(site)
#    np.append(individual['sites'], site)
    return individual

'''
        # Decrease fitness if mutation changes amino acid
        if new_amino_acid != codon['amino_acid']:
            individual['fitness'] *= fit
            individual['nonsyn'] = np.append(individual['nonsyn'], mutation_site)
        # Increase fitness if nonsynonymous mutation gets reversed
        elif new_amino_acid == codon['amino_acid'] and mutation_site in individual['nonsyn']:
            index = np.argwhere(individual['nonsyn'] == mutation_site)
            individual['nonsyn'] = np.delete(individual['nonsyn'], index)
            individual['fitness'] += 0.1
'''

def initialize(size):
    '''
        A function that creates an initial population of individuals.

        Parameters:
            size (int): Population size.
            reference_seq (Bio.Seq.MutableSeq): Nucleotide sequence.

        Returns:
            population (list): A list of dictionaries representing individuals.

                        {
                            'lineage': lineage IDs (list),
                            'sites': array of mutatated sites (numpy.ndarray),
                            'seq': nucleotide sequence (Bio.Seq.MutableSeq),
                            'fitness': individual fitness (float),
                            'nonsynonymous': nonsynonymous mutation sites (numpy.ndarray)
                        }
    '''
    population = []
    for i in range(size):
        population.append(
            {
                "lineage": [],
#                "sites": np.array([]),
                "sites": [],
                'seq': genomeSeq,
                'fitness': 1.0,
                'synon': 0,
                'missense': 0,
                'igs': 0
#                'nonsyn': np.array([])
            })
    return population

def mutate(individual, mutation_sites):
    '''
        A function that mutates a selected individual sequence and checks
        if the mutation was synonymous or nonsynonymous. Individual fitness
        gets decreased for each nonsynonymous mutation.

        Parameters:
            individual (dict): A dictionary representing an individual in
                               a population.
            mutation_sites (list): Sites where mutations occur.
            positions (dict): Dictionary with information about each
                                  coding position in a genome.

        Returns:
            individual (dict): Returns individual dictionary with a
                               mutated sequence and changed fitness
                               value for each nonsynonymous mutation.
    '''
    sequence = individual['seq']

#    gamete_sites = np.concatenate((individual['sites'], mutation_sites), axis=None)
#    individual['sites'] = gamete_sites

    # if multiple sites in a single codon, check one site at a time
    # mutation matrix to be made customizable
    for site in mutation_sites:
        site = int(site)
        mutated_base = sequence[site]
        if mutated_base == 'A':
            sequence[site] = rng.choice(['C', 'G', 'T'], p=[0.15, 0.65, 0.20])
        elif mutated_base == 'C':
            sequence[site] = rng.choice(['A', 'G', 'T'], p=[0.15, 0.05, 0.80])
        elif mutated_base == 'T':
            sequence[site] = rng.choice(['A', 'C', 'G'], p=[0.15, 0.70, 0.15])
        else:
            sequence[site] = rng.choice(['A', 'C', 'T'], p=[0.80, 0.05, 0.15])

        new_base = sequence[site]
        individual['seq'] = sequence # mutated seq, cumulative
        # Check if nonsynonymous mutation (with respect to ref, not to parent)
        individual = fitness(individual, site, new_base)
    return individual

def reproduction(pop, mut_rate):
    '''
        A function that performs reproduction with mutation
        in a population. Lineage IDs and mutated sites are appended
        to the individual dictionary objects.

        Parameters:
            pop (list): A list of dictionaries representing individuals
                        in a population.

            num_gametes (int): Initial number of produced gametes.

            genome_len: (int): Genome size.

            positions (dict): Dictionary containing information about each
                                  coding position in a genome.
        Returns:
            pool (list): A list of dictionaries representing gamete
                         individuals.
    '''

    '''
    for i in range(len(pop)):  # For each individual
        ind = pop[i]
        lineage = ind['lineage']

        # skip nonsense mutations or other way around
        if ind['fitness'] == 0:
            continue
        # Number of gametes based on fitness value
        size = int(num_gametes * ind['fitness'])
        for j in range(size):  # For each gamete
            gamete = copy.deepcopy(ind)
            gamete['lineage'] = lineage + [i] # list concat, not append
            num_mut = rng.poisson(mut_rate)  # Draw a Poisson number, mostly 0, 1
    '''
    pool = []
    gam_count = 0
    while gam_count < len(pop):
        ind = rng.choice(pop)
        lineage = ind['lineage']
        fit_ind = ind['fitness']
        # expm1 = exp(x) - 1; poisson k>0 probs as fitness cutoff:
        cut_off = -1 * np.expm1(-1*fit_ind)
        pick = rng.choice([1,0], p = [cut_off, 1 - cut_off])
        if pick == 1:
#            print("picked gamete no " + str(gam_count) + " for ind " + str(fit_ind))
            gamete = copy.deepcopy(ind)
            gamete['lineage'] = lineage + [gam_count] # list concat, not append
            num_mut = rng.poisson(mut_rate)  # Draw a Poisson number, mostly 0, 1
            gam_count += 1

            if num_mut > 0:
                mut_sites = rng.choice(len(genomeSeq), num_mut) # could be the same site
                # record recurrence
                for mut_site in mut_sites:
                    if mut_site in recurMutSites:
                        recurMutSites[mut_site] += 1
                    else:
                        recurMutSites[mut_site] = 1

                # Mutate individual
                gamete = mutate(gamete, mut_sites)

            pool.append(gamete)
    return pool

# Recombination
def recombination(pool, rec_rate):
    '''
        A function that performs recombination between 2 individual
        sequences.

        Parameters:
            pool (list): A list of dictionaries representing individuals.
            rec_rate (float): Recombination rate
            genome_len (int): Genome length

        Returns:
            pool (list): A modified list where with recombined individuals
    '''
    poisson = rng.poisson(rec_rate, len(pool)) # rate is per individual
    numRec = round(sum(poisson))
    for i in range(numRec):
        x1, x2 = rng.choice(len(pool), 2, replace = False)
        # copy parental gametes
        gam1 = {key: pool[x1][key] for key in ('sites', 'seq', 'lineage', 'fitness')}
        gam2 = {key: pool[x2][key] for key in ('sites', 'seq', 'lineage', 'fitness')}

        # pick breakpoint
        breakup = rng.choice(len(genomeSeq), 1)[0]

        # Make two recombinants
        # add seq:
        left1 = gam1['seq'][:breakup] # not including break point
        right1 = gam1['seq'][breakup:] # including break point
        left2 = gam2['seq'][:breakup]
        right2 = gam2['seq'][breakup:]
        gam1['seq'] = left1 + right2
        gam2['seq'] = left2 + right1

        # add fractional lineage:
        gam1['lineage'] += [ round(breakup/len(genomeSeq), 2) ]
        gam2['lineage'] += [ round(1 - breakup/len(genomeSeq), 2) ]

        # add sites; calculate fitness for each site & each recombinant
        gam1_sites = []
        gam2_sites = []

        for site in pool[x1]['sites']:
            mut_site = site['mut_site']
            if mut_site < breakup:
                gam1_sites.append(site)
                gam1 = fitness(gam1, mut_site, gam1['seq'][mut_site])
            else:
                gam2_sites.append(site)
                gam2 = fitness(gam2, mut_site, gam2['seq'][mut_site])

        for site in pool[x2]['sites']:
            mut_site = site['mut_site']
            if site < breakup:
                gam2_sites.append(site)
                gam2 = fitness(gam2, mut_site, gam2['seq'][mut_site])
            else:
                gam1_sites.append(site)
                gam1 = fitness(gam1, mut_site, gam1['seq'][mut_site])
        gam1['sites'] = gam1_sites
        gam2['sites'] = gam2_sites

        pool.append(gam1)
        pool.append(gam2)
        return pool

'''
        left_sites = gam2['sites'][gam2['sites'] < breakup]
        right_sites = gam1['sites'][gam1['sites'] >= breakup]
        # Nonsynonymous mutation sites
        nonsyn_left = gam2['nonsyn'][gam2['nonsyn'] < breakup]
        nonsyn_right = gam1['nonsyn'][gam1['nonsyn'] >= breakup]
        pool[x2]['sites'] = np.append(left_sites, right_sites)
        pool[x2]['nonsyn'] = np.append(nonsyn_left, nonsyn_right)
        pool[x2]['seq'] = recombined
        # Alter fitness
        pool[x2]['fitness'] = 1 - 0.1 * len(pool[x2]['nonsyn'])
'''

# produce gametes
def wright_fisher(pop, mut_rate, rec_rate):
    '''
        A function that produces gametes and a new population for next
        generation.

        Parameters:
            pop (list): A list of dictionaries representing individuals.
            genome_len (int): Genome length.
            num_gametes (int): Initial number of gametes (if fitness not changed).
            pop_size (int): Population size.
            mut_rate (float): Mutation rate.
            rec_rate (float): Recombination rate.
            positions (dictionary): A dictionary containing information
                                        about each position in a genome.

        Returns:
            gamete_next (list): A new population list of dictionaries
                                representing individuals.
    '''
    gamete_pool = reproduction(pop, mut_rate)
    if rec_rate > 0:
        gamete_pool = recombination(gamete_pool, rec_rate)
    return gamete_pool
#    gamete_next = rng.choice(gamete_pool, pop_size, replace=False)
#    return gamete_next

def outputVariant(gen, population, sample_size, fhInd, fhLine, fhSite, seqs, samp_sites):
    '''
        A function that takes a sample from the final population and outputs
        evolved sequences in a FASTA format.

        Parameters:
            population (list): A list of dictionaries representing individuals.
            sample_size (int): Number of samples.

        Outputs:
            evolved-sequences.fasta: A new FASTA file with sampled sequences.
    '''
    population_sample = rng.choice(population, sample_size, replace = False)
    for i in range(sample_size):
        samId = "_".join(['Seq', str(gen), str(i+1)])
        lineage = [str(l) for l in population_sample[i]['lineage']]
        lineage_info = "|".join(lineage)

        sequence_record = SeqRecord(
            Seq("".join(population_sample[i]['seq'])), id = samId, name = samId)
        seqs.append(sequence_record)

        fhInd.write(tagRun  + "\t" + str(gen) + "\t" + samId + "\t" + str(round(population_sample[i]['fitness'],4)) + "\t" + str(population_sample[i]['igs']) + "\t" + str(population_sample[i]['synon']) + "\t" + str(population_sample[i]['missense']) + "\n")

        if gen == generations:
            fhLine.write(tagRun  + "\t" + str(gen) + "\t"  + samId + "\t" + lineage_info  + "\n")

        sites = []
        for site in population_sample[i]['sites']:
            pos = site['mut_site']
            alt = site['nt_post']
            con = site['conseq']
            if con == 'NA':
                con = 'igs'
            snpID = str(pos) + "_" + alt
            sites.append(snpID + '_' + con[0:3]) # 'igs', 'syn', 'mis'

            if snpID in samp_sites:
                samp_sites[snpID]['count'] += 1
            else:
                samp_sites[snpID] = {'count': 1, 'info': site }

            if pos in posInfo:
                if pos in sample_gene_sites: # count unique gene sites
                    continue
                else:
                    sample_gene_sites[pos] = 1
                    position_info = posInfo[pos]
                    gene = cdsObj[position_info['geneId']]
                    gene['sample_synon'] += population_sample[i]['synon']
                    gene['sample_missense'] += population_sample[i]['missense']

        site_info = "|".join(sites)
        fhSite.write(tagRun  + "\t" + str(gen) + "\t" + samId + "\t" + site_info + "\n")


def simulation(num_gen, pop_size, mut_rate, rec_rate, sample_size):
    '''
        A main function that performs the Wright-Fisher simulation.

        Parameters:
            num_gen (int): Number of generations to run.
            pop_size (int): Population size.
            num_gametes (int): Initial number of gametes (if unchanged fitness)
            mut_rate (float): Mutation rate.
            rec_rate (float): Recombination rate.
            sample_size (int): Number of sequences to sample and write to file.
            reference_seq (Bio.Seq.MutableSeq): Nucleotide sequence used in
                                                the simulation.
            genbank: GenBank file.

        Outputs:
            Outputs a FASTA file containing evolved sequences at the end
            of the simulation.
    '''
    genome_len = len(genomeSeq)
    logging.info("Simulate CoV genome evolution under Wright-Fisher model")
    logging.info("Genome file: -b %s nt", args.genbank)
#    logging.info("Evolution model: -w %s (fitness discount for a missense mutation)", args.fitness)
    logging.info("Genome length: %s nt", genome_len)
    logging.info("Population size: -p %s individuals", pop_size)
    logging.info("Sample size: -s %s individuals per generation", args.sample)
#    logging.info("Gamete size: -n %s per individual", num_gametes)
    logging.info("Mutation rate: -m %s per genome per generation", mutation_rate)
    logging.info("Recombination rate: -r %s per genome per generation", rec_rate)
    logging.info("Total generations: -g %s", num_gen)
    logging.info("Negative selection against missense: -u %s probability with -x %s fitness cost", probNeg, fitNeg)
    logging.info("Positive selection for missense: -v %s probability with -y %s fitness gain", probPos, fitPos)
    logging.info("Neutral missense: %s probability with fitness of 1", 1 - probNeg - probPos)

    # output samples
    pop = initialize(pop_size)
    bar = Bar('Generation', max = num_gen) # Progress bar

    # open output file handles
    hapOut = open("%s-samples.tsv" % tagRun, "w")
    lineageOut = open("%s-lineages.tsv" % tagRun, "w")
    varOut = open("%s-vars.tsv" % tagRun, "w")
    siteOut = open("%s-sites.tsv" % tagRun, "w")

    seq_samples = []
    var_sites = {} # unique mutated sites in samples
    for generation in range(1, num_gen + 1):
        pop = wright_fisher(pop, mut_rate, rec_rate)
        outputVariant(generation, pop, sample_size, hapOut, lineageOut, siteOut, seq_samples, var_sites)
        bar.next()  # Progress

    # close output file handles
    hapOut.close()
    lineageOut.close()
    siteOut.close()

    for snp in var_sites:
        ct = var_sites[snp]['count']
        site = var_sites[snp]['info']
        pos = site['mut_site']
        id = site['gene_id']
        ref_nt = site['nt_ref']
        alt_nt = site['nt_post']
        ref_aa = site['aa_ref']
        alt_aa = site['aa_post']
        ref_codon = site['codon_ref']
        alt_codon = site['codon_post']
        fit = round(site['fit'],4)
        conseq = site['conseq']
        recur = recurMutSites[pos]

        varOut.write(tagRun + "\t" + str(pos) + "\t" + snp + "\t" + id + "\t" + ref_nt + "\t" + alt_nt + "\t" + ref_aa + "\t" + alt_aa + "\t" + ref_codon + "\t" + alt_codon + "\t" + conseq + "\t" + str(fit) + "\t" + str(ct) + "\t" + str(recur) +  "\n")
    varOut.close()

    if args.fasta:
        SeqIO.write(seq_samples, '%s-seqs.fasta' % tagRun, 'fasta')
        logging.info("Genome seq per sample written to file %s-seqs.tsv", tagRun)

    logging.info("Mutation counts per sample written to file %s-samples.tsv", tagRun)
    logging.info("Lineage info per sample written to file %s-lineages.tsv", tagRun)
    logging.info("Mutated sites per sample written to file %s-vars.tsv", tagRun)

    # output genes
    geneOut = open("%s-genes.tsv" % tagRun, "w")
    for geneId in cdsObj: # handles compound location beautifully
        geneOut.write(tagRun + "\t" + geneId + "\t" +  str(len(cdsObj[geneId]['location'])) + "\t" + cdsObj[geneId]['product'] + "\t" + str(cdsObj[geneId]['sample_synon']) + "\t" + str(cdsObj[geneId]['sample_missense']) + "\n")
    geneOut.close()
    logging.info("Mutation counts per gene written to file %s-genes.tsv", tagRun)
    bar.finish()
    # outputFasta(pop, sample_size)

logging.info("Start timestamp: %s", datetime.datetime.now())
generations = args.generations
pop_size = args.population
#num_gametes = args.gametes
sample_size = args.sample
mutation_rate = args.mutation
recombination_rate = args.recombination

recurMutSites = {}
cdsObj = gene_locations(ref_gb)
posInfo = position_info(ref_gb)

sample_gene_sites = {} # record unique sample sites in genes

simulation(
    generations,
    pop_size,
#    num_gametes,
    mutation_rate,
    recombination_rate,
    sample_size
)

# recurFile = open("%s-recur.tsv" % tagRun, "w")
# for site in recurMutSites:
#    if recurMutSites[site] > 1:
#        recurFile.write(str(site) + "\t" + str(recurMutSites[site]) + "\n")
#recurFile.close()
logging.info("End timestamp: %s", datetime.datetime.now())
sys.exit()

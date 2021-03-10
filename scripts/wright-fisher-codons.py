#!/usr/bin/env python
# Coding: utf-8
# Simulate CoV genome evolution under Wright-Fisher model
# Main simulator

import logging
import sys
import argparse
import pprint
from BCBio import GFF
from BCBio.GFF import GFFExaminer
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from progress.bar import Bar
import copy
import numpy as np
from numpy.random import default_rng
rng = default_rng() # Random number generator

# Initialize parser
parser = argparse.ArgumentParser(
    description='Wright-Fisher Simulator for Covid-19. Required input: a Genbank file. Outputs: two files: wf-sample.tsv (sample lineage and muations) & wf-gene.tsv (mutations in samples per gene)')

# Add arguments

parser.add_argument('-b', '--genbank', 
                    help='GenBank file (from NCBI)')

parser.add_argument('-f', '--gff', 
                    help='GFF file (from NCBI)')

parser.add_argument('-l', '--fasta', 
                    help='FASTA file (correponding to GFF)')

parser.add_argument('-w', '--fitness', type = float, default = 1,
                    help='fitness discount of nonsynonymous mutation (with respect to ref genome). Neutral by default. Set < 1 for negative selection')

parser.add_argument('-g', '--generations', type=int, default=100,
                    help='Number of generations (default = 100)')

parser.add_argument('-p', '--population', type=int, default=100,
                    help='Population size (default = 100)')

parser.add_argument('-n', '--gametes', type=int, default=10,
                    help='Number of gametes produced by each genome (default = 10)')

parser.add_argument('-s', '--sample', type=int, default=10,
                    help='Sample size (default = 10, limit = 100)')

parser.add_argument('-m', '--mutation', type=float, default=0.1,
                    help='Mutation rate per genome per generation (default = 0.1)')

parser.add_argument('-r', '--recombination', type=float, default = 0,
                    help='Recombination rate in the form of template switching(default=0.1)')

args = parser.parse_args()
logging.basicConfig(level = logging.DEBUG)

if args.genbank is None:
    logging.info("Need a Genbank file as --genbank NC_045512.2.gb")
    sys.exit()

fit = args.fitness
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
    position_data['seq'] = genbank.seq
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
    # Return individual as it is if a mutation occurs at a non-coding site
    if mutation_site not in posInfo:
        individual['igs'] += 1
        return individual

    # Get information for selected position (codon and translated amino acid)
    position_info = posInfo[mutation_site]
    #    print(position_info)
    #    print(str(mutation_site))
    #    for codon in position_info['location']:
    gene = cdsObj[position_info['geneId']]
    index = position_info['location'].index(mutation_site) # What part of codon
    new_codon = ''
    if index == 0: # First base of codon
        base2 = mutation_site + 1
        base3 = mutation_site + 2
        new_codon = new_base + genomeSeq[base2] + genomeSeq[base3]
    elif index == 1: # Second base of codon
        base1 = mutation_site - 1
        base3 = mutation_site + 1
        new_codon = genomeSeq[base1] + new_base + genomeSeq[base3]
    else: # Third base of codon
        base1 = mutation_site - 2
        base2 = mutation_site - 1
        new_codon = genomeSeq[base1] + genomeSeq[base2] + new_base
    # BioPython Seq object
    new_codon = Seq(new_codon) 
    new_amino_acid = new_codon.translate()
    if new_amino_acid == position_info['aa']: # synonymous (including stop)
        individual['synon'] += 1
        # gene['sample_synon'] += 1
        individual['fitness'] *= 1
    elif new_amino_acid == '*' or position_info['aa'] == '*': # sense <=> stop
        individual['fitness'] = 0 # remove from gametes
    else: # nonsyn
        # gene['sample_missense'] += 1
        individual['missense'] += 1
        individual['fitness'] *= fit                 
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
                "sites": np.array([]), 
                'seq': genomeSeq, 
                'fitness': 1,
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

    # if multiple sites in a single codon, check one site at a time
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

def reproduction(pop, num_gametes, mut_rate):
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
    pool = []
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
            gamete['lineage'] = lineage + [i]
            num_mut = rng.poisson(mut_rate)  # Draw a Poisson number, mostly 0, 1

            if num_mut > 0:
                mut_sites = rng.choice(len(genomeSeq), num_mut)
                # record recurrence
                for mut_site in mut_sites:
                    if mut_site in recurMutSites:
                        recurMutSites[mut_site] += 1
                    else:
                        recurMutSites[mut_site] = 0

                gamete_sites = np.concatenate((gamete['sites'], mut_sites), axis=None)
                gamete['sites'] = gamete_sites
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
            if site < breakup:
                gam1_sites.append(site)
                gam1 = fitness(gam1, site, gam1['seq'][site])
            else:
                gam2_sites.append(site)
                gam2 = fitness(gam2, site, gam2['seq'][site])

        for site in pool[x2]['sites']:
            if site < breakup:
                gam2_sites.append(site)
                gam2 = fitness(gam2, site, gam2['seq'][site])
            else:
                gam1_sites.append(site)
                gam1 = fitness(gam1, site, gam1['seq'][site])
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
def wright_fisher(pop, genome_len, num_gametes, pop_size, mut_rate, rec_rate):
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
    gamete_pool = reproduction(pop, num_gametes, mut_rate)
    if rec_rate > 0:
        gamete_pool = recombination(gamete_pool, rec_rate)
    gamete_next = rng.choice(gamete_pool, pop_size, replace=False)
    return gamete_next

# Take a sample from the final population and output sequences in FASTA format
def outputFasta(population, sample_size):
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
    sequences = []
    for i in range(sample_size):
        info = 'Seq_' + str(i+1)
        lineage = [str(l) for l in population_sample[i]['lineage']]
        lineage_info = "|".join(lineage)
        sites = [str(l) for l in population_sample[i]['sites']]
        sites_info = "|".join(sites)
        sequence_record = SeqRecord(
            Seq("".join(population_sample[i]['seq'])), id = info, name = info, description = lineage_info + "_" + sites_info)
        sequences.append(sequence_record)
    SeqIO.write(sequences, 'evolved-sequences.fasta', 'fasta')

def outputVCF(gen, population, sample_size, fhInd):
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
    sequences = []
    for i in range(sample_size):
        info = 'Seq_' + str(i+1)
        lineage = [str(l) for l in population_sample[i]['lineage']]
        lineage_info = "|".join(lineage)
        sites = [str(l) for l in population_sample[i]['sites']]
        sites_info = "|".join(sites)
        fhInd.write(str(gen) + "\t" + str(i) + "\t"  + info + "\t" + lineage_info + "\t" + sites_info + "\t" + str(population_sample[i]['fitness']) + "\t" + str(population_sample[i]['igs']) + "\t" + str(population_sample[i]['synon']) + "\t" + str(population_sample[i]['missense']) + "\n")

        for site in population_sample[i]['sites']:
            if site in posInfo:
                position_info = posInfo[site]
                gene = cdsObj[position_info['geneId']]
                gene['sample_synon'] += population_sample[i]['synon']
                gene['sample_missense'] += population_sample[i]['missense']

def simulation(num_gen, pop_size, num_gametes, mut_rate, rec_rate, sample_size):
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
    print("Simulate CoV genome evolution under Wright-Fisher model")
    print("* Genome length: %s nt" % genome_len, "\n",
          "* Population size:", pop_size, "\n",
          "* Gamete number:", num_gametes, "\n",
          "* Mutation rate:", mut_rate, "per genome per generation", "\n",
          "* Recombination rate:", rec_rate, "per genome per generation", "\n",
          "* Running generations:", num_gen, "\n",
          "Running......")
    # output samples
    pop = initialize(pop_size)
    bar = Bar('Generation', max = num_gen) # Progress bar
    hapOut = open("wf-sample.tsv", "w")
    for generation in range(1, num_gen + 1):
        pop = wright_fisher(pop, genome_len, num_gametes, pop_size, mut_rate, rec_rate)
        outputVCF(generation, pop, sample_size, hapOut)
        bar.next()  # Progress
    hapOut.close()

    # output genes
    geneOut = open("wf-gene.tsv", "w")
    for geneId in cdsObj:
        geneOut.write(geneId + "\t" + cdsObj[geneId]['product'] + "\t" + str(cdsObj[geneId]['sample_synon']) + "\t" + str(cdsObj[geneId]['sample_missense']) + "\n")
    geneOut.close()
    bar.finish()
    # outputFasta(pop, sample_size)

generations = args.generations
pop_size = args.population
num_gametes = args.gametes
sample_size = args.sample
mutation_rate = args.mutation
recombination_rate = args.recombination

ref_gb = SeqIO.read(args.genbank, 'genbank')
genomeSeq = list(str(ref_gb.seq)) # make a list
recurMutSites = {}
#mutInOrfs = {}
cdsObj = gene_locations(ref_gb)
#print(cdsObj)
posInfo = position_info(ref_gb)


simulation(
    generations, 
    pop_size, 
    num_gametes, 
    mutation_rate, 
    recombination_rate, 
    sample_size
)

for site in recurMutSites:
    if recurMutSites[site] > 1:
        print(recurMutSites[site])
sys.exit()



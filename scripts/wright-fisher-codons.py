#!/usr/bin/env python
# Coding: utf-8
# Simulate CoV genome evolution under Wright-Fisher model
# Main simulator

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
    description='Wright-Fisher Simulator for Covid-19')

# Add arguments

parser.add_argument('-b', '--genbank', 
                    help='GenBank file (from NCBI)')

parser.add_argument('-f', '--gff', 
                    help='GFF file (from NCBI)')

parser.add_argument('-l', '--fasta', 
                    help='FASTA file (correponding to GFF)')

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

parser.add_argument('-r', '--recombination', type=float, default=0.1,
                    help='Recombination rate in the form of template switching(default=0.1)')

args = parser.parse_args()

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
    proteins = []
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
                    proteins.append( { 'id': id,
                                       'gene': feature.qualifiers['gene'][0],
                                       'location': feature.location, # compound for nsp12
                                       'product': feature.qualifiers['product'][0],
                                       'cds_seq': cds,
                                       'pep_seq': pep
                                   })
    return proteins
#    for pep in proteins:
#        print(pep)
#                    print(feature)
#            protein = feature.qualifiers['translation'][0] + '*' # Add 'stop'
#            gene_locations.append( {'location': feature.location, 'protein': protein} )



def get_codons(gene_location):
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
    positions = list(gene_location['location'])
    for index, position in enumerate(positions):
        if index % 3 == 0:
            second_codon_position = positions[index + 1]
            third_codon_position = positions[index + 2]
            codons.append(
                (position, second_codon_position, third_codon_position))
    return codons


def position_info(cds):
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
    coordinates = gene_locations(genbank)
    position_data = {}
    position_data['seq'] = genbank.seq
    for gene in coordinates:
        # Get codon locations
        codons = get_codons(gene)

        for index, codon_location in enumerate(codons):
            # Translated amino acid
            amino_acid = gene['protein'][index]
            codon = str(genbank[codon_location[0]:codon_location[2] + 1].seq)
            # Base position in a codon
            for position in codon_location:
                # Nucleotide at that position
                nucleotide = genbank[position]
                # Create a codon dictionary 
                codon_object = {
                    'codon': codon, 'amino_acid': amino_acid, 'location': codon_location}
                if position not in position_data.keys():
                    position_data[position] = {
                        'nucleotide': nucleotide,
                        'codons': [codon_object]
                    }
                else:
                    # Only add unique codons for overlapping positions
                    duplicates = 0
                    for codon_obj in position_data[position]['codons']:
                        if codon_obj == codon_object:
                            duplicates += 1
                    if not duplicates:
                        position_data[position]['codons'].append(codon_object)
    return position_data


#reference_seq = SeqIO.read('ref.fas', 'fasta').seq.tomutable()
ref_gb = SeqIO.read(args.genbank, 'genbank')
cdsObj = gene_locations(ref_gb)
posInfo = position_info(cdsObj)
sys.exit()

def fitness(individual, mutation_site, new_base, positions):
    '''
        A function that checks if a mutation is synonymous or nonsynonymous
        and tracks nonsynonymous mutation sites. Individual fitness value 
        gets decreased for each nonsynonymous mutation. 
        Increases fitness value if a nonsynonymous mutation gets
        reversed. 

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
    if mutation_site not in positions:
        return individual

    # Get information for selected position (codon and translated amino acid)
    position_info = positions[mutation_site]
    sequence = positions['seq']
    for codon in position_info['codons']:
        index = codon['location'].index(mutation_site) # What part of codon
        if index == 0: # First base of codon
            base2 = mutation_site + 1
            base3 = mutation_site + 2
            new_codon = new_base + sequence[base2] + sequence[base3]
        elif index == 1: # Second base of codon
            base1 = mutation_site - 1
            base3 = mutation_site + 1
            new_codon = sequence[base1] + new_base + sequence[base3]
        else: # Third base of codon
            base1 = mutation_site - 2
            base2 = mutation_site - 1
            new_codon = sequence[base1] + sequence[base2] + new_base
        # BioPython Seq object
        new_codon = Seq(new_codon) 
        new_amino_acid = new_codon.translate()

        # Decrease fitness if mutation changes amino acid
        if new_amino_acid != codon['amino_acid']:
            individual['fitness'] -= 0.1
            individual['nonsyn'] = np.append(individual['nonsyn'], mutation_site)
        # Increase fitness if nonsynonymous mutation gets reversed
        elif new_amino_acid == codon['amino_acid'] and mutation_site in individual['nonsyn']:
            index = np.argwhere(individual['nonsyn'] == mutation_site)
            individual['nonsyn'] = np.delete(individual['nonsyn'], index)
            individual['fitness'] += 0.1
    return individual

def initialize(size, reference_seq):
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
                'seq': reference_seq,
                'fitness': 1,
                'nonsyn': np.array([]) 
            })
    return population

def mutate(individual, mutation_sites, positions):
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

    for site in mutation_sites:
        site = int(site)
        mutated_base = sequence[site]
        if mutated_base == 'A':
            sequence[site] = rng.choice(['C', 'G', 'T'], p=[0.13, 0.67, 0.20])
        elif mutated_base == 'C':
            sequence[site] = rng.choice(['A', 'G', 'T'], p=[0.14, 0.03, 0.83])
        elif mutated_base == 'T':
            sequence[site] = rng.choice(['A', 'C', 'G'], p=[0.16, 0.71, 0.13])
        else:
            sequence[site] = rng.choice(['A', 'C', 'T'], p=[0.41, 0.09, 0.50])
        
        new_base = sequence[site]
        individual['seq'] = sequence
        # Check if nonsynonymous mutation
        individual = fitness(individual, site, new_base, positions)
    return individual

def reproduction(pop, num_gametes, mut_rate, genome_len, positions):
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

        # Number of gametes based on fitness value
        size = int(num_gametes * ind['fitness'])
        for j in range(size):  # For each gamete
            gamete = copy.deepcopy(ind)
            gamete['lineage'] = lineage + [i]
            num_mut = rng.poisson(mut_rate)  # Draw a Poisson number, mostly 0, 1

            if num_mut > 0:
                mut_sites = rng.choice(genome_len, num_mut)
                gamete_sites = np.concatenate((gamete['sites'], mut_sites), axis=None)
                gamete['sites'] = gamete_sites
                # Mutate individual
                gamete = mutate(gamete, mut_sites, positions) 
                
            pool.append(gamete)
    return pool

# Recombination
def recombination(pool, rec_rate, genome_len):
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
    poisson = rng.poisson(rec_rate, len(pool))
    numRec = round(sum(poisson))
    for i in range(numRec):
        x1, x2 = rng.choice(len(pool), 2, replace=False)
        gam1 = {key: pool[x1][key] for key in ('sites', 'seq', 'nonsyn')}
        gam2 = {key: pool[x2][key] for key in ('sites', 'seq', 'nonsyn')}
        breakup = rng.choice(genome_len, 1)[0]
        # Recombine sequences
        leftSeq = gam2['seq'][:breakup]
        rightSeq = gam1['seq'][breakup:]
        recombined = leftSeq + rightSeq
        # Mutation sites
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
    return pool

# produce gametes
def wright_fisher(pop, genome_len, num_gametes, pop_size, mut_rate, rec_rate, positions):
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
    gamete_pool = reproduction(pop, num_gametes, mut_rate, genome_len, positions)
    if rec_rate > 0:
        gamete_pool = recombination(gamete_pool, rec_rate, genome_len)
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
    population_sample = rng.choice(population, sample_size, replace=False)
    sequences = []
    for i in range(sample_size):
        info = 'Seq_' + str(i+1)
        sequence_record = SeqRecord(
            population_sample[i]['seq'].toseq(), id=info, name=info, description=info)
        sequences.append(sequence_record)
    SeqIO.write(sequences, 'evolved-sequences.fasta', 'fasta')


def simulation(sequence, num_gen, pop_size, num_gametes, mut_rate, rec_rate, sample_size, genbank):
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
    genome_len = len(sequence)
    print("Simulate Cov genome evolution under Wright-Fisher model")
    print("* Genome length: %s nt" % genome_len, "\n",
          "* Population size:", pop_size, "\n",
          "* Gamete number:", num_gametes, "\n",
          "* Mutation rate:", mut_rate, "per genome per generation", "\n",
          "* Recombination rate:", rec_rate, "per genome per generation", "\n",
          "* Running generations:", num_gen, "\n",
          "Running......")
    positions = position_info(genbank)
    pop = initialize(pop_size, sequence)
    bar = Bar('Generation', max=num_gen) # Progress bar
    for generation in range(1, num_gen + 1):
        pop = wright_fisher(pop, genome_len, num_gametes, pop_size, mut_rate, rec_rate, positions)
        bar.next()  # Progress

    bar.finish()
    outputFasta(pop, sample_size)


generations = args.generations
pop_size = args.population
num_gametes = args.gametes
sample_size = args.sample
mutation_rate = args.mutation
recombination_rate = args.recombination

simulation(
    reference_seq,
    generations, 
    pop_size, 
    num_gametes, 
    mutation_rate, 
    recombination_rate, 
    sample_size, 
    ref_gb
    )

exit


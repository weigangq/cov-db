from binary_sim_population import *
import argparse
import logging
import sys
import re

######## Setup: arguments, parameters, and logging
parser = argparse.ArgumentParser(
    description="Search of fitness of random binary strings. Author: Winston Koh (Qiu Lab)")
parser.add_argument('-t', '--tag', default='test',
                    help='prefix for output files (default "test")')
                    
parser.add_argument('-land', '--landscape_file', required = True,
    help='landscape file generated by gen-fit-landscape.py. Required')

parser.add_argument('-pop', '--pop_size', type=int, default=100, help='Number of individuals in the population. Default = 100.')

parser.add_argument('-gen', '--generation', type=int, default=100, help='Number of generations to evolve. Default = 100.')

parser.add_argument('-mut', '--mutation_rate', type=int, default=1, help='Number of mutations per generation per haplotype, randomly selected by poisson distribution. Default = 1.')

parser.add_argument('-alg', '--algorithm', choices=['1', '2', '3'], default='1',
                    help='Evolutionary search algorithms. 1. Objective search: Selects the most fit individuals. 2. Novelty search: Selects the most novel individuals. 3. Combo search: Combines fitness and novelty.')

parser.add_argument('-n', '--nearest_neighbors', type=int, default=10,
                    help='Number of nearest neighbors to use in novelty search. Default = 10.')

parser.add_argument('-arch', '--archive_method', type=int, default=1,
                    help='Method to use for selecting what to add to the archive. 1. Fixed size, random archive. 2. Adaptive threshold. Default = 1.')

parser.add_argument('-prob', '--prob_arch', type=float, default=0.1, help='Probability for an individual to get randomly added to the archive during novelty search method 1. Default = 0.1.')

parser.add_argument('-w', '--weight', type=float, default=0.5, help='Weight for the combo algorithm. Weight value is between 0 and 1. Closer to 1 means more bias towards novelty. Default = 0.5.')

args = parser.parse_args()
tagRun = args.tag
logging.basicConfig(level=logging.DEBUG)

# Print simulation parameters
logging.info(f"Simulation parameters:\nLandscape file: {args.landscape_file}\n\tPopulation size: {args.pop_size}\n\tSearch algorithm: {args.algorithm}")

######## Read landscape file
if args.landscape_file is None:
    logging.info("Need a landscape file generated from the script gen-fit-landscape.py")
    sys.exit()

landscape = {}
with open(args.landscape_file, 'r') as fh:
    lines = fh.readlines()
    for line in lines:
        data = line.split()
        if re.match("^H\d+", data[0]): 
            landscape[data[0]] = {
                'hap': str(data[1]),
                'fit': float(data[2]),
                'model': data[3]
            }
            
#print(landscape)
#sys.exit()

############## initialize a population
p = Population(
    pop_size = args.pop_size, 
    landscape = landscape
    )

logging.info(f"starting haplotype on landscape: {p.start_hap}")
#print(p.pop)
#sys.exit()

# Fitness file (elite.tsv): each generation's 10 highest fitness strings and fitness.
#elite = open(tagRun + '-elite.tsv', 'w')
#elite.write('tag\tgen\tclosest_id\tclosest_hap\tclosest_fit\telite_id\telite_hap\tdiff_closest\tdiff_fittest\tlandscape\talgorithm\n')
#print(p.elite)
#sys.exit()
#print(p.pop[0])
############ search by evolution

print(f"Tag\tGen\tcl_id\tcl_fit\te1_id\tdiff_cl\tdiff_fit\tmodel\talgo")

for n in range(args.generation):
    # Mutate 1 site #print(f"{tagRun}\t{p.generation}\t{p.elite1['close_id']}\t{p.elite1['close_hap']}\t{p.elite1['close_fit']}\t{p.elite1['elite_id']}\t{arr_to_str(p.elite1['elite_hap'])}\t{p.elite1['diff_closest']}\t{p.elite1['diff_fittest']}\t{p.land_model}\t{args.algorithm}")
    print(f"{tagRun}\t{p.generation}\t{p.elite1['close_id']}\t{p.elite1['close_fit']}\t{p.elite1['elite_id']}\t{p.elite1['diff_closest']}\t{p.elite1['diff_fittest']}\t{p.land_model}\t{args.algorithm}")

    # end when reaches the peak
    if p.elite1['close_id'] == 'H000':
        break

    p.mutate(args.mutation_rate)

    if args.algorithm == '1':
        # Reproduce the 10 most fit/novel
        p.objective_selection()
    
    if args.algorithm == '2':
        p.novelty_selection(args.nearest_neighbors, args.archive_method, args.prob_arch)

    if args.algorithm == '3':
        p.combo(args.weight, args.nearest_neighbors, args.archive_method, args.prob_arch)
    
#elite.close()
#logging.info("Elite file written")
logging.info("Done")
sys.exit()

from binary_sim_population import *
import argparse
import logging
import sys

parser = argparse.ArgumentParser(
    description="Search of fitness of random binary strings. Author: Winston Koh (Qiu Lab)")
parser.add_argument('-t', '--tag', default='test',
                    help='prefix for output files (default "test")')
parser.add_argument('-r', '--run', default='test',
                    help='Run type: either create or read a landscape')                    
parser.add_argument('-s', '--seed', type=int, default=0, help='Seed for the population and high fitness string. Default = 0.')
parser.add_argument('-l', '--length', type=int, default=20, help='Length of each sequence. Default = 20.')
parser.add_argument('-pop', '--pop_size', type=int, default=100, help='Number of individuals in the population. Default = 100.')
parser.add_argument('-g', '--generations', type=int, default=100, help='Number of generations to evolve. Default = 100.')
parser.add_argument('-m', '--mutation_rate', type=int, default=1, help='Number of mutations per generation. Will eventually be changed to be random via poisson distribution. Default = 1.')
parser.add_argument('-a', '--algorithm', choices=['objective', 'novelty', 'combo', 'all'], default='all',
                    help='Evolutionary algorithm. Objective = Selects for most fit individuals. Novelty = Selects for most novel individuals. Combo = Combines fitness and novelty. Default = all.')
parser.add_argument('-land', '--landscape_model', choices=['1', '2', '3'], default='1',
                    help='Model of fitness landscape. Choices: 1 for additive fitness: -1 for every 1, and '
                         '+args.hfv if the string contains the high fitness string. 2 for fitness decided by a '
                         'normal distribution with mean=0 and sd=1. 3 for fitness decided by an exponential '
                         'distribution with rate=1. Default = 1')
#parser.add_argument('-hfl', '--high_fitness_length', type=int, default=10, help='Length of the high fitness string for the additive fitness landscape. Default = 10.')
#parser.add_argument('-hfv', '--high_fitness_value', type=int, default=20,
#                   help='Fitness value of the high fitness string for the additive fitness landscape. Default = 20.')
parser.add_argument('-n', '--nearest_neighbors', type=int, default=10,
                    help='Number of nearest neighbors to use in novelty search. Default = 10.')
parser.add_argument('-prob', '--archive_prob', type=float, default=0.1, help='Probability for an individual to get randomly added to the archive during novelty search. Default = 0.1.')
parser.add_argument('-w', '--weight', type=float, default=0.5, help='Weight for the combo algorithm. Weight value is between 0 and 1. Closer to 1 means more bias towards novelty. Default = 0.5.')

args = parser.parse_args()

tagRun = args.tag

logging.basicConfig(filename="%s-run.log" % tagRun,
                    filemode="w",
                    level=logging.DEBUG)

# Print simulation parameters
logging.info(f"Simulation parameters:\n Landscape model: {args.landscape_model}\n Population size: {args.pop_size}\n Genome length: {args.length}"
             f"\n Nearest neighbors: {args.nearest_neighbors}\n High "
#             f"fitness string length: {args.high_fitness_length}\n High fitness string value: {args.high_fitness_value}"
             )

#p = Population(length=args.length, pop=args.pop_size, model=int(args.landscape_model), hfs_length=args.high_fitness_length,
#               hfs_value= args.high_fitness_value, seed=args.seed)

p = Population(length=args.length, pop=args.pop_size, model=int(args.landscape_model), seed=args.seed)

# for printing landscape file
if p.model == 1:
    land_model = 'Additive'
elif p.model == 2:
    land_model = 'Normal'
else:
    land_model = 'Exponential'

# landscape.tsv file: All strings in the landscape and their fitness
with open(tagRun + '-' + land_model + '-landscape.tsv', 'w') as landscape_file:
    landscape_file.write('hap\tfit\tmodel\n')
    for n in p.landscape:
        landscape_file.write(f"{arr_to_str(n[1])}\t{n[0]}\t{land_model}\n")
logging.info("Landscape file written")

# Fitness file (elite.tsv): each generation's 10 highest fitness strings and fitness.
elite = open(tagRun + '-elite.tsv', 'w')
elite.write('gen\thap\tfit\tlandscape\talgorithm\n')

if args.algorithm == 'objective':
    for n in range(args.generations):
        # Mutate 1 site
        p.mutate(args.mutation_rate)

        # Reproduce the 10 most fit/novel
        p.objective_selection()
        elite.write(f"{p.generation}\t{arr_to_str(p.elite1[1])}\t{p.elite1[0]}\t{land_model}\t{args.algorithm}\n")
    
elif args.algorithm == 'novelty':
    for n in range(args.generations):
        p.mutate(args.mutation_rate)
        p.novelty_selection(args.nearest_neighbors, args.archive_prob)
        elite.write(f"{p.generation}\t{arr_to_str(p.elite1[1])}\t{p.elite1[0]}\t{land_model}\t{args.algorithm}\n")

elif args.algorithm == 'combo':
    for n in range(args.generations):
        p.mutate(args.mutation_rate)
        p.combo(args.weight, args.nearest_neighbors, args.archive_prob)
        elite.write(f"{p.generation}\t{arr_to_str(p.elite1[1])}\t{p.elite1[0]}\t{land_model}\t{args.algorithm}\n")

elif args.algorithm == 'all':
    algo = ''
    for n in range(args.generations):
        p.mutate(args.mutation_rate)
        p.objective_selection()
        algo = 'objective'
        elite.write(f"{p.generation}\t{arr_to_str(p.elite1[1])}\t{p.elite1[0]}\t{land_model}\t{algo}\n")
        
    for n in range(args.generations):
        algo = 'novelty'
        p.mutate(args.mutation_rate)
        p.novelty_selection(args.nearest_neighbors, args.archive_prob)
        elite.write(f"{p.generation}\t{arr_to_str(p.elite1[1])}\t{p.elite1[0]}\t{land_model}\t{algo}\n")

    for n in range(args.generations):
        algo = 'combo'
        p.mutate(args.mutation_rate)        
        p.combo(args.weight, args.nearest_neighbors, args.archive_prob)
        elite.write(f"{p.generation}\t{arr_to_str(p.elite1[1])}\t{p.elite1[0]}\t{land_model}\t{algo}\n")

else:
    print("search algorithm not implemented: ", args.algorithm)
    sys.exit()
 
    """
    # Write to file the 10 most fit/novel/scoring individuals
    for indiv in p.elite:
        elite.write(f"{p.generation}\t{arr_to_str(indiv[1])}\t{indiv[0]}\t{land_model}\t{args.algorithm}\n")
    """
    # Write to file the top 1 individual.

elite.close()
logging.info("Elite file written")
logging.info("Done")
sys.exit()
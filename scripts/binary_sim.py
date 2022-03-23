import numpy as np
import argparse
import logging
from numpy.random import default_rng
import sys

parser = argparse.ArgumentParser(description="Novelty + objective search of fitness of random binary strings. Author: Winston Koh (Qiu Lab)")
parser.add_argument('-l', '--length', type=int, default=50, help='Length of each sequence.')
parser.add_argument('-p', '--pop_size', type=int, default=100, help='Number of individuals in the population.')
parser.add_argument('-tr', '--trials', type=int, default=100, help='Number of trials to run.')
parser.add_argument('-n', '--nearest_neighbors', type=int, default=5,
                    help='Number of nearest neighbors to use when searching.')
parser.add_argument('-th', '--threshold', type=float, default=18.8, help='Threshold to use for the novelty archive.')
parser.add_argument('-hfl', '--high_fitness_length', type=int, default=10, help='Length of the high fitness string.')
parser.add_argument('-hfv', '--high_fitness_value', type=int, default=20,
                    help='Fitness value of the high fitness string.')
parser.add_argument('-t', '--tag', default = 'test',
                    help='prefix for output files (default "out")')
args = parser.parse_args()
rng = default_rng() # Random number generator
tagRun = args.tag

logging.basicConfig(filename = "%s-run.log" % tagRun,
                    filemode = "w",
                    level = logging.DEBUG)

# Print simulation parameters
logging.info("Simulation parameters:\n Population size: {%s}\n Genome length: {%s}\n Trials: {%s}", args.pop_size, args.length, args.trials)
logging.info(" Nearest neighbors: {%s}\n Threshold: {%s}", args.nearest_neighbors, args.threshold)
logging.info(" High fitness string length: {%s}\n High fitness string value: {%s}", args.high_fitness_length, args.high_fitness_value)

def arr_distances(nparray, compare_pop):
    # Find distance between nparray and every other individual in compare_pop.
    distances = []
    for indiv in range(compare_pop.shape[0]):
        # Hamming distance: number of element-wise differences.
        hamming = np.sum(np.absolute(nparray - compare_pop[indiv]))
        distances.append((hamming, compare_pop[indiv]))
    distances.sort(key=lambda x: x[0])  # Sort by distance.
    return distances

def highest_fitness(population, high_fitness_string):
    # Search for the individual with the highest fitness
    fitness = []
    for indiv in population:
        fitness.append([-1 * sum(indiv), tuple(indiv)])  # Fitness decreases by 1 for every 1 in the haplotype.
        # Fitness increases by 20 for every match in the high fitness sequences.
    for indiv in fitness:
        if high_fitness_string in str(indiv[1]):
            indiv[0] += args.high_fitness_value
    return max(fitness)

class Population:
    def __init__(self, length, pop=20, hfs_length=10, seed=None):
        self.archive = []  # Contains novel sequences in tuple format.
        self.pop_size = pop
        self.genome_length = length
        # self.objective_counter = 0
        self.seen = []
        self.high_fitness_string = str(tuple(np.random.default_rng(seed).integers(2, size=hfs_length)))[1:-1]
        self.genome = np.random.default_rng(seed).integers(2, size=(pop, length))
        self.highest_fitness = highest_fitness([tuple(self.genome[j]) for j in range(self.pop_size)],
                                               self.high_fitness_string)
        # self.genome = np.broadcast_to(np.random.default_rng(seed).integers(2, size=length), (pop, length)).copy()

    def novelty_search(self, threshold, nearest_neighbors=5):
        for num in range(self.pop_size):
            compare_pop = np.delete(self.genome, num, 0)
            distances = arr_distances(self.genome[num], compare_pop)
            sparsity = 0
            for t in distances[:nearest_neighbors]:
                sparsity += t[0]
            sparsity /= nearest_neighbors
            # print("Sparsity:", sparsity)
            if sparsity > threshold and tuple(self.genome[num]) not in self.archive:
                self.archive.append(tuple(self.genome[num]))
        if len(self.archive) == 0:
            raise ValueError('Archive is empty. Set threshold lower.')
        return highest_fitness(self.archive, self.high_fitness_string)

    def objective_fitness(self, nearest_neighbors=5, seed=None):
        # Randomly select 1 individual and recursively search for highest fitness among nearest neighbors.
        rand_num = np.random.default_rng(seed).integers(self.pop_size)
        rand_indiv = self.genome[rand_num]
        return self.objective_nn_search(rand_indiv, nearest_neighbors)

    def objective_nn_search(self, nparray, nearest_neighbors=5):
        # Take an individual as an input and search for the highest fitness among n nearest neighbors
        # self.objective_counter += 1
        compare_pop = self.genome
        distances = arr_distances(nparray, compare_pop)
        n_nearest = []  # Contains the array and the n nearest neighbors
        for d in distances[:nearest_neighbors + 1]:
            n_nearest.append(tuple(d[1]))  # Append the individual in tuple format.
            if hash(tuple(d[1])) not in self.seen:
                self.seen.append(hash(tuple(d[1])))
        highest = highest_fitness(n_nearest, self.high_fitness_string)
        if tuple(nparray) == tuple(highest[1]):
            return [highest[0], tuple(highest[1])]
        else:
            return self.objective_nn_search(highest[1], nearest_neighbors)


# Print results and save to file
#outfile = open(f"binary-sim-tr{args.trials}l{args.length}p{args.pop_size}n{args.nearest_neighbors}th{args.threshold}.tsv", "w")
#outfile.write("Trial\tTarget\tArchive_size\tNovelty_fitness\tSearch_size\tObjective_fitness\n")
print("Trial\tTarget\tArchive_size\tNovelty_fitness\tSearch_size\tObjective_fitness")

results = {"Objective": 0, "Novelty": 0}
for n in range(args.trials):
    p = Population(length=args.length, pop=args.pop_size, hfs_length=args.high_fitness_length, seed=n)

    objective = p.objective_fitness(nearest_neighbors=args.nearest_neighbors, seed=n+1)
    novelty = p.novelty_search(threshold=args.threshold, nearest_neighbors=args.nearest_neighbors)

    # See whether objective/novelty gets the highest fitness.
    if objective[0] == p.highest_fitness[0]:
        results["Objective"] += 1
    if novelty[0] == p.highest_fitness[0]:
        results["Novelty"] += 1
    
#    print(f"{n+1}\t{p.highest_fitness[0]}\t{len(p.archive)}\t{novelty[0]}\t{len(p.seen)}\t{objective[0]}")
#    outfile.write(f"{n+1}\t{p.highest_fitness[0]}\t{len(p.archive)}\t{novelty[0]}\t{len(p.seen)}\t{objective[0]}\n")
    print(f"{n+1}\t{p.highest_fitness[0]}\t{len(p.archive)}\t{novelty[0]}\t{len(p.seen)}\t{objective[0]}")

#outfile.close()

# Print final count of number of trials where objective/novelty gets the highest fitness.
logging.info(results)

sys.exit()

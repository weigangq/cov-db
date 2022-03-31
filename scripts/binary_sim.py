import numpy as np
import argparse
import logging

# from numpy.random import default_rng
# import sys

parser = argparse.ArgumentParser(
    description="Novelty + objective search of fitness of random binary strings. Author: Winston Koh (Qiu Lab)")
parser.add_argument('-l', '--length', type=int, default=50, help='Length of each sequence.')
parser.add_argument('-p', '--pop_size', type=int, default=100, help='Number of individuals in the population.')
parser.add_argument('-tr', '--trials', type=int, default=100, help='Number of trials to run.')
parser.add_argument('-n', '--nearest_neighbors', type=int, default=5,
                    help='Number of nearest neighbors to use when searching.')
parser.add_argument('-th', '--threshold', type=float, default=18.8, help='Threshold to use for the novelty archive.')
parser.add_argument('-hfl', '--high_fitness_length', type=int, default=10, help='Length of the high fitness string.')
parser.add_argument('-hfv', '--high_fitness_value', type=int, default=20,
                    help='Fitness value of the high fitness string.')
parser.add_argument('-t', '--tag', default='test',
                    help='prefix for output files (default "out")')
parser.add_argument('-s', '--seed', type=int, default=0, help='Seed for the population and high fitness string.')
parser.add_argument('-m', '--model', choices=['1', '2', '3'], default='add',
                    help='Model to evaluate fitness. Choices: \'add\' for additive fitness: -1 for every 1, and '
                         '+args.hfv if the string contains the high fitness string. \'norm\' for fitness decided by a '
                         'normal distribution with mean=0 and sd=1. \'exp\' for fitness decided by an exponential '
                         'distribution with rate=1')
"""
parser.add_argument('-add', '--additive_model', type=bool, default='True',
                    help='Evaluate fitness using the additive model.')
parser.add_argument('-norm', '--norm_dist_model', type=bool, default='False',
                    help='Evaluate fitness using a normal distribution with mean=0 and sd=1.')
parser.add_argument('-exp', '--exp_dist_model', type=bool, default='True',
                    help='Evaluate fitness using an exponential distribution with rate=1.')
"""
args = parser.parse_args()
tagRun = args.tag
# rng = default_rng()  # Random number generator

logging.basicConfig(filename="%s-run.log" % tagRun,
                    filemode="w",
                    level=logging.DEBUG)

# Print simulation parameters
logging.info(f"Simulation parameters:\n Fitness model: {args.model}\n Population size: {args.pop_size}\n Genome length: {args.length}\n Trials: "
             f"{args.trials}\n Nearest neighbors: {args.nearest_neighbors}\n Threshold: {args.threshold}\n High "
             f"fitness string length: {args.high_fitness_length}\n High fitness string value: {args.high_fitness_value}"
             )


def arr_distances(nparray, compare_pop):
    # Find distance between nparray and every other individual in compare_pop.
    distances = []
    for indiv in range(compare_pop.shape[0]):
        # Hamming distance: number of element-wise differences.
        hamming = np.sum(np.absolute(nparray - compare_pop[indiv]))
        distances.append((hamming, compare_pop[indiv]))
    distances.sort(key=lambda x: x[0])  # Sort by distance.
    return distances


# Model 1: Additive fitness
def fitness(population, high_fitness_string):
    # Search for the individual with the highest fitness
    fitness_list = []
    for indiv in population:
        fitness_list.append([-1 * sum(indiv), tuple(indiv)])  # Fitness decreases by 1 for every 1 in the haplotype.
        # Fitness increases by args.hfv for every match in the high fitness sequences.
    for indiv in fitness_list:
        if high_fitness_string in str(indiv[1]):
            indiv[0] += args.high_fitness_value
    return fitness_list


# Model 2. Epistasis with normal distribution: assign fitness by normal distribution, with mean=0 and sd=1,
# and pick the highest as the most fit
def norm_fitness(population, pop_size):
    fitness_list = []
    norm_dist = np.random.default_rng().normal(0, 1, pop_size)
    for indiv in range(len(population)):
        # Fitness is assigned randomly by normal distribution
        fitness_list.append([norm_dist[indiv], tuple(population[indiv])])
    return fitness_list


# Model 3. Epistasis, with exponential distribution: assign fitness by exponential distribution,
# with mean=1 and sd=1
def exp_fitness(population, pop_size):
    fitness_list = []
    exp_dist = np.random.default_rng().exponential(1, pop_size)
    for indiv in range(len(population)):
        # Fitness is assigned randomly by normal distribution
        fitness_list.append([exp_dist[indiv], tuple(population[indiv])])
    return fitness_list


class Population:
    def __init__(self, length, pop=20, hfs_length=10):
        self.archive = []  # Contains novel sequences as tuples.
        self.pop_size = pop
        self.genome_length = length
        self.objective_counter = 0  # Count of the path length in objective search
        self.seen = []  # Record the strings seen in objective search.
        self.high_fitness_string = str(tuple(np.random.default_rng(args.seed).integers(2, size=hfs_length)))[1:-1]
        self.genome = np.random.default_rng(args.seed).integers(2, size=(pop, length))
        self.highest_fitness = max(fitness([tuple(self.genome[j]) for j in range(self.pop_size)],
                                           self.high_fitness_string))
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
        return max(fitness(self.archive, self.high_fitness_string))

    def objective_fitness(self, trial, nearest_neighbors=5, seed=None):
        # Randomly select 1 individual and recursively search for highest fitness among nearest neighbors.
        rand_num = np.random.default_rng(seed).integers(self.pop_size)
        rand_indiv = self.genome[rand_num]
        self.seen = []  # Number of different strings encountered during objective search
        return self.objective_nn_search(rand_indiv, trial, nearest_neighbors)

    def objective_nn_search(self, nparray, trial, nearest_neighbors=5):
        # Take an individual as an input and search for the highest fitness among n nearest neighbors
        self.objective_counter += 1
        compare_pop = self.genome
        distances = arr_distances(nparray, compare_pop)
        n_nearest = []  # Contains the array and the n nearest neighbors
        for d in distances[:nearest_neighbors + 1]:
            n_nearest.append(tuple(d[1]))  # Append the individual in tuple format.
            if tuple(d[1]) not in self.seen:
                self.seen.append(tuple(d[1]))
        highest = max(fitness(n_nearest, self.high_fitness_string))

        # Check if the input array has the highest fitness
        if fitness([tuple(nparray)], self.high_fitness_string)[0][0] == highest[0]:
            if self.objective_counter == 1:
                # Write to global variable pathfile
                pathfile.write(f"{trial}\t{self.objective_counter}\t{highest[0]}\n")
            self.objective_counter = 0  # Reset step counter
            return [highest[0], tuple(highest[1])]
        else:
            pathfile.write(f"{trial}\t{self.objective_counter}\t{highest[0]}\n")
            return self.objective_nn_search(highest[1], trial, nearest_neighbors)


# Print results and save to file

# outfile = open(
#    f"binary-sim-tr{args.trials}l{args.length}p{args.pop_size}n{args.nearest_neighbors}th{args.threshold}.tsv", "w")
# outfile.write("Trial\tTarget\tArchive_size\tNovelty_fitness\tSearch_size\tObjective_fitness\n")
print("Trial\tTarget\tArchive_size\tNovelty_fitness\tSearch_size\tObjective_fitness")

p = Population(length=args.length, pop=args.pop_size, hfs_length=args.high_fitness_length)
novelty = p.novelty_search(threshold=args.threshold, nearest_neighbors=args.nearest_neighbors)

# Output files to record simulation results.
if args.model == '1':
    # pop.tsv: pop with additive fitness, a list of binary string, one string per line, with fitness
    pop_fitness = fitness(p.genome, p.high_fitness_string)
    with open("add_pop.tsv", "w") as popfile:
        popfile.write("String\tFitness\n")
        for n in range(len(pop_fitness)):
            # Write the string of digits without other characters.
            popfile.write(f"{''.join(map(str, pop_fitness[n][1]))}\t{pop_fitness[n][0]}\n")
elif args.model == '2':
    norm = norm_fitness(p.genome, p.pop_size)
    with open("norm_pop.tsv", "w") as normfile:
        normfile.write("String\tFitness\n")
        for n in norm:
            normfile.write(f"{''.join(map(str, n[1]))}\t{n[0]}\n")
elif args.model == '3':
    exp = exp_fitness(p.genome, p.pop_size)
    with open("exp_pop.tsv", "w") as expfile:
        expfile.write("String\tFitness\n")
        for e in exp:
            expfile.write(f"{''.join(map(str, e[1]))}\t{e[0]}\n")


# path.tsv: Path of the objective search showing trial, step, and fitness. Written in Population.objective_nn_search
pathfile = open("path.tsv", "w")
pathfile.write("Trial\tStep\tFitness\n")

results = {"Objective": 0, "Novelty": 0}
for n in range(args.trials):

    objective = p.objective_fitness(nearest_neighbors=args.nearest_neighbors, trial=n + 1, seed=n * 10)
    # novelty = p.novelty_search(threshold=args.threshold, nearest_neighbors=args.nearest_neighbors)

    # See whether objective/novelty gets the highest fitness.
    if objective[0] == p.highest_fitness[0]:
        results["Objective"] += 1
    if novelty[0] == p.highest_fitness[0]:
        results["Novelty"] += 1

    # outfile.write(f"{n+1}\t{p.highest_fitness[0]}\t{len(p.archive)}\t{novelty[0]}\t{len(p.seen)}\t{objective[0]}\n")
    print(f"{n + 1}\t{p.highest_fitness[0]}\t{len(p.archive)}\t{novelty[0]}\t{len(p.seen)}\t{objective[0]}")

# Log the final count of number of trials where objective/novelty gets the highest fitness.
logging.info(f"Results: {results}")

# outfile.close()
pathfile.close()

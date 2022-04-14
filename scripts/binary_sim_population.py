import numpy as np


def arr_to_str(nparray):
    # Convert 1D numpy array to a single string.
    return ''.join(str(nparray)[1:-1].split())


def arr_distances(nparray, compare_pop, genome_size):
    # Find distance between nparray and every other individual in compare_pop.
    distances = []
    for indiv in range(compare_pop.shape[0]):
        # Hamming distance: number of element-wise differences.
        hamming = np.sum(np.absolute(nparray - compare_pop[indiv]))
        distances.append((genome_size - hamming, compare_pop[indiv]))
    return distances  # List of [fitness, array]


# Model 1: Additive fitness
def fitness_add(population, high_fitness_string, high_fitness_value):
    # Search for the individual with the highest fitness
    fitness_list = []
    for indiv in range(population.shape[0]):
        fitness_list.append(
            [-1 * np.sum(population[indiv]), population[indiv]])  # Fitness decreases by 1 for every 1 in the haplotype.
        # Fitness increases by args.hfv for every match in the high fitness sequences.
    for indiv in fitness_list:
        if high_fitness_string in arr_to_str(indiv[1]):
            indiv[0] += high_fitness_value
    return fitness_list  # List of [fitness, array]


# Model 2. Epistasis with normal distribution: assign fitness by normal distribution, with mean=0 and sd=1
def fitness_norm(population, pop_size):
    fitness_list = []
    norm_dist = np.random.default_rng().normal(0, 1, pop_size)
    for indiv in range(population.shape[0]):
        # Fitness is assigned randomly by normal distribution
        fitness_list.append([norm_dist[indiv], population[indiv]])
    return fitness_list  # List of [fitness, array]


# Model 3. Epistasis, with exponential distribution: assign fitness by exponential distribution, with rate=1
def fitness_exp(population, pop_size):
    fitness_list = []
    exp_dist = np.random.default_rng().exponential(1, 2*pop_size)
    for indiv in range(population.shape[0]):
        # Fitness is assigned randomly by exponential distribution
        fitness_list.append([exp_dist[indiv], population[indiv]])
    return fitness_list  # List of [fitness, array]


class Population:
    def __init__(self, length=20, pop=20, model='1', hfs_length=10, hfs_value=20, seed=0):
        self.pop_size = pop
        self.genome_length = length
        self.generation = 0
        self.archive = None  # Contains sequences that were novel.
        self.elite = None  # Contains each generation's 10 most fit/novel/combo individuals + their fitness value.
        self.elite1 = None  # Top 1 individual.
        self.model = model  # Landscape fitness model.

        # Initial genome: matrix with dimensions pop_size x genome_length of a single random sequence. 
        self.genome = np.broadcast_to(np.random.default_rng(seed).integers(2, size=length), (pop, length)).copy()

        # Strings to be used for landscape
        self.landscape_strings = np.random.default_rng(seed + 1).integers(2, size=(2 * pop, length))

        # Determine the landscape and highest fitness individual based on the model
        if model == 1:
            # Randomly generated high fitness string for the additive fitness model.
            self.high_fitness_string = arr_to_str(np.random.default_rng(seed + 2).integers(2, size=hfs_length))
            self.landscape = fitness_add(self.landscape_strings, self.high_fitness_string, hfs_value)
        elif model == 2:
            self.landscape = fitness_norm(self.landscape_strings, 2 * pop)
        elif model == 3:
            self.landscape = fitness_exp(self.landscape_strings, 2 * pop)

        self.highest_fitness = max(self.landscape, key=lambda x: x[0])

    def mutate(self, count, gen=1):
        if count > self.genome_length or count < 1:
            raise ValueError('Number of digits to mutate must be between 1 and the string length.')
        for g in range(gen):
            self.generation += 1
            for num in range(self.pop_size):
                # TO-DO: Random number of mutations by poisson distribution

            
                # Random indices to mutate.
                mut_indices = np.random.default_rng().choice(self.genome_length, count, replace=False)
                # Mutate digits
                for i in mut_indices:
                    if self.genome[num, i] == 0:
                        self.genome[num, i] = 1
                    else:
                        self.genome[num, i] = 0
        return

    def replace_pop(self):
        # Get the strings of the elite
        elite10 = [n[1] for n in self.elite]

        # Create the new genome by broadcasting each individual to be 1/10 of the population size, then concatenate.
        self.genome = np.broadcast_to(elite10[0], (self.pop_size // 10, self.genome_length))
        for indiv in elite10[1:]:
            self.genome = np.concatenate(
                (self.genome, np.broadcast_to(indiv, (self.pop_size // 10, self.genome_length))), axis=0)
        return

    def objective_selection(self):
        """
        Get the 10 individuals with the highest fitness and their fitness value.
        Fitness is determined by the hamming distance from the highest fitness string in the landscape.
        More fit means smaller hamming distance to the highest fitness string
        """
        self.elite = sorted(arr_distances(self.highest_fitness[1], self.genome, self.genome_length),
                            key=lambda x: x[0])[-10:]
        self.elite1 = max(self.elite, key=lambda x: x[0])
        # Replace the population with the 10 fittest individuals
        self.replace_pop()
        return

    def novelty_search(self, nearest_neighbors=10, prob=0.10):
        novelty_list = []
        for num in range(self.pop_size):
            # Get the hamming distance for each individual to the rest of the genome and the archive.
            if self.archive is None:
                compare_pop = self.genome
            else:
                compare_pop = np.concatenate((self.archive, np.delete(self.genome, num, 0)), axis=0)
            # Get the distances to the genome+archive and sort to find the k nearest neighbors
            distances = sorted(arr_distances(self.genome[num], compare_pop, self.genome_length), key=lambda x: x[0])

            # Calculate the sparsity: average distance to k-nearest neighbors.
            sparsity = 0
            for t in distances[:nearest_neighbors]:
                sparsity += t[0]
            sparsity /= nearest_neighbors
            novelty_list.append([sparsity, self.genome[num]])

            """
            Using (Lehman and Stanley 2010) fixed size random archive.
            Each individual has a 10% chance of being added to the archive, regardless of their novelty value.
            The archive size stays fixed: when adding one, remove the oldest.
            """
            if np.random.default_rng().random() < prob:
                if self.archive is None:
                    self.archive = np.array([self.genome[num]])
                else:
                    # If the archive is full, remove the oldest individual.
                    if self.archive.shape[0] > 10 * self.pop_size:
                        self.genome = np.delete(self.genome, 0, axis=0)
                    self.archive = np.concatenate((self.archive, np.array([self.genome[num]])), axis=0)
        return novelty_list

    def novelty_selection(self, nearest_neighbors=10, prob=0.10):
        # Get the 10 individuals with the highest novelty
        novel = self.novelty_search(nearest_neighbors, prob)
        most_novel = sorted(novel, key=lambda x: x[0])[-10:]
        most_novel = np.array([n[1] for n in most_novel])  # Convert to numpy array

        # Get the fitness of the most novel
        self.elite = arr_distances(self.highest_fitness[1], most_novel, self.genome_length)
        self.elite1 = max(self.elite, key=lambda x: x[0])

        # Replace the population with the 10 most novel individuals
        self.replace_pop()
        return

    def combo(self, weight=0.5, nearest_neighbors=10, prob=0.10):
        """
        Give each individual a score based on both novelty and fitness.
        score(i) = (1 − ρ) · fit(i) + ρ · nov(i)
        ρ in [0, 1] controls the relative importance of fitness and novelty. High ρ biases towards novelty.
        """
        if weight > 1 or weight < 0:
            raise ValueError('The weight must be between 0 and 1.')
        pop_fitness = arr_distances(self.highest_fitness[1], self.genome, self.genome_length)
        pop_novelty = self.novelty_search(nearest_neighbors, prob)

        score_list = []
        for n in range(self.pop_size):

            fit = pop_fitness[n][0]
            nov = pop_novelty[n][0]
            score = (1 - weight) * fit + weight * nov
            score_list.append([score, self.genome[n]])

        # Get the 10 individuals with the highest scores
        high_score = sorted(score_list, key=lambda x: x[0])[-10:]
        high_score = np.array([n[1] for n in high_score])  # Convert to numpy array

        # Get the fitness of the highest scoring
        self.elite = arr_distances(self.highest_fitness[1], high_score, self.genome_length)
        self.elite1 = max(self.elite, key=lambda x: x[0])

        # Replace the population with the 10 highest scoring individuals
        self.replace_pop()
        return

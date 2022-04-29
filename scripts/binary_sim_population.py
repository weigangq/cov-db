import numpy as np
from numpy.random import default_rng


def arr_to_str(nparray):
    # Convert 1D numpy array to a single string.
    return ''.join(str(nparray)[1:-1].split())


def str_to_arr(string):
    # Convert a binary string into numpy array.
    return np.asarray(list(string), dtype=np.integer)


def dist_to_closest(nparray, land_pop):
    # Get the closest point in land_pop from nparray
    dist_to_land = []
    for id in land_pop:
        hamming = np.sum(np.absolute(nparray - np.asarray(list(land_pop[id]['hap']), dtype=np.integer)))
        dist_to_land.append({
            'fit': hamming,
            'land_id': id,
            'hap': nparray
        })
    closest_id = sorted(dist_to_land, key=lambda x: x['fit'])[0]
    # print(closest_id['hap'], "\t", closest_id['fit'])
    return closest_id


def dist_to_fittest(fittest, compare_pop, genome_size):
    # Find distance between nparray and every other individual in compare_pop.
    distances = []
    for indiv in range(compare_pop.shape[0]):
        # Hamming distance: number of element-wise differences.
        hamming = np.sum(np.absolute(fittest - compare_pop[indiv]))
        distances.append({
            'fit': hamming / genome_size,
            'hap': compare_pop[indiv]
        })
    return distances  # List of mutated haps


def hap_dist_to_fittest(fittest, nparray, genome_size):
    return np.sum(np.absolute(fittest - nparray)) / genome_size


# returns sorted list
def add_fitness(haps, model):
    # print(haps)
    # Search for the individual with the highest fitness
    fitness_list = []

    for indiv in range(haps.shape[0]):
        if model == '1':
            fitness_list.append({'fit': 1 - 1.0 * np.sum(haps[indiv]) / haps.shape[1],
                                 'hap': haps[indiv]})  # Fitness decreases by 1 for every 1 in the haplotype.

        if model == '2':
            norm_dist = np.random.default_rng().normal(0, 1, haps.shape[0])
            fitness_list.append({'fit': norm_dist[indiv], 'hap': haps[indiv]})

        if model == '3':
            exp_dist = np.random.default_rng().exponential(1, haps.shape[0])
            fitness_list.append({'fit': exp_dist[indiv], 'hap': haps[indiv]})

    if model == '4':  # monotonically decrease with the number of 1's
        norm_dist = np.random.default_rng().normal(0, 1, haps.shape[0])
        hap_sorted = sorted(haps, key=lambda x: np.sum(x), reverse=True)
        norm_dist = sorted(norm_dist)
        # print(norm_dist)
        # print(hap_sorted)
        for i in range(haps.shape[0]):
            fitness_list.append({'fit': norm_dist[i], 'hap': hap_sorted[i]})

    # there are ties if not add the zero string for additive model        
    if model == '1':
        fitness_list.append({'fit': 1.0, 'hap': np.zeros(haps.shape[1], dtype=int)})

    # sort by fitness
    fitness_list = sorted(fitness_list, key=lambda x: x['fit'], reverse=True)
    # print(fitness_list)
    return fitness_list  # List of [fitness, array]


class Population:
    def __init__(self, landscape, pop_size):
        self.generation = 0  # initialize when called
        self.elite = []  # Contains each generation's 10 most fit/novel/combo individuals + their fitness value.
        # self.elite1 = []  # Top 1 individual.
        self.pop_size = pop_size
        self.landscape = landscape

        # haps_land = [ landscape[id]['hap'] for id in landscape ] # list comprehension (foreach construct in Perl)
        top = landscape['H000']  # most fit ind in landscape
        # print(haps_land)
        start_ids = [n for n in list(landscape.keys()) if n != 'H000']  # don't allow to start on the peak
        pick = np.random.choice(start_ids, size=1)  # returns a list of ids
        self.genome_length = len(landscape[pick[0]]['hap'])
        self.highest_fitness = np.asarray(list(top['hap']), dtype=np.integer)
        self.land_model = top['model']
        pick_hap = np.asarray(list(landscape[pick[0]]['hap']), dtype=np.integer)
        pick_hap_str = landscape[pick[0]]['hap']
        pick_fit = landscape[pick[0]]['fit']
        self.start_hap = pick[0]

        self.archive = np.empty((0, self.genome_length), dtype=int)  # Contains sequences that were novel.

        # starting point on the landscape
        self.elite1 = {
            'gen': 0,
            'close_id': pick[0],  # closest haplotype on the landscape by hamming distance
            'close_hap': pick_hap_str,
            'close_fit': pick_fit,
            'elite_id': pick[0],  # elite
            'elite_hap': pick_hap,
            'diff_closest': 0,
            'diff_fittest': hap_dist_to_fittest(self.highest_fitness, pick_hap, self.genome_length)
        }

        # print(list(pick[0]))
        # print(len(pick[0]))
        # Initial genome: pick a random hap from the landscape (not a random string)
        # data structure of pop: 2D np array of np.integers (to get hamming distance)
        self.pop = np.broadcast_to(np.asarray(list(pick_hap), dtype=np.integer),
                                   (pop_size, self.genome_length)).copy()  # Population genome
        # print(self.pop.shape)

    def mutate(self, mut_rate):  # mutated a pop
        self.generation += 1
        # print("gen:\t", self.generation)
        for num in range(self.pop_size):  # mutate each hap
            # print("before mutation:\t", self.pop[num])
            num_mut = np.random.default_rng().poisson(mut_rate)  # Draw a Poisson number, mostly 0, 1
            if num_mut > 0:
                # Random indices to mutate.
                mut_indices = np.random.choice(self.genome_length, num_mut, replace=False)  # force change
                # Mutate digits
                for i in mut_indices:
                    if self.pop[num, i] == 0:
                        self.pop[num, i] = 1
                    else:
                        self.pop[num, i] = 0
            # print("after mutation:\t", self.pop[num])
        return

    def replace_pop(self):
        # Get the strings of the elite
        elite10 = [n['elite_hap'] for n in self.elite]
        # Create the new pop by broadcasting each individual to be 1/10 of the population size, then concatenate.
        # self.pop = np.broadcast_to(elite10, (self.pop_size // 10, self.genome_length))
        self.pop = np.empty((0, self.genome_length), dtype=int)  # !!! (do not use row dimension!!!)
        for indiv in elite10:
            self.pop = np.append(self.pop, np.broadcast_to(indiv, (self.pop_size // 10, self.genome_length)), axis=0)
        return

    def get_fitness(self, pop):
        # Get fitness of the all individuals in the pop based on the landscape
        close_pop = []
        id = 0
        for indiv in range(pop.shape[0]):
            closest = dist_to_closest(pop[indiv], self.landscape)  # Get hamming distances with landscape points

            close_pop.append({
                'gen': self.generation,
                'close_id': closest['land_id'],
                'close_fit': self.landscape[closest['land_id']]['fit'],
                'close_hap': self.landscape[closest['land_id']]['hap'],
                'elite_id': f"E{id:03d}",
                'elite_hap': self.pop[indiv],
                'diff_closest': closest['fit'],
                'diff_fittest': hap_dist_to_fittest(self.highest_fitness, self.pop[indiv], self.genome_length)
            })
            id += 1
        return close_pop

    def objective_selection(self):
        """
        Get the 10 individuals with the highest fitness and their fitness value.
        Fitness is determined by the fitness of its closest (most similar genotype) string on the landscape.
        """

        '''
        # select the top 10 closest to the fittest, pop evolves towards the fittest regardless of landscape
        elites = sorted(dist_to_fittest(self.highest_fitness, self.pop, self.genome_length), key=lambda x: x['fit'])[:10]
        for e in elites: # top 10
            closest = dist_to_closest(e['hap'], self.landscape)
            self.elite.append({
                    'gen': self.generation,
                    'close_id': closest['land_id'],
                    'close_fit': self.landscape[closest['land_id']]['fit'],
                    'close_hap': self.landscape[closest['land_id']]['hap'],
                    'elite_id': f"E{id:03d}",
                    'elite_hap': e['hap'],
                    'diff_closest': closest['fit'],
                    'diff_fittest': e['fit'],
                    })
            id += 1
        '''

        self.elite = []
        # Hill-climbing: select top 10 on the landscape (not necessarily close to the fittest in hamming distance
        # this makes the non-additive landscapes deceptive
        close_pop = self.get_fitness(self.pop)

        # Get fitness based on distance to the closest fitness landscape point
        self.elite = sorted(close_pop, key=lambda x: x['close_fit'], reverse=True)[:10]
        self.elite1 = self.elite[0]
        # print(len(self.elite1['close_hap']))
        # Replace the population with the 10 fittest individuals
        self.replace_pop()
        return

    def novelty_search(self, nearest_neighbors=10, method=1, prob=0.10):
        novelty_list = []

        # Parameters for novelty search method 2
        threshold = 10
        add_max = 15
        add_min = 5
        fraction = 0.5

        num_eval = 0
        counter = 0

        for num in range(self.pop_size):
            # Get the hamming distance for each individual to the rest of the genome and the archive.
            compare_pop = np.concatenate((self.archive, np.delete(self.pop, num, 0)), axis=0)

            # Get the distances to the genome+archive and sort to find the k nearest neighbors
            distances = []
            for indiv in range(compare_pop.shape[0]):
                # Hamming distance: number of element-wise differences.
                hamming = np.sum(np.absolute(self.pop[num] - compare_pop[indiv]))
                distances.append([hamming, compare_pop[indiv]])
            distances.sort(key=lambda x: x[0])

            # Calculate the sparsity: average distance to k-nearest neighbors.
            sparsity = 0
            for t in distances[:nearest_neighbors]:
                sparsity += t[0]
            sparsity /= nearest_neighbors
            novelty_list.append([sparsity, self.pop[num]])

            if method == 1:
                """
                Using fixed size random archive.
                (Lehman and Stanley 2010, Efficiently Evolving Programs through the Search for Novelty)
                Each individual has a 10% chance of being added to the archive, regardless of their novelty value.
                The archive size stays fixed: when adding one, remove the oldest.
                """
                if np.random.default_rng().random() < prob:
                    # If the archive is full, remove the oldest individual.
                    if self.archive.shape[0] > 2 * self.pop_size:
                        self.archive = np.delete(self.archive, 0, axis=0)
                    self.archive = np.concatenate((self.archive, np.array([self.pop[num]])), axis=0)

            elif method == 2:
                """
                Adding to archive using adaptive threshold.
                (Kistemaker and Whiteson 2011, Critical Factors in the Performance of Novelty Search)
                If the novelty of a new individual is higher than the threshold (ρ(x) > ρmin), it is added to the 
                archive. To keep the size of the archive approximately constant, ρmin is increased by a fixed fraction 
                if the number of added behaviors exceeds the addmax threshold in a certain number of evaluations. If the
                number of added behaviors is lower than addmin in a certain number of evaluations, ρmin is decreased by 
                a fixed fraction.
                """
                num_eval += 1
                if sparsity > threshold:
                    self.archive = np.concatenate((self.archive, np.array([self.pop[num]])), axis=0)
                    counter += 1
                if num_eval == 20:
                    if counter > add_max:
                        threshold += fraction
                    elif counter < add_min:
                        threshold -= fraction
                    # Reset counters
                    num_eval = 0
                    counter = 0

        return novelty_list

    def novelty_selection(self, nearest_neighbors=10, method=1, prob=0.10):
        self.elite = []  # Reset elite

        # Get the 10 individuals with the highest novelty
        novel = self.novelty_search(nearest_neighbors, method, prob)
        most_novel = sorted(novel, key=lambda x: x[0])[-10:]
        most_novel = np.array([n[1] for n in most_novel])  # Convert to numpy array

        # Get fitness of the 10 most novel individuals
        close_pop = self.get_fitness(most_novel)

        # Get fitness based on distance to the closest fitness landscape point
        self.elite = close_pop
        self.elite1 = sorted(close_pop, key=lambda x: x['close_fit'], reverse=True)[0]

        # Replace the population with the 10 most novel individuals
        self.replace_pop()
        return

    def combo(self, weight=0.5, nearest_neighbors=10, method=1, prob=0.10):
        """
        Give each individual a score based on both novelty and fitness.
        score(i) = (1 − ρ) · fit(i) + ρ · nov(i)
        ρ in [0, 1] controls the relative importance of fitness and novelty. High ρ biases towards novelty.
        """
        if weight > 1 or weight < 0:
            raise ValueError('The weight must be between 0 and 1.')
        pop_fitness = self.get_fitness(self.pop)
        pop_novelty = self.novelty_search(nearest_neighbors, method, prob)

        score_list = []
        for n in range(self.pop_size):
            fit = pop_fitness[n]['close_fit']
            nov = pop_novelty[n][0]
            score = (1 - weight) * fit + weight * nov
            score_list.append([score, self.pop[n]])

        # Get the 10 individuals with the highest scores
        high_score = sorted(score_list, key=lambda x: x[0], reverse=True)[:10]
        high_score = np.array([n[1] for n in high_score])  # Convert to numpy array

        # Get the fitness of the highest scoring
        self.elite = self.get_fitness(high_score)
        self.elite1 = max(self.elite, key=lambda x: x['close_fit'])

        # Replace the population with the 10 highest scoring individuals
        self.replace_pop()
        return

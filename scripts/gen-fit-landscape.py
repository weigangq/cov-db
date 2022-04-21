import argparse
import logging
import sys
import numpy as np
from binary_sim_population import *

### Interface design: arguments, parameters, and logging
parser = argparse.ArgumentParser(
    description="Generate a fitness landscape including haplotypes and their associated fitness values.Author: Winston Koh & Weigang Qiu")
parser.add_argument('-t', '--tag', default='test',
                    help='prefix for output files (default "test")')
parser.add_argument('-l', '--length', type=int, default=50, help='Length of binary haplotypes. Default = 20.')
parser.add_argument('-n', '--num_hap', type=int, default=500, help='Number of haplotypes to make the landscape. Default = 200.')
parser.add_argument('-m', '--landscape_model', choices=['1', '2', '3', '4'], default='1',
                    help=f"Model of fitness landscape.\nChoices: 1 Additive: -1 for every 1\n2 Normal: fitness decided by a normal distribution with mean=0 and sd=1.\n3 Exponental: fitness decided by an exponential distribution with rate=1.\n4 Monotonic: increases with number of 0s\nDefault = 1")
args = parser.parse_args()
tagRun = args.tag
size = args.num_hap
model = args.landscape_model
len = args.length
logging.basicConfig(filename="%s-landscape.log" % tagRun,
                    filemode="w",
                    level=logging.DEBUG)

# Print simulation parameters
logging.info(f"Simulation parameters:\nLandscape model: {args.landscape_model}\nNumber of haplotypes: {args.num_hap}\nHaplotype length: {args.length}")

#### Main body
# Strings to be used for landscape
landscape_strings = np.random.default_rng().integers(2, size=(size, len))
landscape = add_fitness(landscape_strings, model)
#print(landscape_strings.shape)
#sys.exit()

### Output
land_model = ''
if model == '1':
    land_model = 'additive'
if model == '2':
    land_model = 'normal'
if model == '3':
    land_model = 'exponential'
if model == '4':
    land_model = 'monotonic'
#print(landscape)
#highest_fitness = max(landscape, key=lambda x: x[0])

# landscape.tsv file: All strings in the landscape and their fitness
width = 10
precision = 4
with open(tagRun + '-' + land_model + '-landscape.tsv', 'w') as landscape_file:
    landscape_file.write('ranked_id\thap\tfit\tmodel\n')
    id = 0
    for n in landscape:
        landscape_file.write(f"H{id:03d}\t{arr_to_str(n['hap'])}\t{n['fit']:{width}.{precision}}\t{land_model}\n")
        id += 1
logging.info("Landscape file written")
logging.info("Done")
sys.exit()
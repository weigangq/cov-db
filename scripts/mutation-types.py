import csv
import math
import sys

# Error: generation 313 has site 10157 listed twice
# Generation size: 20
genome_size = 29903

with open("../sim-data/neut-sites.tsv", "r") as test_sites:
    sim_data = list(csv.reader(test_sites, delimiter="\t"))
    
    # Find size of generation automatically
    gen_size = 0
    for indiv in sim_data:
        if int(indiv[1]) == 1:
            gen_size += 1
        else:
            break
    
    # Find number of generations automatically
    num_generations = 0
    for indiv in sim_data:
        if int(indiv[1]) > num_generations:
            num_generations = int(indiv[1])
    
    sites = {}
    # indiv: ['neut', '1', 'Seq_1_1', '1202_C_mis|10989_G_syn']
    for indiv in sim_data:
        # Skip unmutated sequences
        if indiv[-1] == '':
            continue

        gen = int(indiv[1])
        if gen not in sites:
            sites[gen] = {}

        mutations = indiv[-1].split("|") # ['1202_C_mis', '10989_G_syn']
        
        filtered_mutations = {}
        for mutation in mutations: 
            # mutation: '1202_C_mis'
            m = mutation.split("_") # ['1202', 'C', 'mis']
            site = int(m[0]) # 1202
            mut = m[1] + "_" + m[2] # 'C_mis'
            
            filtered_mutations[site] = mut


        for site in filtered_mutations:
            if site not in sites[gen]:
                sites[gen][site] = {}
            
            mut = filtered_mutations[site]
            if mut in sites[gen][site]:
                sites[gen][site][mut] += 1
            else:
                sites[gen][site][mut] = 1


#print("Last generation:", num_generations, "\n", sites[num_generations])
    # sites format: {gen 1: {site_num: {mut: #, mut, #}, site_num: {mut: #, mut:#}, ... } 
    #                     gen: {site_num: {mut: #, ...}, site_num: {...} }
    #                     ...
    #                     gen 500: {...}
    #                    }
    
    # e.g., {1: {13090: {'T_mis': 1}, 24637: {'A_mis': 1, 'C_syn': 1}}, 2: {13090: {'T_mis': 2}, 24637: {'A_mis': 2, 'C_syn': 1}}}


#Entropy: for site i, H = summation( p(base) * log p(base)  for base = A, T, C, G)
#math.log(x, base)


#print(sites[313])
#for n in sites[313]:
#    print(n, sites[313][n])

# gen_entropy format: {gen 1: 10 bits, 2: 5 bits, ...}
gen_entropy = {}
for gen in range(1, num_generations+1):
    total_entropy = 0 # sum of entropy for all sites in a generation
    site_entropy = {} # entropy for each mutated site
    for site in sites[gen]: # mutated sites 13090, 9403, ...
        site_entropy[site] = 0
        count_mutated = 0
        for mut in sites[gen][site]: # C_syn, A_mis, ...
            freq = sites[gen][site][mut]
            count_mutated += freq
            site_entropy[site] += -1 * freq/gen_size * math.log(freq/gen_size, 2)
        ancestral_freq = gen_size - count_mutated
        if ancestral_freq != 0:
            site_entropy[site] += -1 * ancestral_freq/gen_size * math.log(ancestral_freq/gen_size, 2)

    for e in site_entropy:
        total_entropy += site_entropy[site]
    #count_nonmutated = genome_size - len(site_entropy)

    gen_entropy[gen] = total_entropy

#print("Total entropy per generation:")

# Save output to file
sys.stdout = open("entropy.tsv", "w")
for n in range(1, num_generations+1):
    print(n, "\t", gen_entropy[n])
sys.stdout.close()

"""
# Old mutations count per site

        for mutation in mutations:
            # mutation: '1202_C_mis'
            m = mutation.split("_") # ['1202', 'C', 'mis']
            site = int(m[0]) # 1202
            mut = m[1] + "_" + m[2] # 'C_mis'

            # Rewrite so that it only counts the latest one.
            if site not in sites[gen]:
                sites[gen][site] = {}
            if mut in sites[gen][site]:
                sites[gen][site][mut] += 1
            else:
                sites[gen][site][mut] = 1


"""

"""
# Count number of each mutation per generation
with open("../sim-data/neut-sites.tsv", "r") as neut_sites:
    sim_data = list(csv.reader(neut_sites, delimiter="\t"))
    
    # Find size of generation automatically
    gen_size = 0
    for indiv in sim_data:
        if int(indiv[1]) == 1:
            gen_size += 1
        else:
            break
    
    # Count number of unique mutations
    mutation_types = {}
    for indiv in sim_data:
        # Skip unmutated sequences
        if indiv[-1] == '':
            continue
        
        gen = int(indiv[1])
        if gen not in mutation_types:
            mutation_types[gen] = {}

        mutations = indiv[-1].split("|")
        for m in mutations:
            if m in mutation_types[gen]:
                mutation_types[gen][m] += 1
            else:
                mutation_types[gen][m] = 1

#print(mutation_types)
# mutation_types: for each generation, count number of each unique mutation.
#{1: {1532_C_syn: #, 52_A_non: #}, 2: {} }


print(mutation_types[313])
# Unique mutations list
#keys = list(mutation_types[313].keys())
#print(keys)

# sites format: { site: {"C_syn": #, "A_mis", #} }
"""
"""
sites = {}
for key in keys:
    k = key.split("_")
    site = k[0]
    if site not in sites:
        sites[site] = {}
    mut = k[1]+"_"+k[2]
    if mut in sites[site]:
        sites[site][mut] += 1
    else:
        sites[site][mut] = 1

print(sites)
"""


"""
# total number of each mutation type per individual
with open("../sim-data/neut-sites.tsv", "r") as sites:
    sim_data = list(csv.reader(sites, delimiter="\t"))
    mutation_types = {}
    for indiv in sim_data:
        gen = int(indiv[1])
        if gen not in mutation_types:
            mutation_types[gen] = []
        mutation_types[gen].append({"mis": 0, "non": 0, "syn": 0})
        mutations = indiv[-1].split("|")
        for m in mutations:
            if m[-3:] == "mis":
                mutation_types[gen][-1]["mis"] += 1
            elif m[-3:] == "non":
                mutation_types[gen][-1]["non"] += 1
            elif m[-3:] == "syn":
                mutation_types[gen][-1]["syn"] += 1

print(mutation_types[500])
"""

"""
# total number of each mutation type per generation 
with open("../sim-data/neut-sites.tsv", "r") as sites:
    sim_data = list(csv.reader(sites, delimiter="\t"))
    mutation_types = {}
    for indiv in sim_data:
        gen = int(indiv[1])
        if gen not in mutation_types:
            mutation_types[gen] = {"mis": 0, "non": 0, "syn": 0}
        mutations = indiv[-1].split("|")
        for m in mutations:
            if m[-3:] == "mis":
                mutation_types[gen]["mis"] += 1
            elif m[-3:] == "non":
                mutation_types[gen]["non"] += 1
            elif m[-3:] == "syn":
                mutation_types[gen]["syn"] += 1

#print(mutation_types)
for n in range(500, 501):
    print(n, mutation_types[n])

print(mutation_types[500]["mis"]+mutation_types[500]["non"]+mutation_types[500]["syn"])
"""

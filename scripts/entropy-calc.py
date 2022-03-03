import csv
import math

neut = "../sim-data/neut-sites.tsv"
adpt = "../sim-data/adpt-sites.tsv"
bkg = "../sim-data/bkg-sites.tsv"
mix = "../sim-data/mix-sites.tsv"

def entropy(file):
    model = file.split("/")[-1].split("-")[0]
    with open(file, "r") as test_sites:
        sim_data = list(csv.reader(test_sites, delimiter="\t"))

    # Find size of generation automatically
    gen_size = 0
    for indiv in sim_data:
        if int(indiv[1]) == 1:
            gen_size += 1
        else:
            break

    # Find number of generations automatically
    num_generations = int(sim_data[-1][1])

    # Parse sites into dictionary
    # Format: {gen 1: {site_num: {mut: #, mut, #}, site_num: {mut: #, mut:#}, ... } 
    #              gen: {site_num: {mut: #, ...}, site_num: {...} }
    #              ...
    #              gen 500: {...}
    #             }

    # e.g., {1: {13090: {'T_mis': 1}, 24637: {'A_mis': 1, 'C_syn': 1}}, 2: {13090: {'T_mis': 2}, 24637: {'A_mis': 2, 'C_syn': 1}}}

    sites = {}
    for gen in range(1, num_generations+1):
        sites[gen] = {}

    # indiv: ['neut', '1', 'Seq_1_1', '1202_C_mis|10989_G_syn']
    for indiv in sim_data:
        # Skip unmutated sequences
        if indiv[-1] == '':
            continue

        gen = int(indiv[1])
        mutations = indiv[-1].split("|") # ['1202_C_mis', '10989_G_syn']
        
        # For each indiv, if a new mutation occurs on a previously mutated site, overwrite the previous mutation
        filtered_mutations = {}
        for mutation in mutations: 
            # mutation: '1202_C_mis'
            
            # Skip NA mutations (non-coding)
            if mutation[-2:] == "NA":
                continue
             
            m = mutation.split("_") # ['1202', 'C', 'mis']
            site = int(m[0]) # 1202
            mut = m[1] + "_" + m[2] # 'C_mis'
            filtered_mutations[site] = mut

        # Add the count to the sites dict
        for site in filtered_mutations:
            if site not in sites[gen]:
                sites[gen][site] = {}
                
            mut = filtered_mutations[site]
            if mut in sites[gen][site]:
                sites[gen][site][mut] += 1
            else:
                sites[gen][site][mut] = 1


    # Entropy calculation: for site i, H = summation( p(base) * log p(base)  for base = A, T, C, G)

    # gen_entropy: total entropy per generation
    # Format: {gen 1: 10 bits, 2: 5 bits, ...}
    gen_entropy = {}

    # Separate entropy into synonymous vs non-synonymous sites
    syn_entropy = {}
    mis_entropy = {}

    for gen in range(1, num_generations+1):
        total_entropy = 0 # Sum of entropy for all sites in a generation
        site_entropy = {} # Entropy for each mutated site
        
        syn_total_entropy = 0 # Sum entropy for syn mutations
        syn_site_entropy = {} # For each syn mutations
        
        mis_total_entropy = 0 # Sum entropy for non-syn mutations
        mis_site_entropy = {} # For non-syn mutations
        
        for site in sites[gen]: # Mutated sites 13090, 9403, ...
            site_entropy[site] = 0
            syn_site_entropy[site] = 0
            mis_site_entropy[site] = 0

            count_mutated = 0 # Total number of individuals with a mutation at that site
            syn_count = 0 # Count number of syn mutations at the site
            mis_count = 0 # Count number of mis mutations.

            for mut in sites[gen][site]: # C_syn, A_mis, ...
                freq = sites[gen][site][mut]
                count_mutated += freq
                if mut[-3:] == "mis":
                    mis_count += freq
                elif mut[-3:] == "syn":
                    syn_count += freq
            ancestral_count = gen_size - count_mutated
             # Calculate syn site entropy with population size only of ancestral + syn
            syn_pop = ancestral_count + syn_count
            mis_pop = ancestral_count + mis_count 
            
            for mut in sites[gen][site]: # C_syn, A_mis, ...
                freq = sites[gen][site][mut]
                count_mutated += freq
                
                # Calculate total entropy for the site
                site_entropy[site] += -1 * freq/gen_size * math.log(freq/gen_size, 2)
                
                # Calculate entropy with population size only of ancestral + mis or ancestral + syn
                if mut[-3:] == "mis":           
                    mis_site_entropy[site] += -1 * freq/mis_pop * math.log(freq/mis_pop, 2)
                elif mut[-3:] == "syn":
                    syn_site_entropy[site] += -1 * freq/syn_pop * math.log(freq/syn_pop, 2)

            if ancestral_count != 0:
                site_entropy[site] += -1 * ancestral_count/gen_size * math.log(ancestral_count/gen_size, 2)
                
                if mut[-3:] == "mis":           
                    mis_site_entropy[site] += -1 * ancestral_count/mis_pop * math.log(ancestral_count/mis_pop, 2)
                elif mut[-3:] == "syn":
                    syn_site_entropy[site] += -1 * ancestral_count/syn_pop * math.log(ancestral_count/syn_pop, 2)


        for s in site_entropy:
            total_entropy += site_entropy[s]
            syn_total_entropy += syn_site_entropy[s]
            mis_total_entropy += mis_site_entropy[s]
        
        gen_entropy[gen] = total_entropy
        syn_entropy[gen] = syn_total_entropy
        mis_entropy[gen] = mis_total_entropy

    print(gen_entropy[3])
    print(syn_entropy[3])
    print(mis_entropy[3])
    print(sites[3])
"""
    # Save output to file
    with open('../sim-data/'+model+'-entropy.tsv', 'w') as outfile:
        outfile.write("Generation\tEntropy\tType\tModel\n")
        for n in range(1, num_generations+1):
            outfile.write("{0}\t{1}\t{2}\t{3}\n".format(n, gen_entropy[n], "Total", model))
            outfile.write("{0}\t{1}\t{2}\t{3}\n".format(n, syn_entropy[n], "Synonymous", model))
            outfile.write("{0}\t{1}\t{2}\t{3}\n".format(n, mis_entropy[n], "Missense", model))
"""
"""
with open('../sim-data/entropy-neut1.tsv', 'w') as outfile:
    outfile.write("Generation\tSynonymous\tNonsynonymous\tEntropy\n")
    for n in range(1, num_generations+1):
        outfile.write("{0}\t{1}\t{2}\t{3}\n".format(n, syn_entropy[n], nsyn_entropy[n], gen_entropy[n]))  
"""
entropy(neut)
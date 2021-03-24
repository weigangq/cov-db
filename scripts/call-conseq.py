#!/usr/bin/env python
# Coding: utf-8

import logging
import sys
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import datetime

# Initialize parser
parser = argparse.ArgumentParser(
    description="Call consequence from dnadiff .snps file (-f) or from a genome position (-p), based on a GenBank genome file (-b)")

################################
# Arguments
##############################

# Arguments for input & outputs
parser.add_argument('-b', '--genbank',
                    help='GenBank file (from NCBI). Required.')

parser.add_argument('-f', '--snp_file',
                    help='SNP file from dnadiff')

parser.add_argument('-p', '--position', type = int,
                    help='genome posiiton (1-based)')

parser.add_argument('-a', '--altNT',
                    help='alternative nucleotide')


######################################
# Program settings
#########################################

args = parser.parse_args()
logging.basicConfig(level = logging.DEBUG)

if args.genbank is None:
    logging.info("Need a Genbank file as -b NC_045512.2.gb")
    sys.exit()

if args.snp_file is None and args.position is None:
    logging.info("Need EITHER a SNP file from DNADIFF as -f test.snps OR a position as -p 300")
    sys.exit()

if args.snp_file is not None and args.position is not None:
    logging.info("Need EITHER a SNP file from DNADIFF as -f test.snps OR a position as -p 300. Not both!")
    sys.exit()

ref_gb = SeqIO.read(args.genbank, 'genbank')
genomeSeq = list(str(ref_gb.seq)) # make a list
if args.position is not None and (args.position < 0 or args.position > len(genomeSeq)):
    logging.info("Position not in genome. Genome length: %s", len(genomeSeq))
    sys.exit()


################################
# Define all functions
##################################
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
    proteins = {}
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
                    proteins[id] = { 'id': id,
                                     'gene': feature.qualifiers['gene'][0],
                                     'location': feature.location, # compound for nsp12
                                     'product': feature.qualifiers['product'][0],
                                     'cds_seq': cds,
                                     'pep_seq': pep,
                                     'sample_synon': 0,
                                     'sample_missense': 0
                                   }
    return proteins

def get_codons(cds_location):
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
    positions = list(cds_location['location'])
    for index, position in enumerate(positions):
        if index % 3 == 0:
            second_codon_position = positions[index + 1]
            third_codon_position = positions[index + 2]
            codons.append(
                (position, second_codon_position, third_codon_position))
    return codons

def position_info(genbank):
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
#    genes = gene_locations(genbank)
    position_data = {}
#    position_data['seq'] = genbank.seq
    seen_position = {}

    for geneId in cdsObj: # for each gene
        # Get codon locations
        # print(gene['location'])
        gene = cdsObj[geneId]
        codons = get_codons(gene)
        # geneId = gene['id']
        geneProd = gene['product']
        geneName = gene['gene']
        # print(codons)
        # continue
        for index, codon_location in enumerate(codons): # for each codon; the location object handles slippage site correctly: site 13467 in two codons "AAC" & "CGG"
            # Translated amino acid
            amino_acid = gene['pep_seq'][index]
            codon = str(genbank[codon_location[0]:codon_location[2] + 1].seq)
            # Base position in a codon
            for position in codon_location: # for each codon position
#                validate slippage position:
#                if position > 13465 and position < 13470:
#                    print(position, end = "\t")
#                    print(codon)
                if position in seen_position: # assuming one position one frame (ribosome slippage site [13468, or 13467 index] actually serves in two reading frames; 2nd reading frame is ignored); also avoid overlapping genes (skipped in the get_locations function)
                    continue
                seen_position[position] = 1
                # Nucleotide at that position
                nucleotide = genbank[position]
                # Create a codon dictionary
#                codon_object = {
#                    'codon': codon, 'amino_acid': amino_acid, 'location': codon_location}
                position_data[position] = {
                    'nt': nucleotide,
                    'aa': amino_acid,
                    'location': codon_location,
                    'codon': codon,
                    'geneId': geneId,
                    'geneName':  geneName,
                    'geneProduct': geneProd
                }
#               position_data[position]['codons'].append(codon_object)
    return position_data

def conseq (mutation_site):
    '''
        Parameters:
            mutation_site (int): A position where mutation occurs in
                                 an individual's sequence.
        Returns:
            individual (dict): Individual with a changed fitness value
                               if there was a nonsynonymous mutation.
    '''
    site = {
        'mut_site': mutation_site,
        'nt_ref': genomeSeq[mutation_site],
        'nt_post': ['NA', 'NA', 'NA'],
        'gene_id': 'NA',
        'codon_ref': 'NA',
        'codon_post': ['NA', 'NA', 'NA'],
        'aa_ref': 'NA',
        'aa_post': ['NA', 'NA', 'NA'],
        'conseq': ['igs', 'igs', 'igs'],
        'codonPos': 0
    }

    # Return individual as it is if a mutation occurs at a non-coding site
    if mutation_site in posInfo:
        position_info = posInfo[mutation_site]
        gene = cdsObj[position_info['geneId']]
        index = position_info['location'].index(mutation_site) # What part of codon
        site['gene_id'] = position_info['geneId']
        site['aa_ref'] = position_info['aa']
        site['codon_ref'] = position_info['codon']
        nt_ref = position_info['nt']
        site['nt_ref'] = nt_ref
        site['codonPos'] = index

        post_nts = []
        new_codons = []
        if index == 0: # First base of codon
            base2 = mutation_site + 1
            base3 = mutation_site + 2
            for nt in bases:
                if nt != nt_ref:
                    post_nts.append(nt)
                    new_codons.append(nt + genomeSeq[base2] + genomeSeq[base3])
        elif index == 1: # Second base of codon
            base1 = mutation_site - 1
            base3 = mutation_site + 1
            for nt in bases:
                if nt != nt_ref:
                    post_nts.append(nt)
                    new_codons.append(genomeSeq[base1] + nt + genomeSeq[base3])
        else: # Third base of codon
            base1 = mutation_site - 2
            base2 = mutation_site - 1
            for nt in bases:
                if nt != nt_ref:
                    post_nts.append(nt)
                    new_codons.append(genomeSeq[base1] + genomeSeq[base2] + nt)

        # BioPython Seq object
        post_codons = []
        post_aas = []
        conseqs = []
        for new_codon in new_codons:
            new_codon = Seq(new_codon)
            post_codons.append(str(new_codon))
            new_amino_acid = str(new_codon.translate()) 
            post_aas.append(new_amino_acid)

            # obtain fitness with respect to ref genome
            if new_amino_acid == position_info['aa']: # synonymous (including stop)
                conseqs.append('syn')
            elif new_amino_acid == '*' or position_info['aa'] == '*': # sense <=> stop
                conseqs.append('non')
            else: # nonsyn; ind fitness is cumulative site/codon fitness
                conseqs.append('mis')
        site['aa_post'] = post_aas
        site['codon_post'] = post_codons
        site['conseq'] = conseqs
        site['nt_post'] = post_nts
    return site

def outputSNPs():
    gb = args.genbank

    for site in snps:
        pos = site['pos']
        refNT = site['refNT']
        altNT = site['altNT']
        
        if site['info'] is not None:
            geneId = site['info']['geneId']
            geneName = site['info']['geneName']
            geneProd = site['info']['geneProduct']
            codonIdx = site['mut']['codonPos']
            codonRef = site['mut']['codon_ref']
            aaRef = site['mut']['aa_ref']
            altCodon = ''
            listIdx = 0
            for i in range(0,3):
                codon = site['mut']['codon_post'][i]
                if codon[codonIdx] == altNT:
                    altCodon = codon
                    listIdx = i
                    altAA = site['mut']['aa_post'][listIdx]
                    conseq = site['mut']['conseq'][listIdx]
            print("\t".join([gb, str(pos), refNT, altNT, geneId, geneName, geneProd, codonRef, altCodon, aaRef, altAA, conseq]))
        else:
            print("\t".join([gb, str(pos), refNT, altNT, 'igs', 'igs', 'igs', 'NA', 'NA', 'NA', 'NA', 'igs']))

def outputPos():
    gb = args.genbank

    site = snps.pop()
    pos = site['pos']
    refNT = site['refNT']
    altNT = site['altNT']
        
    if site['info'] is not None: # an ORF site
        geneId = site['info']['geneId']
        geneName = site['info']['geneName']
        geneProd = site['info']['geneProduct']
        codonIdx = site['mut']['codonPos']
        codonRef = site['mut']['codon_ref']
        aaRef = site['mut']['aa_ref']

        if altNT is None:
            for i in range(0,3):
                altNT = site['mut']['nt_post'][i]
                altCodon = site['mut']['codon_post'][i]
                altAA = site['mut']['aa_post'][i]
                conseq = site['mut']['conseq'][i]
                print("\t".join([gb, str(pos), refNT, altNT, geneId, geneName, geneProd, codonRef, altCodon, aaRef, altAA, conseq]))

        else:
            if refNT == altNT:
                logging.info("Error: refNT and altNT are the same")
            else:
                altCodon = ''
                listIdx = 0
                for i in range(0,3):
                    codon = site['mut']['codon_post'][i]
                    if codon[codonIdx] == altNT:
                        altCodon = codon
                        listIdx = i
                        altAA = site['mut']['aa_post'][listIdx]
                        conseq = site['mut']['conseq'][listIdx]
                print("\t".join([gb, str(pos), refNT, altNT, geneId, geneName, geneProd, codonRef, altCodon, aaRef, altAA, conseq]))

    else: # an IGS site
        if altNT is None:
            for base in bases:
                if refNT == base:
                    continue
                print("\t".join([gb, str(pos), refNT, base, 'igs', 'igs', 'igs', 'NA', 'NA', 'NA', 'NA', 'igs']))
        else:
            if refNT == altNT:
                logging.info("Error: refNT and altNT are the same")
            else:
                print("\t".join([gb, str(pos), refNT, altNT, 'igs', 'igs', 'igs', 'NA', 'NA', 'NA', 'NA', 'igs']))

logging.info("Start timestamp: %s", datetime.datetime.now())
cdsObj = gene_locations(ref_gb)
posInfo = position_info(ref_gb) # Note: zero-based pos index
bases = ['A', 'T', 'C', 'G']

snps = []
if args.snp_file is not None:
    with open(args.snp_file, "r") as fh:
        lines = fh.readlines()
        for line in lines:
            data = line.split()
            pos = int(data[0])
            gIdx = pos - 1
            if data[1] == '.' or data[2] == '.':
                continue
            if gIdx in posInfo:
                snps.append({'pos': pos,
                             'refNT': data[1],
                             'altNT': data[2],
                             'info': posInfo[gIdx],
                             'mut': conseq(gIdx)
                         })
            else:
                snps.append({'pos': pos,
                             'refNT': data[1],
                             'altNT': data[2],
                             'info': None,
                             'mut': None
                         })
    outputSNPs()

if args.position is not None:
    gIdx = args.position - 1
    pos = args.position

    if gIdx in posInfo:
        snps.append({'pos': pos,
                     'refNT': genomeSeq[gIdx],
                     'altNT': args.altNT if args.altNT is not None else None,
                     'info': posInfo[gIdx],
                     'mut': conseq(gIdx)
                 })
    else:
        snps.append({'pos': pos,
                     'refNT': genomeSeq[gIdx],
                     'altNT': args.altNT if args.altNT is not None else None,
                     'info': None,
                     'mut': None
                 })

    outputPos()

#print(snps[:10])
logging.info("End timestamp: %s", datetime.datetime.now())
sys.exit()

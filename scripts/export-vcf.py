#!/usr/bin/env python
# coding: utf-8

from config import config
import re
import numpy as np
import psycopg2
import argparse
import datetime
import sys
import logging
import vcfpy

# Initialize parser
parser = argparse.ArgumentParser(
    description='Export VCF file from two inputs (iso & var) from the export-maf.py')

# Add arguments
parser.add_argument('-i', '--iso', required = True,
                    help = 'provide EPI list in a file')

parser.add_argument('-v', '--var', required = True,
                    help = 'provide variant list in a file')

parser.add_argument('-a', '--syn', action = 'store_true',
                    help = 'synonymous SNPs only')

parser.add_argument('-b', '--nonsyn', action = 'store_true',
                    help = 'Nonsynonymous SNPs only')

parser.add_argument('-s', '--snp', action = 'store_true',
                    help = 'SNPs only')

parser.add_argument('-o', '--vcf', default = 'cov.vcf',
                    help = 'Name of vcf output file (default: cov.vcf)')

parser.add_argument('-d', '--diploid', action = 'store_true',
                    help = 'Diploid Genotype output, for e.g., vcftools --site-pi')

args = parser.parse_args()

logging.basicConfig(level = logging.DEBUG)
refEPI = 'EPI_ISL_406030'
batEPI = 'EPI_ISL_402131'
logging.info("Start timestamp: %s", datetime.datetime.now())

# From https://www.postgresqltutorial.com/postgresql-python/connect/
dbParams = config()
conn = psycopg2.connect(**dbParams)
cur = conn.cursor()

############################################
# get all samples (except the ref and bat)
############################################
isoEPIs = []
isoFile = open(args.iso, "r")
lines = isoFile.readlines()
for line in lines:
    data = line.split()
    isoEPIs.append(data[0])
isoEPIs.sort()
isoFile.close()
logging.info("Isolates collected and sorted: n = %s", len(isoEPIs))


#######################################################
# variant info & create varID (freq > 0.1%)
########################################################
varFile = open(args.var, "r")
snpInfo = {}
indelInfo = {}
lines = varFile.readlines()
varCt = 0
multCt = 0
snpCt = 0
delCt = 0
insCt = 0
seenID = {}
for line in lines:
    varCt += 1
    data = line.split()
    site = int(data[1])
    varType = data[0]
    varID = data[2] # e.g., cv-T29

    if args.snp and varType != 'SNP':
        continue

    if args.syn and data[9] != 'synonymous':
        continue

    if args.nonsyn and data[9] != 'missense':
        continue

    seenID[varID] = 1
    if varType == 'SNP':
        snpCt += 1
        isCoding = 0 if data[9] == 'noncoding'  else 1
        if site in snpInfo:
            snpInfo[site]['varID'].append(varID)
            snpInfo[site]['altNT'].append(data[8])
            snpInfo[site]['conseq'].append(data[9])
            multCt += 1
        else:
            snpInfo[site] = {
                'site': site,
                'varID': [varID],
                'varType': varType,
                'altNT': [data[8]], # could be multiple altNTs
                'refNT': data[3],
                'locus': data[4],
                'codonRef': data[5] if isCoding else 'NA',
                'codonPos': str(data[6]) if isCoding else 'NA',
                'locPos': str(data[7]) if isCoding else 'NA',
                'conseq': [ data[9] ]
            }

    else: # indel
        if varType == 'DEL':
            delCt += 1
        else:
            insCt += 1

        indelInfo[varID] = {
            'site': site,
            'varType': varType,
            'varID': [varID],
            'altNT': [data[8]],
            'refNT': data[3],
            'locus': data[4],
            'codonRef': 'NA',
            'codonPos': 'NA',
            'locPos': str(data[7]),
            'conseq': [ 'NA' ]
        }
varFile.close()

logging.info("Variants read: n = %s", varCt)
logging.info("Number of SNPs: n = %s", snpCt)
logging.info("Number of multi-state SNPs: n = %s", multCt)

sitesSNPs = list(snpInfo.keys())
sitesINDELs = list(indelInfo.keys())
snpDict = dict(sorted(snpInfo.items(), key=lambda item: item[1]['site']))
indelDict = dict(sorted(indelInfo.items(), key=lambda item: item[1]['site']))

logging.info("high freq SNPs info collected and sorted: n = %s", len(sitesSNPs))
logging.info("high freq INDELs info collected and sorted: n = %s", len(sitesINDELs))
logging.info("changed skipped: n = %s",  varCt - snpCt - delCt - insCt)

##############################
# get genotype for samples
###############################
genoSample = {}
sampleCt = 1

tup_acc = tuple(isoEPIs)
par_acc = {'l': tup_acc}
cur.execute('select acc, chg from human_anno a, hap_var b where a.hid = b.hid and acc in %(l)s', par_acc)
hap = cur.fetchall()
genoSample = {}
# initialize double dict !!!!!!!
for acc in isoEPIs:
    genoSample[acc] = {}

for line in hap:
    acc = line[0]
    change = line[1]
    varID = "cv-" + change
    if varID not in seenID:
        continue
    if re.match("[ATCG]\d+", change): # a SNP, with id "T1234"
        site = int(change[1:])
        alt = change[0]
        genoSample[acc][site] = alt
# note mixed keys for genoSample: str for indels (to avoid conflict)
    elif re.search("\d+-\d+", change): # a deletion at n-1, with id "123-45"
        x = change.split("-")
        site = int(x[0]) - 1
        genoSample[acc][varID] = indelInfo[varID]['altNT'][0] # 'N'
    else:  # insertion, e.g., "234_ATG"
        x = change.split("_")
        site = int(x[0]) - 1
        genoSample[acc][varID] = indelInfo[varID]['altNT'][0] # 'NAAAA'
#    logging.info("Genotypes collected for sample %s at site %s with %s is %s", acc, site, varID, genoSample[acc][site])
logging.info("genotypes collected for : n = %s isolates, at %s", len(list(genoSample.keys())), datetime.datetime.now())
filteredEPIs = list(genoSample.keys())
filteredEPIs.sort() # filter out EPIs contains only low-freq variants
#for acc in ['ISL_700228', 'ISL_539719', 'ISL_539706', 'ISL_539708']:
#    x = genoSample[acc]
#    print(acc, x)
#sys.exit()
#####################
# construct VCF records
#############################
timePre = datetime.datetime.now()
varCt = 0

header = vcfpy.Header(
    samples = vcfpy.SamplesInfos(filteredEPIs),
    lines = [
        vcfpy.HeaderLine('fileformat', 'VCFv4.0'),
        vcfpy.HeaderLine('fileDate', str(datetime.datetime.now())),
        vcfpy.HeaderLine('source', parser.prog),
        vcfpy.ContigHeaderLine('contig', '<ID=String,Length=Integer>', {'ID': 'EPI_ISL_406030', 'length': 29903}),
        vcfpy.InfoHeaderLine('INFO', '<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">', {'ID': 'NS', 'Number': 1, 'Type': 'Integer', 'Description': 'Number of Samples With Data'}),
        vcfpy.FormatHeaderLine('FORMAT', '<ID=GT,Number=1,Type=String,Description="Genotype">', {'ID': 'GT', 'Number': 1, 'Type': 'String', 'Description': 'Genotype'})
    ]
)

with vcfpy.Writer.from_path(args.vcf, header) as writer:
    for site in snpDict:
        geno = {}
        genoCalls = []
        refNT = snpInfo[site]['refNT']
        altNTs = snpInfo[site]['altNT']
        geno[refNT] = 0
        altCount = 1
        subs = [] # substitutions
        for alt in altNTs:
            geno[alt] = altCount
            altCount += 1
            subs.append(vcfpy.Substitution(type_ = 'SNV', value = alt))

        for acc in filteredEPIs:
            allele = 0
            if site in genoSample[acc]: # is mutated
 #               print(genoSample[acc][site])
                if genoSample[acc][site] in geno: # alt is valid
                    allele = geno[genoSample[acc][site]]
#                    logging.info("alt assigned for %s at %s: %s", acc, site, allele)
                else: # alt is singleton/discarded
                    logging.warning("alt is singleton for %s at %s: assign ref allele", acc, site)
#            else:
#                logging.info("ref alleles assigned for %s at %s", acc, site)
            gt = str(allele) + "|" + str(allele) if args.diploid else str(allele)
            sampleCall = vcfpy.Call(
                sample = acc,
                data = {'GT': gt }, # has to be string; diploid
#                data = {'GT': str(allele) }, # has to be string
                site = site
            )
            genoCalls.append(sampleCall)

        record = vcfpy.Record(
            CHROM = refEPI,
            POS = site,
            ID = snpInfo[site]['varID'],
            REF = snpInfo[site]['refNT'],
            ALT = subs,
            QUAL = None,
            FILTER = [], # PASS
            INFO = {},  # consequence calls, locus, etc; a dict
            FORMAT = ['GT'], # a list
            calls = genoCalls
        )
        varCt += 1
        writer.write_record(record)
    logging.info("SNPs records written to file: n = %s at %s", len(sitesSNPs), datetime.datetime.now())

    for change in indelDict: # change is "cv-" varID
        geno = {}
        genoCalls = []
        refNT = indelInfo[change]['refNT'] # 'NAAAAA'
        altNTs = indelInfo[change]['altNT'] # [ 'N' ]
        subs = [] # substitutions
        geno[refNT] = 0
        altCount = 1
        varType = indelInfo[change]['varType']
        for alt in altNTs:
            geno[alt] = altCount
            altCount += 1
            subs.append(vcfpy.Substitution(type_ = varType, value = alt))

        for acc in filteredEPIs:
            allele = 0
            if change in genoSample[acc]: # is mutated
                if genoSample[acc][change] in geno: # alt is valid
                    allele = geno[genoSample[acc][change]]
                else: # alt is singleton/discarded
                    logging.warning("alt is singleton for %s at %s: assign ref allele", acc, site)

            gt = str(allele) + "|" + str(allele) if args.diploid else str(allele)

            sampleCall = vcfpy.Call(
                sample = acc,
                data = {'GT': gt }, # has to be string, diploid
#                data = {'GT': str(allele) }, # has to be string
                site = site
            )
            genoCalls.append(sampleCall)

        record = vcfpy.Record(
            CHROM = refEPI,
            POS = indelInfo[change]['site'],
            ID = indelInfo[change]['varID'],
            REF = indelInfo[change]['refNT'],
            ALT = subs,
            QUAL = None,
            FILTER = [], # PASS
            INFO = {},  # consequence calls, locus, etc; a dict
            FORMAT = ['GT'], # a list
            calls = genoCalls
        )
        varCt += 1
        writer.write_record(record)
    logging.info("INDELs records written to file: n = %s at %s", len(sitesINDELs), datetime.datetime.now())

logging.info("End timestamp: %s", datetime.datetime.now())
sys.exit();

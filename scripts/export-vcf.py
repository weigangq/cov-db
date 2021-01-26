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
    description='Export isolate IDs from database')

# Add arguments
parser.add_argument('-i', '--iso', required = True,
                    help = 'provide EPI list in a file')

parser.add_argument('-v', '--var', required = True,
                    help = 'provide variant list in a file')

parser.add_argument('-o', '--vcf', required = True,
                    help = 'Name of vcf output file')

args = parser.parse_args()

#logging.basicConfig(filename = 'export-vcf.log', filemode = 'w', level = logging.DEBUG)
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
#print(isoEPIs[:9])
#sys.exit()


#######################################################
# variant info & create varID (freq > 0.1%)
########################################################
varFile = open(args.var, "r")
snpInfo = {}
delInfo = {}
insInfo = {}
lines = varFile.readlines()
changes = {} # only the high frequency vars
varCt = 0
multCt = 0
snpCt = 0
snpList = []
delList = []
insList = []
for line in lines:
    varCt += 1
    data = line.split()
    site = int(data[0])
    change = data[1]
    changes[change] = 1
#    print(change, end = '\t')
    if re.match("[ATCG]\d+", change): # a SNP, combine multi-allelic
        snpList.append(change)
    elif re.search("\d+-\d+", change): # deletion, e.g., "4560-45"
        delList.append(change)
    else:
        insList.append(change)
logging.info("Variants read: n = %s", varCt)

# collect SNPs        
par_snp = { 'l': tuple(snpList) }  
cur.execute("select * from cv_snp where concat(alt, site) in %(l)s", par_snp) # PK: site + alt
snp = cur.fetchall()
#if len(snp) == 0:
for dataSNP in snp:
    site = dataSNP[0]
    change = dataSNP[1] + str(dataSNP[0])
    snpCt += 1
    isCoding = 0 if dataSNP[3] is None else 1
    conseq = ''
    if isCoding:
        if dataSNP[4] == dataSNP[5]:
            conseq = 'synonymous'
        else:
            conseq = 'missense'
    else:
        conseq = 'noncoding'
            
    if site in snpInfo:
        snpInfo[site]['varID'].append("cv" + "-" + change)
        snpInfo[site]['altNT'].append(dataSNP[1])
        snpInfo[site]['conseq'].append(conseq)
        multCt += 1
    else:
        snpInfo[site] = {
            'site': site,
            'varID': ["cv" + "-" + change],
            'vartype': 'SNP',
            'altNT': [dataSNP[1]], # could be multiple altNTs
            'refNT': dataSNP[2][dataSNP[3]-1] if isCoding else dataSNP[2], 
            'locus': dataSNP[6],
            'codonRef': dataSNP[2] if isCoding else 'NA',
            'codonPos': str(dataSNP[3]) if isCoding else 'NA',
            'locPos': str(dataSNP[7]) if isCoding else 'NA',
            'conseq': [ conseq ]
        }

logging.info("Number of SNPs: n = %s", snpCt)
logging.info("Number of multi-state SNPs: n = %s", multCt)

# deletion
delCt = 0
par_del = { 'l': tuple(delList) }
 #   elif re.search("\d+-\d+", change): # deletion, e.g., "4560-45"
 #       print("Not matching")
cur.execute("select * from cv_del where del in %(l)s", par_del) 
DEL = cur.fetchall()
for dataDEL in DEL:
    delCt += 1
    change = dataDEL[0]
    site = int(change.split("-")[0])
    delInfo[change] = {
        'site': site - 1,
        'vartype': 'DEL',
        'varID': ["cv" + "-" + change],
        'altNT': ['N'], # A 
        'refNT': 'N' + dataDEL[3], # ATTTTTTTTTTTTTTTTT
        'locus': dataDEL[1],
        'codonRef': 'NA',
        'codonPos': 'NA',
        'locPos': str(dataDEL[2]),
        'conseq': [ 'NA' ]
    }

insCt = 0    
par_ins = { 'l': tuple(insList) } 
cur.execute("select * from cv_ins where ins in %(l)s", par_ins) 
INS = cur.fetchall()
for dataINS in INS:
    insCt += 1
    ins = dataINS[0].split("_")
    site = int(ins[0])
    alt = ins[1]
    change = dataINS[0]
    refBase = dataINS[4][dataINS[3]-1] if dataINS[3] else dataINSp[4]
    insInfo[change] = {
        'site': site - 1,
        'vartype': 'INS',
        'varID': ["cv" + "-" + change],
        'altNT': [refBase + alt], 
        'refNT': refBase, 
        'locus': dataINS[1],
        'codonRef': 'NA',
        'codonPos': 'NA',
        'locPos': str(dataINS[2]),
        'conseq': [ 'NA' ]
    }
    
sitesSNPs = list(snpInfo.keys())
#sitesSNPs.sort()
sitesDELs = list(delInfo.keys())
sitesINSs = list(insInfo.keys())
#sitesDELs.sort()
#print(delInfo)

#print(snpInfo)
# sort by site (for printing)
snpDict = dict(sorted(snpInfo.items(), key=lambda item: item[1]['site']))
delDict = dict(sorted(delInfo.items(), key=lambda item: item[1]['site']))
insDict = dict(sorted(insInfo.items(), key=lambda item: item[1]['site']))
#print(snpDict)
#sys.exit()

varOut = open("var.tsv", "w")
for site in snpInfo.keys(): # int key
    varOut.write(snpInfo[site]['vartype'] + "\t"+ str(snpInfo[site]['site']) + "\t" + ",".join(snpInfo[site]['varID']) + "\t" + snpInfo[site]['refNT'] + "\t" + snpInfo[site]['locus'] + "\t" + snpInfo[site]['codonRef'] + "\t" + snpInfo[site]['codonPos'] + "\t" + snpInfo[site]['locPos'] + "\t" + ",".join(snpInfo[site]['altNT']) + "\t" + ",".join(snpInfo[site]['conseq']) + "\n")

for site in delInfo.keys(): # str key!!
    varOut.write(delInfo[site]['vartype'] + "\t"+ str(delInfo[site]['site']) + "\t" + ",".join(delInfo[site]['varID']) + "\t" + delInfo[site]['refNT'] + "\t" + delInfo[site]['locus'] + "\t" + delInfo[site]['codonRef'] + "\t" + delInfo[site]['codonPos'] + "\t" + delInfo[site]['locPos'] + "\t" + ",".join(delInfo[site]['altNT']) + "\t" + ",".join(delInfo[site]['conseq']) + "\n")

for site in insInfo.keys(): # str key!!
    varOut.write(insInfo[site]['vartype'] + "\t"+ str(insInfo[site]['site']) + "\t" + ",".join(insInfo[site]['varID']) + "\t" + insInfo[site]['refNT'] + "\t" + insInfo[site]['locus'] + "\t" + insInfo[site]['codonRef'] + "\t" + insInfo[site]['codonPos'] + "\t" + insInfo[site]['locPos'] + "\t" + ",".join(insInfo[site]['altNT']) + "\t" + ",".join(insInfo[site]['conseq']) + "\n")
   
varOut.close()
varFile.close()
logging.info("high freq SNP info collected and sorted: n = %s", len(sitesSNPs))
logging.info("high freq DEL info collected and sorted: n = %s", len(sitesDELs))
logging.info("high freq INS info collected and sorted: n = %s", len(sitesINSs))
logging.info("Var file created: var.tsv, at %s", datetime.datetime.now())
logging.info("change not found in database: n = %s",  varCt - snpCt - delCt - insCt)

#sys.exit()
##############################
# get genotype for samples
###############################
genoSample = {}
sampleCt = 1

tup_acc = tuple(isoEPIs)
par_acc = {'l': tup_acc}
#tup_chg = tuple(changes)
#par_chg = {'m': tup_chg}
cur.execute('select acc, chg from acc_hap a, hap_chg b where a.hid = b.hid and acc in %(l)s', par_acc)
hap = cur.fetchall()

genoSample = {} 
for line in hap:
    acc = line[0]
    change = line[1]
    if change not in changes: # didn't make it freq >= 0.1% cutoff
        continue
    genoSample[acc] = {}
    if re.match("[ATCG]\d+", change): # a SNP, with id "cv1234"
        site = int(change[1:])
        alt = change[0]
        genoSample[acc][site] = alt
    elif re.search("\d+-\d+", change): # a deletion at n-1, with id "cv123-45"
        x = change.split("-")
        site = int(x[0])
        genoSample[acc][site] = delInfo[change]['altNT'][0] # 'N'
    else:  # insertion, e.g., "234_ATG"
        x = change.split("_")
        site = int(x[0])
        genoSample[acc][site] = insInfo[change]['altNT'][0] # 'NAAAA'
#    logging.info("Genotypes collected for sample %s at %s", acc, sampleCt)
#    sampleCt += 1
#print(genoSample)
logging.info("genotypes collected for : n = %s isolates, at %s", len(list(genoSample.keys())), datetime.datetime.now())
filteredEPIs = list(genoSample.keys())
filteredEPIs.sort() # filter out EPIs contains only low-freq variants
#print(filteredEPIs)
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
                if genoSample[acc][site] in geno: # alt is valid
                    allele = geno[genoSample[acc][site]]
#                    logging.info("alt assigned for %s at %s: %s", acc, site, allele)
                else: # alt is singleton/discarded
                    logging.warning("alt is singleton for %s at %s: assign ref allele", acc, site)

            sampleCall = vcfpy.Call(
                sample = acc,
                data = {'GT': str(allele)}, # has to be string
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
        
    for change in delDict:
        geno = {}
        genoCalls = []
        refNT = delInfo[change]['refNT'] # 'NAAAAA'
        altNTs = delInfo[change]['altNT'] # [ 'N' ]
        subs = [] # substitutions
        geno[refNT] = 0
        altCount = 1
        for alt in altNTs:
            geno[alt] = altCount
            altCount += 1
            subs.append(vcfpy.Substitution(type_ = 'DEL', value = alt))

        for acc in filteredEPIs:
            allele = 0
            if site in genoSample[acc]: # is mutated
                if genoSample[acc][site] in geno: # alt is valid
                    allele = geno[genoSample[acc][site]]
                else: # alt is singleton/discarded
                    logging.warning("alt is singleton for %s at %s: assign ref allele", acc, site)

            sampleCall = vcfpy.Call(
                sample = acc,
                data = {'GT': str(allele)}, # has to be string
                site = site - 1
            )
            genoCalls.append(sampleCall)
            
        record = vcfpy.Record(
            CHROM = refEPI, 
            POS = delInfo[change]['site'], 
            ID = delInfo[change]['varID'], 
            REF = delInfo[change]['refNT'], 
            ALT = subs,
            QUAL = None, 
            FILTER = [], # PASS
            INFO = {},  # consequence calls, locus, etc; a dict
            FORMAT = ['GT'], # a list
            calls = genoCalls
        )
        varCt += 1
        writer.write_record(record)
    logging.info("DELs records written to file: n = %s at %s", len(sitesDELs), datetime.datetime.now())
    
    for change in insDict:
        geno = {}
        genoCalls = []
        refNT = insInfo[change]['refNT'] # 'A'
        altNTs = insInfo[change]['altNT'] # [ 'AG' ]
        subs = [] # substitutions
        geno[refNT] = 0
        altCount = 1
        for alt in altNTs:
            geno[alt] = altCount
            altCount += 1
            subs.append(vcfpy.Substitution(type_ = 'INS', value = alt))

        for acc in filteredEPIs:
            allele = 0
            if site in genoSample[acc]: # is mutated
                if genoSample[acc][site] in geno: # alt is valid
                    allele = geno[genoSample[acc][site]]
                else: # alt is singleton/discarded
                    logging.warning("alt is singleton for %s at %s: assign ref allele", acc, site)

            sampleCall = vcfpy.Call(
                sample = acc,
                data = {'GT': str(allele)}, # has to be string
                site = site - 1
            )
            genoCalls.append(sampleCall)
            
        record = vcfpy.Record(
            CHROM = refEPI, 
            POS = insInfo[change]['site'], 
            ID = insInfo[change]['varID'], 
            REF = insInfo[change]['refNT'], 
            ALT = subs,
            QUAL = None, 
            FILTER = [], # PASS
            INFO = {},  # consequence calls, locus, etc; a dict
            FORMAT = ['GT'], # a list
            calls = genoCalls
        )
        varCt += 1
        writer.write_record(record)
    logging.info("INSs records written to file: n = %s at %s", len(sitesINSs), datetime.datetime.now())       
 
logging.info("End timestamp: %s", datetime.datetime.now())
sys.exit();






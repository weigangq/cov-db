#!/usr/bin/env python
# coding: utf-8

# output iso & var files from DB
#
from config import config
import re
import psycopg2
import argparse
import datetime
import sys
import logging

if len(sys.argv) != 2:
    print("Usage: snp-freq.py country_name", file == sys.stderr)
    sys.exit()

freqCut = 0.001    
refEPI = 'EPI_ISL_406030'
batEPI = 'EPI_ISL_402131'
logging.basicConfig(level = logging.DEBUG)

# From https://www.postgresqltutorial.com/postgresql-python/connect/
dbParams = config()
conn = psycopg2.connect(**dbParams)
cur = conn.cursor()

############################################
# get all samples of a country (except the ref and bat)
############################################
logging.info("getting all samples of %s from the database...", sys.argv[1])

cur.execute("select acc, col_date, country, state from vhuman_anno where acc != %s and acc != %s and country = %s order by col_date", [refEPI, batEPI, sys.argv[1]])
accData = cur.fetchall()            
isoEPIs = [ iso[0] for iso in accData]

############################################
# define a function
############################################
def get_variant(changeList):
    snpInfo = {}
    delInfo = {}
    insInfo = {}
    varChange = {} # only the high frequency vars
    varCt = 0
    multCt = 0
    lowFreq = 0
    snpCt = 0
    snpList = []
    delList = []
    insList = []
    seenSite = {}
    for change in changeList:
        if varFreq[change] < freqCut:
            lowFreq += 1
            continue
        varCt += 1
        varChange[change] = 1
        if re.match("[ATCG]\d+", change): # a SNP, combine multi-allelic
            snpList.append(change)
        elif re.search("\d+-\d+", change): # deletion, e.g., "4560-45"
            delList.append(change)
        else: # insertion, e.g., "4560_AGT"
            insList.append(change)
    logging.info("Variants read: n = %s", varCt)

    # collect SNPs 
    par_snp = { 'l': tuple(snpList) }  
    cur.execute("select * from cv_snp where concat(alt, site) in %(l)s", par_snp) # PK: site + alt
    snp = cur.fetchall()
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

        if site in seenSite:
            multCt += 1
        else:
            seenSite[site] = 1

        snpInfo[change] = {
            'site': site,
            'varID': "cv" + "-" + change,
            'vartype': 'SNP',
            'altNT': dataSNP[1], # could be multiple altNTs
            'refNT': dataSNP[2][dataSNP[3]-1] if isCoding else dataSNP[2], 
            'locus': dataSNP[6],
            'codonRef': dataSNP[2] if isCoding else 'NA',
            'codonPos': str(dataSNP[3]) if isCoding else 'NA',
            'locPos': str(dataSNP[7]) if isCoding else 'NA',
            'conseq': conseq,
            'freq': str(round(varFreq[change],6)),
            'count': str(varAccCt[change]),
            'aaID': dataSNP[4] + str(int((dataSNP[7]-1)/3)+1) + dataSNP[5] if conseq == 'missense' else 'NA'
        }

    logging.info("Number of SNPs: n = %s", snpCt)
    logging.info("Number of multi-state SNPs: n = %s", multCt)

    # deletion
    delCt = 0
    par_del = { 'l': tuple(delList) }
    if len(delList):
        cur.execute("select * from cv_del where del in %(l)s", par_del) 
        DEL = cur.fetchall()
        for dataDEL in DEL:
            delCt += 1
            change = dataDEL[0]
            site = int(change.split("-")[0])
            delInfo[change] = {
                'site': site - 1,
                'vartype': 'DEL',
                'varID': "cv" + "-" + change,
                'altNT': 'N', # A 
                'refNT': 'N' + dataDEL[3], # ATTTTTTTTTTTTTTTTT
                'locus': dataDEL[1],
                'codonRef': 'NA',
                'codonPos': 'NA',
                'locPos': str(dataDEL[2]),
                'conseq': 'NA',
                'freq': str(round(varFreq[change],6)),
                'count': str(varAccCt[change]),
                'aaID': 'NA'
            }
            
    insCt = 0    
    par_ins = { 'l': tuple(insList) } 
    if len(insList):
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
                'varID': "cv" + "-" + change,
                'altNT': refBase + alt, 
                'refNT': refBase, 
                'locus': dataINS[1],
                'codonRef': 'NA',
                'codonPos': 'NA',
                'locPos': str(dataINS[2]),
                'conseq': 'NA',
                'freq': str(round(varFreq[change],6)),
                'count': str(varAccCt[change]),
                'aaID': 'NA'
            }

    sitesSNPs = list(snpInfo.keys())
    sitesDELs = list(delInfo.keys())
    sitesINSs = list(insInfo.keys())

    # sort by site (for printing)
    snpDict = dict(sorted(snpInfo.items(), key=lambda item: item[1]['site']))
    delDict = dict(sorted(delInfo.items(), key=lambda item: item[1]['site']))
    insDict = dict(sorted(insInfo.items(), key=lambda item: item[1]['site']))

    varOut = open('var_%s.tsv'% sys.argv[1], 'w')
    logging.info("Excluding low freq vars: n = %s", lowFreq)
    logging.info("exporting genetic changes with freq >= 0.1%...")
    for site in snpDict.keys(): # str key
        varOut.write(snpInfo[site]['vartype'] + "\t"+ str(snpInfo[site]['site']) + "\t" + snpInfo[site]['varID'] + "\t" + snpInfo[site]['refNT'] + "\t" + snpInfo[site]['locus'] + "\t" + snpInfo[site]['codonRef'] + "\t" + snpInfo[site]['codonPos'] + "\t" + snpInfo[site]['locPos'] + "\t" + snpInfo[site]['altNT'] + "\t" + snpInfo[site]['conseq'] + "\t"  + snpInfo[site]['freq'] + "\t" + snpInfo[site]['count'] + "\t" + snpInfo[site]['aaID'] +  "\n")

    for site in delDict.keys(): # str key!!
        varOut.write(delInfo[site]['vartype'] + "\t"+ str(delInfo[site]['site']) + "\t" + delInfo[site]['varID'] + "\t" + delInfo[site]['refNT'] + "\t" + delInfo[site]['locus'] + "\t" + delInfo[site]['codonRef'] + "\t" + delInfo[site]['codonPos'] + "\t" + delInfo[site]['locPos'] + "\t" + delInfo[site]['altNT'] + "\t" + delInfo[site]['conseq'] + "\t"  + delInfo[site]['freq'] + "\t" + delInfo[site]['count'] + "\t" + delInfo[site]['aaID'] +"\n")

    for site in insDict.keys(): # str key!!
        varOut.write(insInfo[site]['vartype'] + "\t"+ str(insInfo[site]['site']) + "\t" + insInfo[site]['varID'] + "\t" + insInfo[site]['refNT'] + "\t" + insInfo[site]['locus'] + "\t" + insInfo[site]['codonRef'] + "\t" + insInfo[site]['codonPos'] + "\t" + insInfo[site]['locPos'] + "\t" + insInfo[site]['altNT'] + "\t" + insInfo[site]['conseq'] + "\t"  + insInfo[site]['freq'] + "\t" + insInfo[site]['count'] + "\t" + insInfo[site]['aaID'] +"\n")

############################################
# get all genetic changes of each isolate
############################################
logging.info("getting all genetic changes: SNPs and indels...")
l = tuple(isoEPIs)
params = {'l': l}
cur.execute('select acc, chg from acc_hap a, hap_chg b where a.hid = b.hid and acc in %(l)s', params)
acc_snp = cur.fetchall()

allSamples = {}
var_count = {}
for line in acc_snp:
    acc = line[0]
    chg = line[1]
    if chg in allSamples:
        allSamples[chg].append(acc)
    else:
        allSamples[chg] = [acc]
 
    if acc in var_count:
        var_count[acc] += 1
    else:
        var_count[acc] = 1
    
changes = list(allSamples.keys())
logging.info("total genetic changes: n = %s", len(changes))

# get frequency & count for each variant
varFreq = {}
varAccCt = {}
for change in changes:
    iso = allSamples[change]
    freq = float(len(iso))/float(len(isoEPIs))
    varFreq[change] = freq
    varAccCt[change] = len(iso)

get_variant(changes)

# write iso file
acc_output = open('iso_%s.tsv'% sys.argv[1], 'w')
for iso in accData:
    iso = list(iso) # make it list, since tuple is not editable
    if iso[3] is None:
        iso[3] = 'NA'
    iso[1] = re.sub(r"^(\d{4}-\d{2})$", r"\1-15", iso[1])
    ct = var_count[iso[0]]
    acc_output.write("\t".join(iso) + "\t" + str(ct) + "\n")
acc_output.close()

logging.info("total isolates: %s", len(isoEPIs))
logging.info("Done!")
sys.exit



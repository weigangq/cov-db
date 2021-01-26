#!/usr/bin/env python
# coding: utf-8

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

# From https://www.postgresqltutorial.com/postgresql-python/connect/
dbParams = config()
conn = psycopg2.connect(**dbParams)
cur = conn.cursor()

############################################
# get all samples of a country (except the ref and bat)
############################################
print("getting all samples of %s from the database..." % sys.argv[1])

cur.execute("select a.acc, a.col_date, b.area, b.country from human_anno a, vgeo b where a.need is true and a.geo_id=b.geo_id and a.acc != %s and a.acc != %s and b.country = %s order by col_date;", [refEPI, batEPI, sys.argv[1]])
acc = cur.fetchall()            

acc_output = open('iso_%s.tsv'% sys.argv[1], 'w')
isoEPIs = []
for data in acc:
#    data = list(acc[i]) # make it list, since tuple is not editable
    isoEPIs.append(data[0])
    acc_output.write("\t".join(data) + "\n")
    #print("\t".join(data))
#isoEPIs.sort()
acc_output.close()
print("total isolates: %s" % len(isoEPIs))

############################################
# define a function
############################################
def get_variant():
    snpInfo = {}
    delInfo = {}
    insInfo = {}
    varChange = {} # only the high frequency vars
    varCt = 0
    multCt = 0
    snpCt = 0
    snpList = []
    delList = []
    insList = []
    for change in changes:
        varCt += 1
        #data = line.split()
        #site = int(data[0])
        #change = data[1]
        varChange[change] = 1
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

############################################
# get all genetic changes of each isolate
############################################
print("getting all genetic changes (SNPs and indels) of the isolates...")

l = tuple(isoEPIs)
params = {'l': l}
cur.execute('select acc, chg from acc_hap a, hap_chg b where a.hid = b.hid and acc in %(l)s', params)
acc_snp = cur.fetchall()

allSamples = {}
var_count = {}
for line in acc_snp:
    acc = line[0]
    snp = line[1]
    if snp in allSamples:
        allSamples[snp].append(acc)
    else:
        allSamples[snp] = [acc] 
    if acc in var_count:
        var_count[acc] += 1
    else:
        var_count[acc] = 1
    
changes = allSamples.keys()
print("total genetic changes:%s" % len(changes))

#get a variant consequence
get_variant()
    
############################################
# output mutation position and freq
############################################
output_chg = open('var_%s.tsv'% sys.argv[1], 'w')
print("exporting genetic changes with freq >= 0.1%...")

for change in changes:
    iso = allSamples[change]
    freq = float(len(iso))/float(len(isoEPIs))
    if freq >= freqCut:
        site = ''
        if re.search("\d+-\d+", change): # a deletion
            site = change.split("-")[0]
        elif re.search("\d+_[ATCG]+", change): # an insertion
            site = change.split("_")[0]
        else: # a snv
            site = re.sub('[A-Z]', '', change)
        output_chg.write(site + "\t" + change + "\t" + str(len(iso)) + "\t" + str(round(freq,6)) + "\n")
output_chg.close()

print("Done!")
#sys.exit



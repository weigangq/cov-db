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
import numpy
from numpy.random import choice

# Initialize parser
parser = argparse.ArgumentParser(
    description='Sample isolates and filter high-frequency variants')

# Add arguments
parser.add_argument('-t', '--topmost', type = int, nargs = '?', const = 1000,
                    help = 'list countries with at least (default 1000) sequenced genomes')

parser.add_argument('-c', '--country', 
                    help = 'Country name. Quote for names containing blanks, e.g., "South Africa"')

parser.add_argument('-f', '--freq_cut', type = float,
                    help = 'provide minimum variant frequency (default 0.005)', default = 0.005)

parser.add_argument('-l', '--locus', 
                    help = 'select a locus')

parser.add_argument('-m', '--missense', action = 'store_true',
                    help = 'missense SNPs only')

parser.add_argument('-r', '--r4s', 
                    help = 'add rate4site output to var output (require --snp and --locus options)')

parser.add_argument('-s', '--snp', action = 'store_true',
                    help = 'SNPs only')

parser.add_argument('-p', '--per_month', type = int, default = 100,
                    help = 'Sample size per month (default 100)')

args = parser.parse_args()

freqCut = args.freq_cut    
refEPI = 'ISL_406030'
batEPI = 'ISL_402131'
logging.basicConfig(level = logging.DEBUG)

# From https://www.postgresqltutorial.com/postgresql-python/connect/
dbParams = config()
conn = psycopg2.connect(**dbParams)
cur = conn.cursor()

def countByCountry():
    logging.info("Count isolates for topmost countries ...")
    cur.execute("select country, count(*) as count_total from vhuman_anno where acc != %s and acc != %s group by country having count(*) > %s order by count_total desc", [refEPI, batEPI, args.topmost])
    ctryData = cur.fetchall()            
    for ctry in ctryData:
        print(ctry[0] + "\t" + str(ctry[1]))


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

        if args.locus and dataSNP[6] != args.locus:
            continue

        if args.missense and conseq != 'missense':
            continue

        aaCode = dataSNP[4] + str(int((dataSNP[7]-1)/3)+1) if isCoding else 'NA'
            
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
            'aaID': aaCode + dataSNP[5] if conseq == 'missense' else 'NA',
            'r4s': str(rate4site[aaCode]) if aaCode in rate4site else 'NA',
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
            
            if args.locus and dataDEL[1] != args.locus:
                continue

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
                'aaID': 'NA',
                'r4s': 'NA'
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
            refBase = dataINS[4][dataINS[3]-1] if dataINS[3] else dataINS[4]

            if args.locus and dataINS[1] != args.locus:
                continue

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
                'aaID': 'NA',
                'r4s': 'NA'
            }

    sitesSNPs = list(snpInfo.keys())
    sitesDELs = list(delInfo.keys())
    sitesINSs = list(insInfo.keys())

    # sort by site (for printing)
    snpDict = dict(sorted(snpInfo.items(), key=lambda item: item[1]['site']))
    delDict = dict(sorted(delInfo.items(), key=lambda item: item[1]['site']))
    insDict = dict(sorted(insInfo.items(), key=lambda item: item[1]['site']))

    varOut = open('var_%s.tsv'% ctryName, 'w')
    logging.info("Excluding low freq vars: n = %s", lowFreq)
    logging.info("exporting genetic changes with freq >= %s ...", freqCut)
    for site in snpDict.keys(): # str key
        varOut.write(snpInfo[site]['vartype'] + "\t"+ str(snpInfo[site]['site']) + "\t" + snpInfo[site]['varID'] + "\t" + snpInfo[site]['refNT'] + "\t" + snpInfo[site]['locus'] + "\t" + snpInfo[site]['codonRef'] + "\t" + snpInfo[site]['codonPos'] + "\t" + snpInfo[site]['locPos'] + "\t" + snpInfo[site]['altNT'] + "\t" + snpInfo[site]['conseq'] + "\t"  + snpInfo[site]['freq'] + "\t" + snpInfo[site]['count'] + "\t" + snpInfo[site]['aaID'] + "\t" + snpInfo[site]['r4s'] + "\t" + ctryName + "\n")

    if args.missense or args.snp:
        return

    for site in delDict.keys(): # str key!!
        varOut.write(delInfo[site]['vartype'] + "\t"+ str(delInfo[site]['site']) + "\t" + delInfo[site]['varID'] + "\t" + delInfo[site]['refNT'] + "\t" + delInfo[site]['locus'] + "\t" + delInfo[site]['codonRef'] + "\t" + delInfo[site]['codonPos'] + "\t" + delInfo[site]['locPos'] + "\t" + delInfo[site]['altNT'] + "\t" + delInfo[site]['conseq'] + "\t"  + delInfo[site]['freq'] + "\t" + delInfo[site]['count'] + "\t" + delInfo[site]['aaID'] + "\t" + delInfo[site]['r4s'] + ctryName + "\n")

    for site in insDict.keys(): # str key!!
        varOut.write(insInfo[site]['vartype'] + "\t"+ str(insInfo[site]['site']) + "\t" + insInfo[site]['varID'] + "\t" + insInfo[site]['refNT'] + "\t" + insInfo[site]['locus'] + "\t" + insInfo[site]['codonRef'] + "\t" + insInfo[site]['codonPos'] + "\t" + insInfo[site]['locPos'] + "\t" + insInfo[site]['altNT'] + "\t" + insInfo[site]['conseq'] + "\t"  + insInfo[site]['freq'] + "\t" + insInfo[site]['count'] + "\t" + insInfo[site]['aaID'] +  "\t" + insInfo[site]['r4s'] + ctryName + "\n")
    varOut.close()

def main():
############################################
# get all samples of a country (except the ref and bat)
############################################
    logging.info("getting all samples of %s from the database...", args.country)
    cur.execute("select acc, col_date, country, state, area_id from vhuman_anno where acc != %s and acc != %s and country = %s order by col_date", [refEPI, batEPI, args.country])
    accData = cur.fetchall()            
    isoPerMonth = {}
    isoData = {}
    for iso in accData:
        iso = list(iso) # make it list, since tuple is not editable
        if iso[3] is None:
            iso[3] = 'NA'

        iso[1] = re.sub(r"^(\d{4}-\d{2})$", r"\1-15", iso[1])
            #    colDate = datetime.datetime.strptime(iso[1], "%Y-%m-%d")

        yearMonth = re.sub(r"-\d{2}$", r"", iso[1])
        if re.match(r"2019", yearMonth):
            yearMonth = '2020-01'

        if yearMonth in isoPerMonth:
            isoPerMonth[yearMonth].append(iso[0])
        else:
            isoPerMonth[yearMonth] = [iso[0]]

        isoData[iso[0]] = { 'col_date': iso[1],
                            'country': iso[2],
                            'state': iso[3],
                            'area': "Area_" + str(iso[4]),
                            'var_ct': 0 # initialize
                        }

###########################
# sample isolates evenly among months
###############################
    isoEPIs = []
    for yearMonth in isoPerMonth:
        isoCt = isoPerMonth[yearMonth]
        if len(isoCt) > args.per_month:
            isoChoose = numpy.random.choice(isoCt, size = args.per_month, replace = False)
            isoPerMonth[yearMonth] = isoChoose
        logging.info("sample for month %s: %s", yearMonth, str(len(isoPerMonth[yearMonth])))

    for yearMonth in isoPerMonth:
        isoCt = isoPerMonth[yearMonth]
        for iso in isoCt:
            isoEPIs.append(iso)

############################################
# get all genetic changes of each isolate
############################################

    logging.info("getting all genetic changes: SNPs and indels...")
    l = tuple(isoEPIs)
    params = {'a': l}
    cur.execute('select acc, chg, a.hid from acc_hap a, hap_chg b where a.hid_full = b.hid and acc in %(a)s', params)
    acc_lines = cur.fetchall() 

    hidMajor = {}
    cur.execute("select * from hap_lineage")
    lineage_lines = cur.fetchall() 
    for line in lineage_lines:
        hidMajor[line[0]] = line[1]

    allSamples = {}
    var_count = {}
    lineage = {}
    for line in acc_lines:
        acc = line[0]
        chg = line[1]
        lineage[acc] = hidMajor[line[2]]
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
    global varFreq # needs to be accessed in get_variant
    global varAccCt
    global ctryName
    ctryName = args.country.replace(' ', '_')
    varFreq = {}
    varAccCt = {}
    for change in changes:
        iso = allSamples[change]
        freq = float(len(iso))/float(len(isoEPIs))
        varFreq[change] = freq
        varAccCt[change] = len(iso)

    get_variant(changes)

# write iso file
    acc_output = open('iso_%s.tsv'% ctryName, 'w')
    for iso in isoEPIs:
        if iso in var_count: # hid = 1 is not collected
            isoData[iso]['var_ct'] = var_count[iso]
        acc_output.write(iso + "\t" + 
                         isoData[iso]['col_date'] + "\t" + 
                         isoData[iso]['country'] + "\t" + 
                         isoData[iso]['state'] + "\t" + 
                         isoData[iso]['area'] + "\t" + 
                         str(isoData[iso]['var_ct']) + "\t" +
                         lineage[iso] + "\n"
        )
    acc_output.close()

    logging.info("total isolates: %s", len(isoEPIs))
    logging.info("Done!")
#    sys.exit()

rate4site = {}
if args.r4s is not None:
    if args.snp is None or args.locus is None or args.country is None:
        logging.info("--r4s works only with --locus, --snp, and --country")
        sys.exit()
    else:
        logging.info("reading rate4site file....")
        with open(args.r4s) as r4s:
            lines = r4s.readlines()
            for line in lines:
                lineMatch = re.match(r"(\d+)\s+(\S)\s+(-*\d\.\d+)\s+", line)
                if lineMatch:
                    aaSite = str(lineMatch.group(1))
                    refAA = lineMatch.group(2)
                    rate = lineMatch.group(3)
                    rate4site[refAA + aaSite] = rate
    
if args.topmost is not None:
    countByCountry()
else:
    if args.country is None:
        logging.info("Provide a country name: --country <name>")
    else:
        main()

sys.exit()



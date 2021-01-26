#!/usr/bin/env python
# coding: utf-8

import re
import psycopg2
import argparse
#import datetime
import sys

if len(sys.argv) != 2:
    print("Usage: snp-freq.py country_name", file == sys.stderr)
    sys.exit()

freqCut = 0.001    
refEPI = 'EPI_ISL_406030'
batEPI = 'EPI_ISL_402131'

conn = psycopg2.connect("host=borreliabase.org dbname=bb3-dev user=lab password=homology")
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
# get all genetic changes of each isolate
############################################
print("getting all genetic changes (SNPs and indels) of the isolates...")

l = tuple(isoEPIs)
params = {'l': l}
cur.execute('select acc, chg from acc_hap a, hap_chg b where a.hid = b.hid and acc in %(l)s', params)
acc_snp = cur.fetchall()

allSamples = {}
for line in acc_snp:
    acc = line[0]
    snp = line[1]
    if snp in allSamples:
        allSamples[snp].append(acc)
    else:
        allSamples[snp] = [acc] 

snps = allSamples.keys()
#snps.sort()
print("total genetic changes:%s" % len(snps))
    
############################################
# output mutation position and freq
############################################
snp_output = open('var_%s.tsv'% sys.argv[1], 'w')
print("exporting genetic changes with freq >= 0.1%...")

for snp in snps:
    iso = allSamples[snp]
    freq = float(len(iso))/float(len(isoEPIs))
    if freq >= freqCut:
        site = ''
        if re.search("\d+-\d+", snp): # a deletion
            site = snp.split("-")[0]
        elif re.search("\d+_[ATCG]+", snp): # an insertion
            site = snp.split("_")[0]
        else: # a snv
            site = re.sub('[A-Z]', '', snp)
        snp_output.write(site + "\t" + snp + "\t" + str(len(iso)) + "\t" + str(round(freq,6)) + "\n")
snp_output.close()

print("Done!")
sys.exit



# Workflow for COV mutation analysis
* Updated on: Feb 7, 2021 *

## Haplotype pipeline
- Command interface:
```
 pick-iso-var.py --help
 export-vcf.py --help
```
- Show contries with 1000 or more genomes: 
`pick-iso-var.py --topmost`
- Export iso & var for a country:
```
pick-iso-var.py --country India # default sample size & var freq: at most 100 per month; at least 0.5% var frequency 
pick-iso-var.py --country India --freq_cut 0.01 # change minimum freq
pick-iso-var.py --country India --per_month 200 # change most number per month
```
- Filter on a locus:
`pick-iso-var.py --country India --locus S`
-Filter on missense:
`pick-iso-var.py --country India --missense`
-Combine options:
`pick-iso-var.py --country Inida --missense --locus S`


## Step 1. GISAID download, parsing into database, update browser
- Personnel: Edgar & Lia
- Tools: Python, Perl, SQL & JavaScripts
- Outputs: 
  - Relational database populated
  - Genome browser updated

- Main tables & schema:
  * Annotation: human_anno; vgeo
  * Variants: cv_snp; cv_del; cv_ins
  * Genotype: acc_hap; hap_chg

## Step 2. Export isolate EPIs & high-frequency variants (cutoff: derived allele frequency >= 0.1%)
- Personnel: Lily
- Tools: 
  - `export-maf.py <country name>`

- Main queries:
 ```
select a.acc, a.col_date, b.area, b.country from human_anno a, vgeo b where a.need is true and a.geo_id=b.geo_id and a.acc != %s and a.acc != %s and b.country = %s order by col_date;`
select acc, chg from acc_hap a, hap_chg b where a.hid = b.hid and acc in %(l)s;
```
- Outputs:
  - iso.tsv (isolates, geo, collection dates)
  - var.tsv (variants with consequence calls and frequency/counts)

- TO DO:	
  - [ ] implement command-line options (for continent & country)	
  - [ ] standardize dates

## Step 3. Export vcf.gz file
- Personnel: Weigang
- Script: 
  	  `export-vcf.py --iso <iso.tsv> --var <snp.tsv>`
- Main queries:
```
select * from cv_snp where concat(alt, site) = %s;
select * from cv_del where del = %s;     
select * from cv_ins where ins = %s;
select acc, chg from acc_hap a, hap_chg b where a.hid = b.hid and acc in %(l)s;
```
- Outputs:
  - country.vcf (gzip to vcf.gz)

## Step 4. Mutation analysis
- Personnel: Saymon, Lily & Desiree
- Scripts & tools: vcftools, bcftools, ldhat, hapview, TCS, RStudio
- Outputs:
  - time series of maf (> 1% in a country)	
  - syn & nonsyn diversity (ps & ps) per gene	
  - recombination analysis: linkage maps, ldhat, hapview	
  - haplotype networks: TCS & tcsBU


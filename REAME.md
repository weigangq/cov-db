Workflow for COV mutation analysis
Updated on: January 24, 2021

########################
# Step 1. GISAID download, parsing into database, update browser
######################################
Personnel: Edgar & Lia
Python, Perl, SQL & JavaScripts
Output: 
	Relational database populated
	Genome browser updated
Main tables & schema:
     (1) Annotation: human_anno; vgeo
     (2) Variants: cv_snp; cv_del; cv_ins
     (3) Genotype: acc_hap; hap_chg

##################################
# Step 2. Export isolate EPIs & high-frequency variants (cutoff: derived allele frequency >= 0.1%)
####################################
Personnel: Lily
MScript: export-maf.py <country name>
Main queries:
     select a.acc, a.col_date, b.area, b.country from human_anno a, vgeo b where a.need is true and a.geo_id=b.geo_id and a.acc != %s and a.acc != %s and b.country = %s order by col_date;
     select acc, chg from acc_hap a, hap_chg b where a.hid = b.hid and acc in %(l)s;
Outputs:
	iso.tsv
	var.tsv

##############################
# Step 3. Export vcf.gz file
###############################
Personnel: Weigang
Script: export-vcf.py --iso <iso.tsv> --var <snp.tsv> --vcf country.vcf
Main queries:
     select * from cv_snp where concat(alt, site) = %s;
     select * from cv_del where del = %s;
     select * from cv_ins where ins = %s;
     select acc, chg from acc_hap a, hap_chg b where a.hid = b.hid and acc in %(l)s;
Outputs:
	var.tsv (contains consequence calls)
	country.vcf (gzip to vcf.gz)

##############################
# Step 4. Mutation analysis
##################################
Personnel: Saymon, Lily & Desiree
Scripts & tools: vcftools, bcftools, ldhat, hapview, TCS, RStudio
Outputs:
	time series of maf (> 1% in a country)
	syn & nonsyn diversity (ps & ps) per gene
	recombination analysis: linkage maps, ldhat, hapview
	haplotype networks: TCS	& tcsBU


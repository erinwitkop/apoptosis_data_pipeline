#!/usr/bin/env python
#extract gene ID and gene attributes from subsetted GFF3 file based on significantly differentially expressed gene accessions

import sys

for line in sys.stdin:
    attr = dict(item.strip().split(' ') for item in line.split('\t')[8].strip('\n').split('') if item)
    sys.stdout.write("%s\n" % (attr['gene_id'].strip('\"') + '\t' + attr['gene_name'].strip('\"')))
    
#Code reference: Reynolds, Alex. 2015. Forum Answer to "Question: extract only geneID and gene symbol from GTF file".https://www.biostars.org/p/140471/.
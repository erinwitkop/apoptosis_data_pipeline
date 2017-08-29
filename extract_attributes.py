#!/usr/bin/env python3

#To be used as a GFF3 file parser in conjunction with fetchEnsembl_ID.sh
#Extract gene_ID, transcript_ID, UniProtKB Accession from subsetted GFF3 file based on significantly differentially expressed gene accessions

from collections import namedtuple
import gzip
import urllib.request, urllib.parse, urllib.error

#Initialized GeneInfo named tuple. Note: namedtuple is immutable
gffInfoFields = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes", "gene_id", "transcript_id", "description"]
GFFRecord = namedtuple("GFFRecord", gffInfoFields)

def parseGFFAttributes(attributeString):
    """Parse the GFF3 attribute column and return a dict"""#
    if attributeString == ".": return {}
    ret = {}
    if attributeString == "-": return {} #added because my file has this character in a column
    ret = {}
    for attribute in attributeString.split(";"):
        key, value = attribute.split("=")
        ret[urllib.parse.unquote(key)] = urllib.parse.unquote(value)
    return ret

def parseGFF3(filename):
    """
    A minimalistic GFF3 format parser.
    Yields objects that contain info about a single GFF3 feature.
    
    Supports transparent gzip decompression.
    """
    #Parse with transparent decompression
    openFunc = gzip.open if filename.endswith(".gz") else open
    with openFunc(filename) as infile:
        for line in infile:
            if line.startswith("#"): continue
            parts = line.strip().split("\t")
            #If this fails, the file format is not standard-compatible
            #assert len(parts) == len(gffInfoFields)
            #Normalize data
            normalizedInfo = {
                "seqid": None if parts[0] == "." else urllib.parse.unquote(parts[0]),
                "source": None if parts[1] == "." else urllib.parse.unquote(parts[1]),
                "type": None if parts[2] == "." else urllib.parse.unquote(parts[2]),
                "start": None if parts[3] == "." else int(parts[3]),
                "end": None if parts[4] == "." else int(parts[4]),
                "score": None if parts[5] == "." else float(parts[5]),
                "strand": None if parts[6] == "." else urllib.parse.unquote(parts[6]),
                "phase": None if parts[7] == "." else urllib.parse.unquote(parts[7]),
                "attributes": parseGFFAttributes(parts[8])
            }
          
            #Alternatively, you can emit the dictionary here, if you need mutability:
            #    yield normalizedInfo
            #yield GFFRecord(**normalizedInfo)
          
           #add dictionaries for specific parts of the attributes field I want, can change and add, make sure to add all titles to gffInfoFields
            #for k in normalizedInfo["attributes"]:
             #   if k == 'gene_id':
              #          normalizedInfo['gene_id'] = k['gene_id']
               # if k == 'transcript_id':
                #        normalizedInfo['transcript_id'] = k['transcript_id']
                #if k == 'description':
                 #       normalizedInfo['description'] = k['description']

          
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("file", help="The GFF3 input file (.gz allowed)")
    parser.add_argument("--print-records", action="store_true", help="Print all GeneInfo objects, not only")
    parser.add_argument("--filter-type", help="Ignore records not having the given type")
    args = parser.parse_args()
    outfile = open('resoshv1.csv', 'w') #create a python file named for the data you're running, in write mode
    file = open("resoshv1_05_dfSig_matching_gff_lines.txt", 'a')
    #Execute the parser
    recordCount = 0 
    for record in parseGFF3(args.file):
        #Apply filter, if any
        if args.filter_type and record.type != args.filter_type: continue
        #Print record if specified by the user
        if args.print_records: 
        	print(record)
        	outfile.write(record) #write all the records to the file name specified by the outfile variable
        #Access attributes like this: my_strand = record.strand
        recordCount += 1
    print("Total records: %d" % recordCount)
    outfile.close()

#Adapted from code by Uli Kohler, Apache License v2.0, Version 1.1. Accessed via BioStars Forum Post by Reynolds, Alex. 2015. Forum Answer to "Question: extract only geneID and gene symbol from GTF file".https://www.biostars.org/p/140471/.
#https://techoverflow.net/2013/11/30/a-simple-gff3-parser-in-python/ 



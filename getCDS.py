#!/usr/bin/python
import sys
import re
import shutil
from Bio import SeqIO, SeqFeature
from Bio.SeqRecord import SeqRecord
import argparse

"""This script reads a genbank or embl and extract CDS features that are not tagged as pseudogenes."""


def search_CDS(embl,form,fout):

    fhandoutcds = open(fout + ".fna","w")
    fhandoutfaa = open(fout + ".faa","w")
    genome = open(fout +".fasta","w")
    for seq_record in SeqIO.parse(open(embl),iform):
        for feature in seq_record.features:
            if feature.type == "CDS" and "pseudo"not in feature.qualifiers.keys():
                header = []
                locus = "".join(feature.qualifiers["locus_tag"])
                protein = feature.location.extract(seq_record).seq.translate(table=11,to_stop=True)
                cds = feature.location.extract(seq_record).seq
	        if "gene" in feature.qualifiers.keys():
                    header.append("".join(feature.qualifiers["gene"]))
                    header.append("".join(feature.qualifiers["product"]))
		    fhandoutcds.write(">{} {}\n{}\n".format(locus," ".join(header),cds))
                    fhandoutfaa.write(">{} {}\n{}\n".format(locus," ".join(header),protein))
		else:
                    header.append("|".join(feature.qualifiers["product"]))                   
		    fhandoutcds.write(">{} {}\n{}\n".format(locus," ".join(header),cds))
                    fhandoutfaa.write(">{} {}\n{}\n".format(locus," ".join(header),protein))
        genome.write(">{}\n{}\n".format(seq_record.id,seq_record.seq))
    genome.close()	
    fhandoutcds.close()
    fhandoutfaa.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="This script reads a genbank or embl and extract CDS features.")
    parser.add_argument("-a","--annotation",action='store',required = True, help = "Annotation file.")
    parser.add_argument("-if","--informat",action='store',required = True, help = "Input format: genbank or embl.")
    parser.add_argument("-o","--outfile",action='store',required = True, help = "Output file basename")
    args = parser.parse_args()
    annot = args.annotation
    iform = args.informat 
    fout = args.outfile
    embl = annot
    search_CDS(embl,iform,fout)
    print """\nWork is done, take a beer!\n

         _~~-_~~--_
        (          )
         |~-~-~-~~|
       	=|) :( ) (|
       | |): ( ):(|
       	=|)::( ).(|
         |):;( ) (|
          ~~----~~
\n
"""

#!/usr/bin/python
from Bio import SeqIO
from sys import argv
from Bio.Alphabet import IUPAC

"""Extracts 1+2 and 3 codon positions from a multifasta or fasta alingment"""

fhand = open(argv[1],"r") #Fasta alingment
fout12 = open(argv[1] + "_12.fasta","w")
fout3 = open(argv[1] + "_3.fasta","w")

for record in SeqIO.parse(fhand,"fasta"):
	fout12.write(">" + record.id + "\n")
	fout3.write(">" + record.id + "\n")
	firstsecondc = ''
	thirdc = ''
	first = str(record.seq[0::3]) #Slice the sequence
	second = str(record.seq[1::3])
	third = str(record.seq[2::3])
	for i in range(len(third)):
		firstsecondc += first[i] #Join 1+2 position
		firstsecondc += second[i]#Join 3 position
		thirdc += third[i]
	fout12.write(firstsecondc + "\n")
	fout3.write(thirdc + "\n")


	

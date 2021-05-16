#!/usr/bin/python

from Bio import SeqIO
import argparse

"""This script reads a genbank file and generates the karyotipe, several feature files, calculates GC/AT content and GC Skew to use with Circos."""

def gbk2genelist(gbk_file, out_file,out_file2,out_karyo,tRNA,rRNA,miscRNA,colors):
	colors2 = list(colors)
	with open(out_file, 'w') as f1, open(out_file2,'w') as f2, open(out_karyo,'w') as f3, open(tRNA,"w") as f4, open(rRNA,"w") as f5, open(miscRNA,"w") as f6:
		for record in SeqIO.parse(gbk_file, "genbank"):
			name = record.id
			chstart = 0 
			chend = len(record.seq)
			if len(colors2) == 0:
				colors2 = list(colors)
			color = colors2.pop(0)
				
			f3.write("chr - {} {} {} {} {}\n".format(name, name, chstart, chend, color))
			for feature in record.features:	   	
				if feature.type == 'gene':#in the original is gene to plot also pseudogenes   
					location = feature.location.parts
					for part in location:
						strand, start, end = part.strand, part.start, part.end
						if strand == -1:
							colorstrand = "fill_color=vdred"
							start, end = end, start
							if 'pseudo' in feature.qualifiers.keys():
								colorstrand = "fill_color=gray"
							if "EC_number" in feature.qualifiers.keys():
								colorstrand = "fill_color=vvdpurple"
							f2.write("{}\t{}\t{}\t{}\n".format(name, start, end, colorstrand).replace(' ', '_'))
						else:
							colorstrand = "fill_color=vdblue"
							if 'pseudo' in feature.qualifiers.keys():
								colorstrand = "fill_color=gray"
							if "EC_number" in feature.qualifiers.keys():
								colorstrand = "fill_color=vvdpurple"
							f1.write("{}\t{}\t{}\t{}\n".format(name, start, end, colorstrand).replace(' ', '_'))
				if feature.type == 'tRNA':   
					location = feature.location.parts
					for part in location:
						strand, start, end = part.strand, part.start, part.end
						if strand == -1:
							start, end = end, start
						colorstrand = "fill_color=black"
						f4.write("{}\t{}\t{}\t{}\n".format(name, start, end, colorstrand).replace(' ', '_'))
				if feature.type == 'rRNA':   
					location = feature.location.parts
					for part in location:
						strand, start, end = part.strand, part.start, part.end
						if strand == -1:
							start, end = end, start
						colorstrand = "fill_color=vvdyellow"
						f5.write("{}\t{}\t{}\t{}\n".format(name, start, end, colorstrand).replace(' ', '_'))
				if feature.type == 'misc_RNA' or  feature.type == "tmRNA" or feature.type == "ncRNA":
					location = feature.location.parts
					for part in location:
						strand, start, end = part.strand, part.start, part.end
						if strand == -1:
							start, end = end, start
						colorstrand = "fill_color=vvdpurple"
						f6.write("{}\t{}\t{}\t{}\n".format(name, start, end, colorstrand).replace(' ', '_'))			
	f1.close()
	f2.close()
	f3.close()
	f4.close()
	f5.close()
	f6.close()

def nt_skew(seq, window=100): 
	"""Calculates GC skew (G-C)/(G+C) for multiple windows along the sequence. Returns a list of ratios (floats), controlled by the length of the sequence and the size of the window. Returns 0 for windows without any G/C by handling zero division errors. Does NOT look at any ambiguous nucleotides. Pairs could be GC or AT Modified from SeqtUils.""" 
	gcvalues = []
	atvalues = []
	gcpercentage = [] 
	for i in range(0, len(seq), window): 
		s = seq[i: i + window] 
		g = s.count('G') + s.count('g') 
		c = s.count('C') + s.count('c') 
		try: 
			gcskew = (g - c) / float(g + c) 
		except ZeroDivisionError: 
			gcskew = 0.0 
		gcvalues.append(gcskew)
		s = seq[i: i + window] 
		a = s.count('A') + s.count('a') 
		t = s.count('T') + s.count('t') 
		try: 
			atskew = (a - t) / float(a + t) 
		except ZeroDivisionError: 
			atskew = 0.0 
		atvalues.append(atskew)
		percentage = (g+c)/float((g+c+a+t))
		gcpercentage.append(percentage)
	return gcvalues,atvalues,gcpercentage 


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="This script reads a genbank file and generates the karyotipe, several feature files, calculates GC/AT content and GC Skew to use with Circos.")
	parser.add_argument("-i","--input",action='store',required = True, help = "Genbank file with annotations.")
	parser.add_argument("-b","--basename",action='store',required = True, help = "Basename for Circos files.")
	args = parser.parse_args()
	gbk_file = args.input
	basename = args.basename
	#Output files
	out_file = basename + "_genes_pos.txt"
	out_file2 = basename + "_genes_neg.txt"
	out_karyo = basename + "_karyotype.txt"
	tRNA = basename + "_tRNAs.txt"
	rRNA = basename + "_rRNAs.txt"
	miscRNA = basename + "_miscRNA.txt"
	#enzymes = "enzymes.txt"
	out_file_posgc = basename + "_gcskew_posgc.txt"
	out_file_neggc = basename + "_gcskew_neggc.txt"
	out_file_ccggc = basename + "_cum_gcskewgc.txt"
	out_file_posat = basename + "_gcskew_posat.txt"
	out_file_negat = basename + "_gcskew_negat.txt"
	out_file_ccgat = basename + "_cum_gcskewat.txt"
	out_file_gcpercentage = basename + "_gcpercentage.txt"
	#Colors to use in the karyotype
	colors = ["paired-12-qual-1","paired-12-qual-2","paired-12-qual-3","paired-12-qual-4","paired-12-qual-5","paired-12-qual-6","paired-12-qual-7","paired-12-qual-8","paired-12-qual-9","paired-12-qual-10","paired-12-qual-11","paired-12-qual-12"]

	#Generate feature files
	gbk2genelist(gbk_file, out_file,out_file2,out_karyo,tRNA,rRNA,miscRNA,colors)

	#Calculates the GC SKew and cumulative GC Skew
	with open(out_file_posgc, 'w') as f1, open(out_file_neggc, 'w') as f2, open(out_file_ccggc, 'w') as f3, open(out_file_posat, 'w') as f4, open(out_file_negat, 'w') as f5, open(out_file_ccgat, 'w') as f6, open(out_file_gcpercentage, 'w') as f7 :
		for record in SeqIO.parse(gbk_file, "genbank"):
			name = record.id
			sequence = record.seq
			lenseq = len(record.seq)
			gcvalues,atvalues,gcpercentage = nt_skew(sequence, 500)
			ccgskew = 0
			catskew = 0
			for i in range(0,lenseq,500):
				start = i
				end = i + 500
				if end > lenseq:
					end = lenseq
				gc = float(gcvalues[i/500])
				at = float(atvalues[i/500])
				#Get GC skew and CGCSkww
				if gc >= 0:
					f1.write("{}\t{}\t{}\t{}\n".format(name, start,end, gc))
				else:
					f2.write("{}\t{}\t{}\t{}\n".format(name, start,end, gc))
				ccgskew += gc
				f3.write("{}\t{}\t{}\t{}\n".format(name, start,end, ccgskew))
				#Get AT skew and CATSkww
				if at >= 0:
					f4.write("{}\t{}\t{}\t{}\n".format(name, start,end, gc))
				else:
					f5.write("{}\t{}\t{}\t{}\n".format(name, start,end, gc))
				catskew += at
				f6.write("{}\t{}\t{}\t{}\n".format(name, start,end, catskew))
				#Percentage
				percentage = float(gcpercentage[i/500])
				f7.write("{}\t{}\t{}\t{}\n".format(name, start,end, percentage))		

	f1.close()
	f2.close()
	f3.close()
	f4.close()
	f5.close()
	f6.close()
	f7.close()


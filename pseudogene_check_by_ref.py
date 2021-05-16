#!/usr/bin/python
import sys,re,argparse,subprocess
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
#===============================================================================
#    Author: Diego Santos-Garcia
#    Contact: Contact the author at diego.santos@mail.huji.ac.il or diego.santos.garcia@protonmail.com
#
#    COPYRIGHT: Copyright (C) 2018  Diego Santos-Garcia.
#
#    LICENCE: This program is free software: you can redistribute it and/or modify it under the terms
#             of the GNU General Public License as published by the Free Software Foundation, either
#             version 3 of the License, or (at your option) any later version.
#             This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
#             without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#             See the GNU General Public License for more details.
#             You should have received a copy of the GNU General Public License along with this program.
#             If not, see <http://www.gnu.org/licenses/>.
#===============================================================================

def prepare_reference(genbanks):
	'''Extract proteins, CDS and genome in fasta format from genbank/embl file. Fasta header is locus_tag, gene (if exist) and product'''
	flag = len(genbanks) > 1#Record if we are working with refernece or query
	for genbank in genbanks:
		if flag:
			fout = "Reference"
			modewrite = "a"
		else:
			fout = genbank
			modewrite = "w"
		fhandoutcds = open(fout + ".fna",modewrite)
		fhandoutfaa = open(fout + ".faa",modewrite)
		genome = open(fout +".fasta",modewrite)
		for seq_record in SeqIO.parse(open(genbank),"genbank"):
			for feature in seq_record.features:
				if feature.type == "CDS" and "pseudo" not in feature.qualifiers.keys() and len(feature.location.extract(seq_record).seq)%3 == 0:
					header = []
					locus = "".join(feature.qualifiers["locus_tag"])
					protein = feature.location.extract(seq_record).seq.translate(table=11,to_stop=True)
					cds = feature.location.extract(seq_record).seq
					if "gene" in feature.qualifiers.keys():
							header.append("".join(feature.qualifiers["gene"]))
							print locus
							header.append("".join(feature.qualifiers["product"]))
							fhandoutcds.write(">{}\t{}\n{}\n".format(locus,"|".join(header),cds))
							fhandoutfaa.write(">{}\t{}\n{}\n".format(locus,"|".join(header),protein))
					else:
							header.append("|".join(feature.qualifiers["product"]))
							fhandoutcds.write(">{}\t{}\n{}\n".format(locus,"".join(header),cds))
							fhandoutfaa.write(">{}\t{}\n{}\n".format(locus,"".join(header),protein))
			genome.write(">{}\n{}\n".format(seq_record.id,seq_record.seq))
	genome.close()
	fhandoutcds.close()
	fhandoutfaa.close()

def run_last(querygenome):
	'''Run LAST aligner.'''
	querygenome = querygenome
	queryprot = querygenome + ".faa"
	querygenomefas = querygenome + ".fasta"
	blastproteome = querygenome + ".faa_last.tab"
	blastgenome = querygenome + ".genome_last.tab"
	try:#Collapse reference proteomes
		call = 'cd-hit -i Reference.faa -o Reference.cdhit.faa -c 0.95 -M 0'
		subprocess.call(call,shell=True)
	except:
		print "Check if CD-HIT is installed or accesible."
	try:
		call = 'lastdb -p RefDB Reference.cdhit.faa'
		subprocess.call(call,shell=True)
	except:
		print "Check if LAST is installed or accesible."

	try:
		call = 'lastal -f BlastTab+ -P 0 RefDB -p MIQS -K 1 {} > {}'.format(queryprot,blastproteome)
		subprocess.call(call,shell=True)
	except:
		print "Check if LAST is installed or accesible."

	try:
		call = 'lastal -f BlastTab+ -P 0 RefDB -F15 {} > {}'.format(querygenomefas,blastgenome)
		subprocess.call(call,shell=True)
	except:
		print "Check if LAST is installed or accesible."
	return blastproteome,blastgenome

def read_last(lastfile,blasttype,evaluelast):
	'''Read last BlastTab+ into a dictionary. What about paralogs'''
	lastdict= {}
	lastfile = open(lastfile,"r")
	for line in lastfile:
		if line.startswith("#"):
			pass
		else:
			line = line.strip().split("\t")
			if blasttype == "protein":
				query,subject,alingment,qstart,qend,sstart,send,qlength,slength,evalue,score,coverage = line[0],line[1],line[3],line[6],line[7],line[8],line[9],line[12],line[13],line[10],line[11],float(line[3])/float(line[13])
			if 	blasttype == "genome":#Just change query by subject. Multiple hits will overwrite the dictionary
				if int(line[6]) > int(line[7]):
					query = "{}_{}_{}".format(line[0],line[7],line[6])#Make ids from genomic regions
				else:
					query = "{}_{}_{}".format(line[0],line[6],line[7])#Make ids from genomic regions
				subject,alingment,qstart,qend,sstart,send,qlength,slength,evalue,score,coverage = line[1],line[3],line[6],line[7],line[8],line[9],line[12],line[13],line[10],line[11],float(line[3])/float(line[13])
			if float(evalue) >=  evaluelast:#Don't record low quality alingmenst
				continue
			elif query in lastdict.keys():#Choose the best alingment (score based) if paralogs are found in teh reference
				if float(score) > float(lastdict[query]['score']):#Replace Values in dictionary
					if 'paralogs' in lastdict[query].keys():
							lastdict[query]['paralogs'].append(lastdict[query]['subject'])
					else:
						lastdict[query]['paralogs'] = []
					lastdict[query]['subject'] =  subject
					lastdict[query]['alingment'] =  alingment
					lastdict[query]['qstart'] =  qstart
					lastdict[query]['qend'] =  qend
					lastdict[query]['sstart'] =  sstart
					lastdict[query]['send'] =  send
					lastdict[query]['qlength'] =  qlength
					lastdict[query]['slength'] =  slength
					lastdict[query]['score'] =  score
					lastdict[query]['coverage'] =  coverage

				else:
					if 'paralogs' in lastdict[query].keys():
						lastdict[query]['paralogs'].append(subject)
					else:
						lastdict[query]['paralogs'] = []
			else:
				lastdict[query] = {'subject':subject,'alingment':alingment,'qstart':qstart,'qend':qend,'sstart':sstart,'send':send,'qlength':qlength,'slength':slength,'score':score,'coverage':coverage,'paralogs':[]}
	return lastdict

def full_genes(lastfile,coveragethreshold,evaluelast,blasttype="protein"):
	'''Prints a list of genes which coverage to the reference proteins are at least X%'''
	lastdict = read_last(lastfile,blasttype,evaluelast)
	foutcomplete = open("complete_genes.tab","w")
	#Make a list of complete and imcomplete genes to use later their positions in the genome. It can be duplicated genes pseudogenized
	complete = []
	incomplete = []
	foutincomlete = open("icomplete_genes.tab","w")
	foutcomplete.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("Completeness","Query","Query Length","Subject","Subject Length","Coverage","Paralogs"))
	foutincomlete.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("Completeness","Query","Query Length","Subject","Subject Length","Coverage","Paralogs"))
	for query in sorted(lastdict.keys()):
		if float(lastdict[query]['coverage']) >= coveragethreshold:
			if float(lastdict[query]['coverage']) > 1.1:#Reference in 10% smaller than our protein. It can be a pseudogene in the refernece
				foutincomlete.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("Warning",query,lastdict[query]['qlength'],lastdict[query]['subject'],lastdict[query]['slength'],lastdict[query]['coverage'],",".join(lastdict[query]['paralogs'])))
				incomplete.append(query)
			else:
				foutcomplete.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("Complete",query,lastdict[query]['qlength'],lastdict[query]['subject'],lastdict[query]['slength'],lastdict[query]['coverage'],",".join(lastdict[query]['paralogs'])))
				complete.append(query)
		else:
				foutincomlete.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("Incomplete",query,lastdict[query]['qlength'],lastdict[query]['subject'],lastdict[query]['slength'],lastdict[query]['coverage'],",".join(lastdict[query]['paralogs'])))
				incomplete.append(query)
	return complete,incomplete

def genome_positions(blastgenome,querygenome,evaluelast,blasttype="genome"):
	'''Extract for each protein in the reference genome its orthologous position in the genome to curate.'''
	lastdict = read_last(blastgenome,blasttype,evaluelast)
	refgenomepositions = {}
	for key in lastdict.keys():
		if int(lastdict[key]['qstart']) > int(lastdict[key]['qend']):
			strand = -1
			start,end = lastdict[key]['qend'],lastdict[key]['qstart']
		else:
			strand = +1
			start,end = lastdict[key]['qstart'],lastdict[key]['qend']
		refgenomepositions[key] = {'subject':lastdict[key]['subject'],'start':int(start),'end':int(end),'strand':strand}#Reduce the information
#	for i in refgenomepositions.keys():
#		print i,refgenomepositions[i]
	return refgenomepositions

def get_pseudo_desc(refgenomepositions):
	pseudodesc = {}
	referenceprot = "Reference.cdhit.faa"
	fhand = open(referenceprot,"r")
	for line in fhand:
		if line.startswith(">"):
			line = line.replace(">","").strip().split('\t')
			desc = line[1].split("|")
			locus = line[0]
			if len(desc) == 2:
				pseudodesc[locus] = {'gene':desc[0],'product':desc[1]}
			else:
				product = desc[0]
				pseudodesc[locus] = {'gene':"",'product':desc[0]}
	fhand.close()
	for position in refgenomepositions.keys():
		subject = refgenomepositions[position]['subject']
		if subject in pseudodesc.keys():
			refgenomepositions[position].update(pseudodesc[subject])
	return pseudodesc,refgenomepositions

def tidy_pseudos(refgenomepositions,querygenome,complete,pseudodesc):
	'''Check all possible hits and delete those that are overlapping with complete genes. Mix overlapping/near hits with the same product that are putative pseudogenes.'''
	completegenes = {}
	for seq_record in SeqIO.parse(open(querygenome),"genbank"):
		for feature in seq_record.features:#Capture complete genes coordinates to delete pseudos
			if feature.type == "CDS":
				if "".join(feature.qualifiers['locus_tag']) in complete:
					completegenes["".join(feature.qualifiers['locus_tag'])] = {'start':feature.location.start.position,'end':feature.location.end.position,'contig':seq_record.id}
		for gene in completegenes.keys():
			for pseudo in refgenomepositions.keys():
				contigid = re.sub("_[0-9]+_[0-9]+$","",pseudo)
				if max(0,min(completegenes[gene]['end'],refgenomepositions[pseudo]['end']) - max(completegenes[gene]['start'],refgenomepositions[pseudo]['start'])) >= 100 and contigid == completegenes[gene]['contig']:#Eliminar pseudo si tiene un overlap de mas de 100 nt
					del refgenomepositions[pseudo]
				elif completegenes[gene]['start'] <= refgenomepositions[pseudo]['start'] and completegenes[gene]['end'] >= refgenomepositions[pseudo]['end'] and contigid == completegenes[gene]['contig']:
					del refgenomepositions[pseudo]
	lastiter = 1
	currentiter = 0
	while currentiter != lastiter:	#Iterates until no more pseudos are mixed
		lastiter = len(refgenomepositions.keys())
		for key_a in refgenomepositions.keys():#Code to mix hits with same product
			contigid = re.sub("[0-9]+_[0-9]+$","",key_a)#Keep track of the contig to mix only genes in same contig
			for key_b in refgenomepositions.keys():
				if key_a == key_b:
					break
				elif re.match(contigid,key_b):
					overlap = max(0,min(refgenomepositions[key_a]['end'],refgenomepositions[key_b]['end']) - max(refgenomepositions[key_a]['start'],refgenomepositions[key_b]['start']))
					if refgenomepositions[key_a]['subject'] == refgenomepositions[key_b]['subject'] or refgenomepositions[key_a]['product'] == refgenomepositions[key_b]['product']:#check by blast hit, gene or product
							if overlap > 0:#Mix hits with same last hit that are overlapping and keep only one. modifuy to product or gene name?
								if refgenomepositions[key_a]['start'] > refgenomepositions[key_b]['start']:
									refgenomepositions[key_a]['start'] = refgenomepositions[key_b]['start']
								if 	refgenomepositions[key_a]['end'] < refgenomepositions[key_b]['end']:
									refgenomepositions[key_a]['end'] = refgenomepositions[key_b]['end']
								del refgenomepositions[key_b]
								break
							if refgenomepositions[key_a]['end'] < refgenomepositions[key_b]['start'] and abs(refgenomepositions[key_b]['start'] - refgenomepositions[key_a]['end']) <= 500:#Mix adjacent pseudos closet than 500 nt
								refgenomepositions[key_a]['end'] = refgenomepositions[key_b]['end']
								if refgenomepositions[key_a]['start'] > refgenomepositions[key_b]['start']:
									refgenomepositions[key_a]['end'] = refgenomepositions[key_b]['end']
								del refgenomepositions[key_b]
								break
					else:
						if overlap/float(refgenomepositions[key_a]['end'] - refgenomepositions[key_a]['start']) >= 0.3 or overlap/float(refgenomepositions[key_b]['end'] - refgenomepositions[key_b]['start']) >= 0.3:#Mix overlaping pseudos with different product if overlap 30%
							refgenomepositions[key_a]['start'] = min(refgenomepositions[key_a]['start'],refgenomepositions[key_b]['start'])
							refgenomepositions[key_a]['end'] = max(refgenomepositions[key_a]['end'],refgenomepositions[key_b]['end'])
							del refgenomepositions[key_b]
							break
		currentiter = len(refgenomepositions.keys())
	return	refgenomepositions

def compare_features(querygenome,refgenomepositions,pseudodesc):
	'''Iterate over genes/CDS and check if they overlap with putative pseudogenes. Capture gene/product if different from hypothetical porteins or gene name is not present in pseudodesc'''
	delgenes = []
	for seq_record in SeqIO.parse(open(querygenome),"genbank"):
		for feature in seq_record.features:
			if feature.type == 'source' or feature.type == 'assembly_gap':
				continue
			genestart = feature.location._start.position
			geneend = feature.location._end.position
			locus = "".join(feature.qualifiers["locus_tag"])
			if feature.type == 'gene' or feature.type == 'CDS':#Just in case some gene/CDS are orphans: Locus is the smae for CDS/gene. Only one dict is needed
				for pseudo in refgenomepositions.keys():
					contigid = re.sub("_[0-9]+_[0-9]+$","",pseudo)#Keep track of the contig to mix only genes in same contig
					if contigid == seq_record.id:
						refpseudo = refgenomepositions[pseudo]['subject']
						overlap = max(0,min(geneend,refgenomepositions[pseudo]['end']) - max(genestart,refgenomepositions[pseudo]['start']))
						if overlap/float(geneend - genestart) >= 0.4:
							delgenes.append(locus)
							if 'gene' in feature.qualifiers.keys() and pseudodesc[refpseudo]['gene'] == "":#Append gene names an product names if were annotated in the genome
								pseudodesc[refpseudo]['gene'] = re.sub("_[0-9]+","","".join(feature.qualifiers['gene']))#Prokka appends _[0-9] when more than one gene is found, including pseudogenes. I took it out
							if 'product' in feature.qualifiers.keys():
								if pseudodesc[refpseudo]['product'] == "hypothetical protein" and feature.qualifiers['product'] != "hypothetical protein":
									pseudodesc[refpseudo]['product'] = "".join(feature.qualifiers['product'])
	delgenes = set(delgenes)#Remove duplicated locus
	return delgenes,pseudodesc

def add_to_gbk(querygenome,refgenomepositions,pseudodesc):
	counter = 0#make locus for pseudo genes
	records = []
	for seq_record in SeqIO.parse(open(querygenome),"genbank"):
		newseq_record = SeqRecord(seq=seq_record.seq,id=seq_record.id,name=seq_record.name,description=seq_record.description,annotations=seq_record.annotations,features=seq_record.features)
		for query in refgenomepositions.keys():#Add pseudogenes
			if re.search(seq_record.id + "_",query):
				locus = "PDS_" + str(counter)
				pstart = int(refgenomepositions[query]['start'])
				pend = int(refgenomepositions[query]['end'])
				reflocus = refgenomepositions[query]['subject']
				location = FeatureLocation(pstart,pend,refgenomepositions[query]['strand'])
				inference = "similar to DNA sequence:RefSeq:{}".format(reflocus)
				if pseudodesc[reflocus]['gene'] != "":
					gene = pseudodesc[reflocus]['gene']
					product = pseudodesc[reflocus]['product']
					feature = SeqFeature(location=location,type="gene",qualifiers={'pseudo':"",'locus_tag':locus,'product':product,'gene':gene,'inference':inference})
				else:
					product = pseudodesc[reflocus]['product']
					if re.search("[A-Z][a-z]{2,3}[A-Z]$",product):#add regex to catch protein names ad put as gene name
						gene = list(re.search("[A-Z][a-z]{2,3}[A-Z]$",product).group())
						gene[0] = gene[0].lower()
						gene = "".join(gene)
					else:
						gene = locus
					feature = SeqFeature(location=location,type="gene",qualifiers={'pseudo':"",'locus_tag':locus,'product':product,'gene':gene,'inference':inference})
				newseq_record.features.append(feature)
				counter += 1
		records.append(newseq_record)
	return records

def delete_genes(records,delgenes,fout,incomplete):
	finalrecord = []
	genescheck = open(fout + "_genes_manual_check.txt","w")
	for gene in incomplete:
		if gene not in delgenes:
			genescheck.write(gene + "\n")
	genescheck.close()
	for seq_record in records:
		newseq_record = SeqRecord(seq=seq_record.seq,id=seq_record.id,name=seq_record.name,description=seq_record.description,annotations=seq_record.annotations)
		for feature in seq_record.features:
			if feature.type == "source":
				newseq_record.features.append(feature)
			elif feature.type == "gene" or feature.type == "CDS":#Delete gene/CDS
				if "".join(feature.qualifiers['locus_tag']) in delgenes:
					continue
				else:
					newseq_record.features.append(feature)
			else:
				newseq_record.features.append(feature)
		finalrecord.append(newseq_record)
	SeqIO.write(finalrecord,fout,"genbank")

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="This script takes a reference genome (genbnak format) and a query genome (genbank fromat) and uses LAST to annotate putative pseudogenes assuming a desired coverage cutoff (overlap between reference and query proteins).")
	parser.add_argument("-i","--input",action='store',required = True, help = "Query genome (genbank)")
	parser.add_argument("-r","--reference",action='store',required = True, help = "List of reference genomes (genbank), comma separated")
	parser.add_argument("-o","--output",action='store',required = True, help = "Query genome genbak file with pseudogenes")
	parser.add_argument("-t","--threshold",action='store',required = True, help = "Proein overlap (query against reference) to consider a CDS as a pseudogene")
	parser.add_argument("-e","--evalue",action='store',required = True, help = "Evalue threshold for LAST (float)")
	args = parser.parse_args()
	querygenome = [args.input]
	refgenome = args.reference.split(",")
	coveragethreshold = float(args.threshold)
	fout = args.output
	evaluelast = float(args.evalue)
	print "\nFormatting files.\n"
	prepare_reference(querygenome)
	prepare_reference(refgenome)
	querygenome = "".join(querygenome)
	print "\nBlasting.\n"
	blastproteome,blastgenome = run_last(querygenome)
	print "\nGetting list of complete/incomplete genes.\n"
	complete,incomplete = full_genes(blastproteome,coveragethreshold,evaluelast)
	print "\nGetting reference proteome genome positions.\n"
	refgenomepositions = genome_positions(blastgenome,querygenome,evaluelast)
	print "\nGetting pseudognenes descriptions.\n"
	pseudodesc,refgenomepositions = get_pseudo_desc(refgenomepositions)
	print "\nTidying up pseudogene list.\n"
	refgenomepositions = tidy_pseudos(refgenomepositions,querygenome,complete,pseudodesc)
	print "\nComparing reference and query genes.\n"
	delgenes,pseudodesc =compare_features(querygenome,refgenomepositions,pseudodesc)
	print "\nAdding pseudogenes.\n"
	records = add_to_gbk(querygenome,refgenomepositions,pseudodesc)
	print "\nWriting output..\n"
	delete_genes(records,delgenes,fout,incomplete)

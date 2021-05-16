#!/usr/bin/python
# -*- coding: utf-8 -*

from Bio.Phylo.PAML import codeml
from Bio.Seq import Seq
from Bio import SeqIO,AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC, Gapped
from scipy.stats import chisqprob
import argparse,re,sys,os,subprocess

#===============================================================================
#   Author: Diego Santos-Garcia
#
#   Date: 05-12-2016
#   Version: 0.3
#
#This is the pipeline used on the paper for codeML analysis. The purpouse of it
#is to allow other researchers to check and repeat the analysis.
#
#===============================================================================


#Read orthoMCL groups
def read_groups(group):
	'''Read once the group file and store it in a dictionary
	GROUP_2689: PaAn|PANA5342_RS39575 SoHS|SANT_RS07665 SoMe|SOME_04470'''
	fhand = open(group,'r')
	groups = {}
	for line in fhand:
		line = line.strip()
		header = line.split(":")[0]
		species = line.split(":")[1].split(" ")[1:]
		groups[header] = species
	return groups

#Iterate the dictionary and get species names
def get_organism(groups):
	organism = []
	#Get organism names
	for key in groups.keys():
		for gene in groups[key]:
			organism.append(re.sub("\|.+","",gene))#Extract organism prefix
	organism = set(organism)
	return organism

#Iterate the dictionary, delete non-single copy orthogroups and select core genome
def get_core(groups,organism):
	for key in groups.keys():
		genes = {}
		#Delete non-single copy groups
		for specie in organism:
			for item in groups[key]:
				if specie in item:
					if specie in genes.keys():
						genes[specie].append(item)
					else:
						genes[specie] = []
						genes[specie].append(item)
		for item in genes.keys():
			if len(genes[item]) > 1:
				if key in groups.keys():
					groups.pop(key)
		#Keep only the core
		if len(genes) < len(organism):
			if key in groups.keys():
				groups.pop(key)
	return groups

#Select genes to run codeML
def select_genes(groups,organism,infile):
	results = {}
	sort_keys = sorted(groups.keys())
	for key in sort_keys:
		if key not in results.keys():
			results[key] = {}
		best_model,order = select_fasta(key,groups[key],infile)
		if best_model == None:
			print "Alingment contains a protein that is less than 60% of the largest protein.\n"
		else:
			results[key].update(best_model['NSsites'][0]['parameters']['branches'])
			for branch in results[key].keys():
				results[key][branch].update({'model':best_model['model']})
	return results

#Search for the genes in a concatenated fasta file of all the species, select them and tranlate them to amino acids
def select_fasta(key,genes,infile):
	"""Search for the genes in a concatenated fasta file of all the species, writes the nucleotides and the protein to DoEn.rary files to continue the pipeline."""
	order = []
	geneslength = []
	handle = open(infile,"r")
	fout1 = open("genes.fna","wa")
	print "Extracting genes " + " ".join(genes) + " of group " + key + ".\n"  #Print for checking point
	for record in SeqIO.parse(handle,"fasta"):
		for gene in genes:
			locus = gene.split("|")
			specie = locus[0]
			locus = locus[1]
			if locus in record.description:
				record.id = gene
				record.description = ""
				record.seq = del_stop(record.seq)
				SeqIO.write(record,fout1,"fasta")
				order.append(record.id)
				geneslength.append(len(record.seq))
	if (float(min(geneslength))/float(max(geneslength))) < 0.6:#Stops calculation if the smallest gene is less than 60%  of the biggest gene
		return None,None
	else:
		fout1.close()
		best_model = codon_alingment(key,order)
		handle.close()
		cmd = 'rm genes.*'
		subprocess.check_call(cmd,shell=True)
		return best_model,order


#Make codon alingment with  mafft and pal2nal and prepare unrooted tree for codeml
def codon_alingment(key,order):
	print "Prepraring codon alingment for codeML.\n"
	order = run_macse()
	cmd = "sed -i 's/!/-/g' genes.codon.fasta"  #Delete unclear alingments from macse
	subprocess.call(cmd,shell=True)
	cmd = 'Gblocks genes.codon.fasta -t=c -b5=n'#add gblocks
	subprocess.call(cmd,shell=True)
	print "Codon alingment finished.\n"
	cmd = 'cp genes.codon.fasta macse_alingments/' + key + '.codon.fasta'
	subprocess.check_call(cmd,shell=True)
	cmd = 'cp genes.codon.fasta-gb macse_alingments/' + key + '.codon.fasta-gb'
	subprocess.check_call(cmd,shell=True)
	#Prepare tree
	tree = ["","",""]
	for item in order:
		if "MyTu" in item:
			tree[0] = item
		elif "MyHa" in item:
			tree[1] = item
		elif "MyBr" in item:
			tree[2] = item
		else:
			tree.append(item)
	temptree = "(" + ",".join(tree[0::]) + ");"
	tree = temptree
	fhand = open("genes.newick","w")
	fhand.write(tree)
	fhand.close()
	best_model = run_codeml(order,key,tree)
	return best_model

#MACSE alignment
def run_macse():
	order = []
	genetic_code = '11'
	#No pseudo genes included
	macse_cmd = ['java', '-jar', '/home/diego/Software/macse_v2.03.jar', '-prog', 'alignSequences', '-gc_def', genetic_code, '-out_AA', 'genes.faa', '-out_NT', 'genes.codon.fasta', '-seq', 'genes.fna']
	macse_cmd = " ".join(macse_cmd)
	macse_results = subprocess.check_output(macse_cmd,shell=True)
	with open('genes.codon.fasta','r') as alingment:
		for line in alingment:
			if line.startswith('>'):
				line = line.strip().replace('>',"")
				order.append(line)
	return order


#Delete stop codons for codeml
def del_stop(sequence):
	codon_stop_array = ["TAG", "TGA", "TAA", "UGA", "UAA", "UAG"]
	tempseq = ""
	for index in range(0, len(sequence), 3):
		codon = str(sequence[index:index+3]).upper()
		if codon not in codon_stop_array:
			tempseq += codon
	return Seq(tempseq)

#Run codeML with BioPython
def run_codeml(order,key,tree):
	print "Running codeML for " + str(key) + ".\n"
	cml = codeml.Codeml(alignment = "genes.codon.fasta-gb", tree = "genes.newick", out_file = "genes.codeml")
	define_options(cml)

	#Run model m0 (one ω ratio for all the branches and 5 free parameters)
	cml.set_options(model = 0)
	print "Running codeML m0 for " + str(key) + ".\n"
	resultsm0 = best_log(cml)

	#Run model m1 - free ω ratios for branches and 7 free parameters
	cml.set_options(model = 1)
	print "Running codeML m1 for " + str(key) + ".\n"
	resultsm1 = best_log(cml)

	#Run model m2 - one ω ratio for background branches and one ω for foreground (intracellular), 6 free parameters
	add_background1(order,tree)
	cml = codeml.Codeml(alignment = "genes.codon.fasta-gb", tree = "genes.newick", out_file = "genes.codeml")
	define_options(cml)
	cml.set_options(model = 2)
	print "Running codeML m2 for " + str(key) + ".\n"
	resultsm2 = best_log(cml)

	#Compare models by LTR test
	best_model = LTR_test(resultsm0,resultsm1,resultsm2)
	#print best_model
	#Added for rooted tree
	branches = best_model['NSsites'][0]['parameters']['branches'].keys()
	for item in branches:#Because the internal nodes 5 and 4 can be changed depending on the input order, we use the key to selec only the leafs and not the internal nodes independently of the numbering
		if "..1" in item:
			best_model['NSsites'][0]['parameters']['branches'][(item.replace("1","") + order[0])] = best_model['NSsites'][0]['parameters']['branches'].pop(item)#Replace default branch names by species|gene convention
		if "..2" in item:
			best_model['NSsites'][0]['parameters']['branches'][item.replace("2","") + order[1]] = best_model['NSsites'][0]['parameters']['branches'].pop(item)
		if "..3" in item:
			best_model['NSsites'][0]['parameters']['branches'][item.replace("3","") + order[2]] = best_model['NSsites'][0]['parameters']['branches'].pop(item)
		if "..4" in item:
			best_model['NSsites'][0]['parameters']['branches'][item.replace("4","") + order[3]] = best_model['NSsites'][0]['parameters']['branches'].pop(item)
	return best_model


def add_background1(order,tree):#Add leprae and lepromatosis as foreground branch
	tree = ["","",""]
	for item in order:
		if "MyTu" in item:
			tree[0] = item
		elif "MyHa" in item:
			tree[1] = item
		elif "MyBr" in item:
			tree[2] = item
		else:
			tree.append(item)
	temptree = "(" + "".join(tree[1:2]) + ",(" + ",".join(tree[2::]) + ") $1));"
	temptree = "(" + str(tree[0]) + "," + temptree
	tree2 = temptree
#	tree2 = tree.replace(")))",") $1 ))")
	fhand = open("genes.newick","w")
	fhand.write(tree2)
	fhand.close()

def best_log(cml):#run codeML three times, store the runs in different dictionaries, compare log-likelihoods and kept best of them
	i = 0
	run =  {}
	best_run = {}
	while i < 3:
		run[i] =  cml.run(verbose=False)
		i += 1
	best_lnL = 0
	best_lnL_key = 0
	for key in run:
		if run[key]['NSsites'][0]['lnL'] < best_lnL:
			best_lnL = run[key]['NSsites'][0]['lnL']
			best_lnL_key = key
	best_run = run[best_lnL_key]
	return best_run

def LTR_test(resultsm0,resultsm1,resultsm2):
	model0 = float(resultsm0['NSsites'][0]['lnL'])#5 free parameters
	model1 = float(resultsm1['NSsites'][0]['lnL'])#7 free parameters
	model2 = float(resultsm2['NSsites'][0]['lnL'])#6 free paremeters
	df = 1 #7-6 free parameters
	#Compare fisrt model1 against model2
	value = 2*(model1 - model2)
	pvalue = chisqprob(value, df)
	if pvalue <= 0.025:#Bonferroni correction (0.05/2) for two tests
		#Model 1 is better, check against model 0
		df = 2
		value = 2*(model1 - model0)
		pvalue = chisqprob(value, df)
		if pvalue <= 0.025:#Bonferroni correction (0.05/2) for two tests
			return resultsm1
		else: # One ω ratio for all the branches can't be discarded or model 0 (H0)
			return resultsm0
	else:
		#Model 2 is better, check against model 0
		df = 1
		value = 2*(model2 - model0)
		pvalue = chisqprob(value, df)
		if pvalue <= 0.025:#Bonferroni correction (0.05/2) for two tests
			return resultsm2
		else: # One ω ratio for all the branches can't be discarded or model 0 (H0)
			return resultsm0


def print_general(results):#Prints the matrix for all genes by groups, winning omega model and all codeml parameters
	fout = open("codeml_gene_all_parameters_matrix.tab","w")
	parameters = ["model","N","S","dS","S*dS","dN","N*dN","omega","t"]
	array = [parameters]
	for key in sorted(results.keys()):
		for species in results[key]:
			item = []
			for parameter in parameters:
				item.append(results[key][species][parameter])
			item.insert(0,species)
			item.insert(0,key)
			array.append(item)
	array[0].insert(0,"Organism|Gene")
	array[0].insert(0,"Group")
	fout.write("\t".join(map(str,array[0])) + "\n")
	for i in range(1,len(array)):
		fout.write("\t".join(map(str,array[i])) + "\n")
	fout.close()


def define_options(cml):##Define the options for codeML
	cml.set_options(verbose = 1)
	cml.set_options(CodonFreq = 2)
	cml.set_options(cleandata = 1)
	cml.set_options(fix_blength = 1)
	cml.set_options(NSsites = [0])#Should be a list if you want more than one NSites analisys [0,1,2]
	cml.set_options(fix_omega = 0) # Ha: fix_omega = 0 ; Ho: fix_omega = 1
	cml.set_options(clock = 0)
	cml.set_options(ncatG = 8)
	cml.set_options(runmode = 0)
	cml.set_options(fix_kappa = 0)
	cml.set_options(fix_alpha = 1)
	cml.set_options(Small_Diff = 5e-6)
	cml.set_options(method = 1)
	cml.set_options(Malpha = 0)
	cml.set_options(aaDist = 0)
	cml.set_options(RateAncestor = 0)
	cml.set_options(aaRatefile = "/usr/share/paml/dat/wag.dat")#Select the desired aa rate: dayhoff.dat, jones.dat, wag.dat, mtmam.dat, or your own aa rate file
	cml.set_options(icode = 0)#Genetic code:  0:universal, 1:mammalian mt., 2:yeast mt., 3:mold mt.,* 4: invertebrate mt., 5: ciliate nuclear,6: echinoderm mt., 7: euplotid mt., 8: alternative yeast nu. 9: ascidian mt.,  10: blepharisma nu.
	cml.set_options(alpha = 0)
	cml.set_options(seqtype = 1)# 1:codons; 2:AAs; 3:codons-->AAs
	cml.set_options(omega = 1)
	cml.set_options(getSE = 0)
	cml.set_options(noisy = 0)
	cml.set_options(Mgene = 0)
	cml.set_options(kappa = 2)
	cml.set_options(ndata = 1)
	cml.ctl_file = "codeml.ctl"#Print las control file used for checking purpouses
	cml.write_ctl_file()

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Run codeML and extract desired information from it.")
	parser.add_argument("-g","--group",action='store',required = True, help = "File with ortrogroups in tabular format: \"GROUP_2689: PaAn|PANA5342_RS39575 SoHS|SANT_RS07665\"")
	parser.add_argument("-i","--infile",action='store',required = True, help = "Concatenated files of all genes from all the species under analysis.")
	args = parser.parse_args()
	group = args.group
	infile = args.infile
	#Read orthogroups from tabular file
	if not os.path.isdir("macse_alingments"):
		cmd = 'mkdir macse_alingments'
		subprocess.check_call(cmd,shell=True)
	groups = read_groups(group)
	#Collect all organism and select single copy genes orthogroups (core genome)
	organism = get_organism(groups)
	print "\n%s species and %s orthogroups detected.\n" %(str(len(organism)),str(len(groups)))
	groups = get_core(groups,organism)
	print "%s selected core orthogroups for running codeML.\n" %(str(len(groups)))
	#Get sequences from multifasta file
	results = select_genes(groups,organism,infile)
	print_general(results)
	#cladeAresults = select_triplets(cladeA,fixbranchA,groups)
	#print_general(cladeAresults,"DoEn")

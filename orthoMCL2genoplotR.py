#!/usr/bin/python
import re
import argparse
#===============================================================================
#    Author: Diego Santos-Garcia
#    Contact: Contact the author at diego.santos@mail.huji.ac.il or diego.santos.garcia@protonmail.com
#
#    COPYRIGHT: Copyright (C) 2020  Diego Santos-Garcia.
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

def get_orgs(groups):
	"""(str) -> dict

	Read the output of OrthoMCL and return a dictionary of organims.

	>>>POR0001: BtBH|B186_001 BtQV|PAQ_001 CaCD|CRDC_00415 CaHC|A353_086
	BtBH, BtQV, CaCD, CaHC
	"""
	fhand = open(groups,"r")
	orgs = {}
	for line in fhand:
		line = line.split()
		for element in line:
			if re.search('\w+\|',element):
				element = element.split("|")
				orgs[element[0]] = 1
	return orgs


def get_core(orgs,order):
	"""(dict) -> str

	Return core orthologous groups from OrthoMCL output for all the organims.
	"""
	fhand = open(groups,"r")
	fout = open(output,"aw")
	orgs_list = order
	for line in fhand:
		if all(letter in line for letter in orgs_list):
			fout.write(line)
	fout.close()

def get_dict(ortho_tab,artemis_file,genes_dict):
	"""Return a dictionary with the locus tag of each organism for every orthologous group from a list of Core genes from 	getCorepy.
	>>>POR0001: BtBH|B186_001 BtQV|PAQ_001 CaCD|CRDC_00415 CaHC|A353_086 CaPV|CRP_062 PADV|PAD_001 PLFV|PLF_001 PTVM|PalTV_001 PTVV|PTV_001
	"""

	for line in ortho_tab:
		line = line.strip().split(" ")
		#print line
		#print "Generating comparision table for gene" + line[0]
		group = line[0].replace(':',"")
		line.pop(0)
		genes_dict[group]={}
		for item in line:
			item = item.split("|")
			#print item
			genes_dict[group][item[0]] = {item[1]:[]}
	return genes_dict





def get_location(genes_dict,artemis_file):
	"""Read a concatenated file with the fasta headers produced by Artemis for all the CDS used in OrthoMCL and fill the genes_dict.
	>>>B186_001 dnaK chaperone protein DnaK 1:1926 forward MW:70656
	"""

	for line in artemis_file:
		locus = re.search("^\w+\s",line)
#		print locus.group()
		locus = locus.group().strip()
		for key in genes_dict:
			for item in genes_dict[key]:
				for orgs in genes_dict[key][item]:
					if locus in orgs:
						genes_dict[key][item][orgs] = get_start_end(line)
	#print genes_dict
	return genes_dict


def get_start_end(line):
	"""" Return a list with the location of the gene in the chormosome"""

	location = re.search('[0-9]+:[0-9]+',line)
	location = (location.group()).split(":")
	if 'reverse' in line:
		location = [location[1],location[0]]
	else:
		pass
	return location


def print_table(genes_dict,order):
	"""Print a fancy table suitable for genoPlotR"""

	sorted_keys = sorted(genes_dict.keys())
	first = 0
	second = 1
	for item in order:
		if second == len(order):
			break
		else:
			fout = open(order[first] + 'vs' + order[second],"aw")
			fout.write('locus1\tlocus2\tstart1\tend1\tstart2\tend2\tdirection\tcol\n')
			for key in sorted_keys:
				locus1 = genes_dict[key][order[first]].values()
				start1 = int(locus1[0][0])
				end1 = int(locus1[0][1])
				locus2 = genes_dict[key][order[second]].values()
				#print "".join(genes_dict[key][order[first]].keys()),"".join(genes_dict[key][order[second]].keys())
				start2 = int(locus2[0][0])
				end2 = int(locus2[0][1])
				if (start1 < end1 and start2 < end2) or (start1 > end1 and start2 > end2):
					frame = 1
					color = 'grey'
				else:
					frame = -1
					color = 'blue'
				if start1 > end1:
					temp = end1
					end1 = start1
					start1 = temp
				if start2 > end2:
					temp = end2
					end2 = start2
					start2 = temp  
				fout.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %("".join(genes_dict[key][order[first]].keys()),"".join(genes_dict[key][order[second]].keys()),start1,end1,start2,end2,frame,color))
			first += 1
			second += 1
		fout.close()

def print_dna_seg(genes_dict):

	sorted_keys = sorted(genes_dict.keys())
	for item in order:
		fout = open(item + '_dna_seg',"aw")
		fout.write('name	start	end	strand	col	lty	lwd	pch	cex	gene_type	fill\n')
		for key in sorted_keys:
			location = genes_dict[key][item].values()
			if int(location[0][0]) < int(location[0][1]):
				frame = 1
			else:
				frame = -1
			if frame == 1:
				color = 'orange'
				fout.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %(genes_dict[key][item].keys()[0],location[0][0],location[0][1],frame,color,'1','1','8','1','side_blocks',color))
			if frame == -1:
				color = 'red'
				fout.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %(genes_dict[key][item].keys()[0],location[0][1],location[0][0],frame,color,'1','1','8','1','side_blocks',color))

		fout.close()

def get_grimm(genes_dict,order,grimm,numbers):
	for item in order:
		grimm[item] = {}
		sorted_keys = sorted(genes_dict.keys())
		for key in sorted_keys:
			ID = genes_dict[key][item].keys()
			#gene_number = re.search('[0-9]+$',key)  #####Modify for renumbering clusters name, only can run wiht core genome
			#gene_number = int(gene_number.group())
			location = genes_dict[key][item].values()[0]
			if int(location[0]) < int(location[1]):
				grimm[item][ID[0]] = str(numbers[key])
			else:
				grimm[item][ID[0]] = str(-numbers[key])
	return grimm

def get_numbers(genes_dict):
	i = 1
	numbers = {}
	sorted_keys = sorted(genes_dict.keys())
	for key in sorted_keys:
		numbers.update({key:''})
	numbers_sort = sorted(numbers.keys())
	for key in numbers_sort:
		numbers[key] = i
		i += 1
	return numbers

def print_grimm(grimm):

	fout = open("grimm_core.txt","aw")
	for key in grimm.keys():
		temp = []
		for item in sorted(grimm[key].keys()):
			temp.append(grimm[key][item])
		temp = ' '.join(temp)
		fout.write('>' + key + '\n' + temp + '$' + '\n')
	fout.close()


if __name__ == "__main__":

	parser = argparse.ArgumentParser(description="Give the core genes from the OrthoMCL output.")
	parser.add_argument("-g","--groups",action='store',required = True, help = "OrthoMCL groups output file, i.e. groups_1.5.txt")
	parser.add_argument("-o","--output",action='store',required = True, help ="Name for the output file containing the core genes")
	parser.add_argument("-l","--list",action='store',required = True, help ="Tabular file with four letters name of the desired organism to extract the core genome and in the desired order for plotting with genoPlotR")
	parser.add_argument("-a","--artemis",action='store',required = True, help ="Headers with gene names extracted with artemis, i.e. CRDC_00955 rpsG 30S ribosomal protein S7 152259:152729 reverse MW:18566")
	args = parser.parse_args()
	groups = args.groups
	output = args.output
	artemis = args.artemis
	order_list = open(args.list,"r")
	for line in order_list:
		order = line.strip().split("\t")
	orgs = get_orgs(groups)
	get_core(orgs,order)
	ortho_tab = open(output,"r")
	artemis_file = open(artemis,"r")  ####File with header of a fasta generated with artemis
	genes_dict = {}
	genes_dict = get_dict(ortho_tab,artemis_file,genes_dict)
	genes_dict = get_location(genes_dict,artemis_file)
	print_table(genes_dict,order)
	print_dna_seg(genes_dict)
	numbers = get_numbers(genes_dict)
	grimm = {}
	grimm = get_grimm(genes_dict,order,grimm,numbers)
	print_grimm(grimm)

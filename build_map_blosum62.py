# coding=utf-8
import os
import sys
import numpy as np
import math
import re
import xlsxwriter as xw
import matplotlib
import requests
import Bio.PDB as biopdb
from lxml import etree
from matplotlib import pyplot as plt
from matplotlib import colors

#execution modifiers
normalize = True #reports frequencies instead of absolute numbers of mutations
normalization_multiplier = 10000 #frequencies are multiplied by this number for their visualization over a "standard" sample of this number of mutations
integer_normalization = True #decides if to approximate to integer after normalizing and multiplitying by "normalization_multiplier"

#mutation object class
class Mutation:
	
	#constructor
	def __init__(self, SequentialID, name, nucl_var, prot_var):
		self.SeqID = SequentialID			#sequential ID of mutation assigned by ClinVar in the output file relative to our search (Ours)
		self.name = name					#name of mutation (ClinVar)
		self.variation_nucleic = nucl_var	#nucleic acid variation (ClinVar)
		self.variation_protein = prot_var	#protein variation (ClinVar)
		self.genes = None					#genes involved in the mutation (ClinVar)
		self.protein_change = None			#type of protein change (ClinVar)
		self.condition = None				#condition associated to this mutation (ClinVar)
		self.significance = None			#clinical significance of this mutation (ClinVar)
		self.review = None					#review status of the mutation entry (ClinVar)
		self.chromosome = None				#chromosome on which the mutation occurs (ClinVar)
		self.location = None				#location in the chromosome sequence (ClinVar)
		self.accession = None				#accession key (ClinVar)
		self.ClinVarID = None				#ID of mutation entry (ClinVar)
		self.canonical_spdi = None			#canonical SPDI notation for this mutation (ClinVar)
		self.pdb = None						#Protein Data Bank (PDB) structure file associated to this mutation (VarMap)
		self.chain = None					#chain on which the mutation occurs in the PDB structure (VarMap)
		self.res_name = None				#residue name of the aminoacid found at the mutation point in the PDB structure (VarMap)
		self.res_number = None				#residue number at which the mutation occurs in the PDB structure (VarMap)
		self.dna_complex_type = None		#type of DNA complex, only for mutations on protein-DNA complexes (Ours)	
		self.interface_type = None			#list of type of interfaces this mutation occurs in if any, empty or None otherwise (Ours)
		self.entries = list()				#PISA interface entries this mutation was found in (Ours)
		self.is_on_surface = None			#boolean variable telling if the mutation occurs on the protein surface (Ours)

	#method to add an interface entry to this mutation
	def add_interface_entry(self, interface_type, interacting_chain):
		self.entries.append({"it": interface_type, "ic": interacting_chain})

	#method to refine an interface entry, adding information regarding the interacting residue and the interacting atom
	def refine_interface_entry(self, index, interacting_residue, interacting_atom):
		self.entries[index]["ir"] = interacting_residue
		self.entries[index]["ia"] = interacting_atom

	#method to add an interface type to the list of interface types this mutation is involved in
	def add_interface_type(self, it):
		if self.interface_type is None:
			self.interface_type = [it]
		else:
			self.interface_type.append(it)

	#method to create a summary string that describes the types of interface the mutation is involved in
	def get_interface_type_string(self):
		if self.interface_type is None:
			return '-'
		res = ""
		if "DNA" in self.interface_type:
			res = res+"D"
		else:
			res = res+"-"
		if "RNA" in self.interface_type:
			res = res+"R"
		else:
			res = res+"-"
		if "Ligand" in self.interface_type:
			res = res+"L"
		else:
			res = res+"-"
		if "Protein" in self.interface_type:
			res = res+"P"
		else:
			res = res+"-"
		return res

	#method to create a summary string of the protein-dna complex types of the mutations
	def get_dnac_type_string(self):
		if self.dna_complex_type is None:
			return '-'
		res = ""
		if "SP" in self.dna_complex_type:
			res = res+"S"
		else:
			res = res+"-"
		if "TF" in self.dna_complex_type:
			res = res+"T"
		else:
			res = res+"-"
		if "DE" in self.dna_complex_type:
			res = res+"E"
		else:
			res = res+"-"
		return res
	
	#method to format the residue number string
	def get_res_num_string(self):
		if int(self.res_number) > 999:
			return self.res_number+"\t"
		else:
			return self.res_number+"\t\t" 
	
	#update input dicts accounting for Arg and Lys mutations
	def UpdatePisaDicts(self, pisa_dicts, current_key):
		#retrieve base aminoacid for this mutation
		base_aminoacid =  self.variation_protein[:3]
		#if the three letter "base aminoacid" code is not represented in the dict keys, report an error
		if base_aminoacid not in pisa_dicts.keys():
			sys.exit("ERROR: unrecognized aminoacid "+base_aminoacid)
		#otherwise update the corresponding class count
		pisa_dicts[base_aminoacid][current_key] = pisa_dicts[base_aminoacid][current_key] + 1   

#parameters
path_input = "Data/clinvar_result.txt"
path_plot = "Data/mutation_plot.png"
path_dist_matrix =  "Data/dist_matrix.txt"
path_mutation_matrix = "Data/mutation_matrix.txt"
path_output = "Data/mutation_map.xlsx"
path_blosum62 = "Data/Blosum62/Blosum62.csv"
path_report = "Data/report.txt"
colours = np.array([ [1.0, 1.0, 1.0], [0.0, 0.9, 0.9], [0.9, 0.9, 0.0], [0.9, 0.0, 0.9] ])
aminoacids = ["A", "C", "D", "E", "F", "G", "H",  "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
aminoacid_lookup = {"A":0, "C":1, "D":2, "E":3, "F":4, "G":5, "H":6,  "I":7, "K":8, "L":9, "M":10, "N":11, "P":12, "Q":13, "R":14, "S":15, "T":16, "V":17, "W":18, "Y":19}
dict_1L = {"A": "Ala", "C": "Cys", "D": "Asp", "E": "Glu", "F": "Phe", "G": "Gly", "H": "His",  "I": "Ile", "K": "Lys", "L": "Leu", "M": "Met", "N": "Asn", "P": "Pro", "Q": "Gln", "R": "Arg", "S": "Ser", "T": "Thr", "V": "Val", "W": "Trp", "Y": "Tyr"}
dict_3L = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V"}
codons = {"TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "TAT": "Y", "TAC": "Y", "TAA": "STOP", "TAG": "STOP", "TGT": "C", "TGC": "C", "TGA": "STOP", "TGG": "W", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P", "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M", "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R", "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"}
codon_lists = { "A": ["GCT", "GCC", "GCA", "GCG"],
				"C": ["TGT", "TGC"],
				"D": ["GAT", "GAC"],
				"E": ["GAA", "GAG"],
				"F": ["TTT", "TTC"],
				"G": ["GGT", "GGC", "GGA", "GGG"],
				"H": ["CAT", "CAC"],
				"I": ["ATT", "ATC", "ATA"],
				"K": ["AAA", "AAG"],
				"L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
				"M": ["ATG"],
				"N": ["AAT", "AAC"],
				"P": ["CCT", "CCC", "CCA", "CCG"],
				"Q": ["CAA", "CAG"],
				"R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
				"S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
				"T": ["ACT", "ACC", "ACA", "ACG"],
				"V": ["GTT", "GTC", "GTA", "GTG"],
				"W": ["TGG"],
				"Y": ["TAT", "TAC"] }

#execution start

#build distance matrix between aminoacids
blosum62 = np.loadtxt(path_blosum62, dtype=int, delimiter=",")
print("Building matrix of Blosum62 distances between aminoacids")
dist_matrix = np.zeros((len(aminoacids),len(aminoacids)), dtype = int)
for i in range(len(aminoacids)):
	for j in range(len(aminoacids)):
		dist = 0
		if i == j:
			continue
		#determine distance according to blosum62 scores
		if blosum62[i][j] > 0:
			dist = 1
		elif blosum62[i][j] in [-1, 0]:
			dist = 2
		else:
			dist = 3
		#update dist_matrix
		dist_matrix[i][j] = dist
#save matrix to file
np.savetxt(path_dist_matrix, dist_matrix, fmt="%d")
#build custom colormap
cm = matplotlib.colors.ListedColormap(colours , "traffic_light")
#plot distance matrix with matplotlib
figure, ax = plt.subplots()
plt.imshow(dist_matrix, cmap = cm)
plt.colorbar()
ax.set_xticks(np.arange(len(aminoacids)))
ax.set_yticks(np.arange(len(aminoacids)))
ax.set_xticklabels(aminoacids)
ax.set_yticklabels(aminoacids)
plt.savefig(path_plot)

#load input txt file
print("Reading list of ClinVar mutations from file")
in_file = open(path_input, 'r')
in_text = in_file.read()
in_file.close()
lines = in_text.splitlines(in_text.count("\n"))
#build mutation list
mutations = list()
j = -1
#parse input lines
for i in range(1, len(lines)):
	#skip empty lines 
	if len(lines[i])<=1:
		continue
	#check for Sequential Number (declaration line)
	if re.match(r"[\d]+", lines[i]):
		#split the line in three around ":" charachters
		cells = re.split(":", lines[i])
		#name corresponds to the second cell
		name = cells[1].strip()
		#split third cell in three around "." charachters
		subcells = re.split("\.", cells[2])
		#match the nucleic variation in the second subcell: if a match is found, save the nucleic acid and protein variations
		if re.match(r"\d+\w.\w", subcells[1]):
			m = re.match(r"\d+\w.\w", subcells[1])
			nucl_var = m.group(0)
			#the protein variation corresponds to the last subcell (without last parenthesis and \n), if present
			if len(subcells)>=3:
				prot_var = subcells[2]
				prot_var = prot_var[:-2]
			else:
				prot_var = None
		#otherwise copy the nucleic variation string and skip the protein variation
		else:
			nucl_var = subcells[1].strip()
			prot_var = None
		#increase j
		j = j+1
		#create new mutation object	and append it to the list of mutations
		mutations.append( Mutation( j, name, nucl_var, prot_var ) )
		#skip other possible operations
		continue
	#split the line in cells around ":" charachters
	cells = re.split(":", lines[i])
	#switch on field name
	if cells[0] == "Gene(s)":
		mutations[j].genes = cells[1].strip()
	elif cells[0] == "Protein change":
		mutations[j].protein_change = cells[1].strip()
	elif cells[0] == "Condition(s)":
		mutations[j].condition = cells[1].strip()
	elif cells[0] == "Clinical significance":
		subcells = re.split("\(", cells[1])
		mutations[j].significance = subcells[0].split()
	elif cells[0] == "Review status":
		mutations[j].review = cells[1].strip()
	elif cells[0] == "Chromosome":
		mutations[j].chromosome = cells[1].strip()
	elif cells[0] == "Location  (GRCh38)":
		mutations[j].location = cells[1].strip()
	elif cells[0] == "Accession":
		mutations[j].accession = cells[1].strip()
	elif cells[0] == "ID":
		mutations[j].ClinVarID = cells[1].strip()
	elif cells[0] == "Canonical SPDI":
		mutations[j].canonical_spdi = cells[1].strip()+":"+cells[2]+":"+cells[3]+":"+cells[4]

#build filtered list of mutations
print("Analyzing mutations")
filtered_mutations = list()
#build count matrix
mutation_matrix = np.zeros((len(aminoacids),len(aminoacids)), dtype = int)
#iterate over all the parsed mutations
for m in mutations:
	#check for missing or invalid protein change descriptions
	if m.variation_protein is None:
		continue
	#check review status (uncomment the two lines below to discard all the mutations without review)
	#if re.match("no assertion", m.review):
	#	continue
	#check for unwanted annotations
	if m.variation_protein[-1] == "=":
		continue
	#otherwise extract source and destination
	src = m.variation_protein[:3]
	dst = m.variation_protein[-3:]
	#skip stops, unknown symbols or special aminoacids
	if src in ["Ter", "Xaa"] or dst in ["Ter", "Xaa"]:
		continue
	#update mutation matrix
	mutation_matrix[aminoacid_lookup[dict_3L[src]]][aminoacid_lookup[dict_3L[dst]]] += 1
	#update filtered list of mutations
	filtered_mutations.append(m)

#save mutation matrix to file
np.savetxt(path_mutation_matrix, mutation_matrix, fmt="%d")

#calculate total mutations from each aminoacid to others
total_mutations_from = np.sum(mutation_matrix, axis=1)
#calculate total mutations to each aminoacid from others
total_mutations_to = np.sum(mutation_matrix, axis=0)
#calculate global total
number_of_mutations = np.sum(total_mutations_from)

#normalize quantities of mutations to a fixed standard sample size
if normalize:
	#transform the quantities into frequencies
	mutation_matrix = np.divide(mutation_matrix, number_of_mutations)
	total_mutations_from = np.divide(total_mutations_from, number_of_mutations)
	total_mutations_to = np.divide(total_mutations_to, number_of_mutations)
	#multiply frequencies by the standard sample size
	mutation_matrix = np.multiply(mutation_matrix, normalization_multiplier)
	total_mutations_from = np.multiply(total_mutations_from, normalization_multiplier)
	total_mutations_to = np.multiply(total_mutations_to, normalization_multiplier)
	number_of_mutations = normalization_multiplier
	#approximate to closest integer if needed
	if integer_normalization:
		mutation_matrix = mutation_matrix.astype(int)
		total_mutations_from = total_mutations_from.astype(int)
		total_mutations_to = total_mutations_to.astype(int)
	
#save the mutation matrix to .xlsx
print("Writing matrix of aminoacid mutations to file")
workspace = xw.Workbook(path_output)
page = workspace.add_worksheet()
format_bold = workspace.add_format( {"bold": True} )
format_colors = [workspace.add_format(), workspace.add_format({"bg_color": "#00FF00"}), workspace.add_format({"bg_color": "#FFFF00"}), workspace.add_format({"bg_color": "#FF0000"})]
#define column and row locations
col_head = "B"
cols = ["C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V"]
col_tail = "X"
row_head = "2"
rows = [str(i) for i in range(3,23)]
row_tail = "24"
#build header row
for i in range(len(aminoacids)):
	page.write(cols[i]+row_head, aminoacids[i], format_bold)
page.write(col_tail+row_head, "total", format_bold)
#build header column
for i in range(len(aminoacids)):
	page.write(col_head+rows[i], aminoacids[i], format_bold)
page.write(col_head+row_tail, "total", format_bold)
#build table
for i in range(len(aminoacids)):
	for j in range(len(aminoacids)):
		page.write(cols[j]+rows[i], mutation_matrix[i][j], format_colors[dist_matrix[i][j]])
#build tail row
for i in range(len(aminoacids)):
	page.write_formula(cols[i]+row_tail, '=SUM('+cols[i]+rows[0]+':'+cols[i]+rows[-1]+')', format_colors[0], total_mutations_to[i])
#build tail column
for i in range(len(aminoacids)):
	page.write_formula(col_tail+rows[i], '=SUM('+cols[0]+rows[i]+':'+cols[-1]+rows[i]+')', format_colors[0], total_mutations_from[i])
#write sum of totals
page.write_formula(col_tail+row_tail, '=SUM('+col_tail+rows[0]+':'+col_tail+rows[-1]+')', format_colors[0], np.sum(number_of_mutations))
#save .xlsx file
workspace.close()

#write report
print("Writing final report")
out_file = open(path_report, "w")
out_file.write("Number of mutations from ClinVar : "+str(len(mutations))+"\n")
out_file.write("Valid mutations analyzed : "+str(len(filtered_mutations))+"\n")
out_file.write("Of which R/X : "+str(total_mutations_from[aminoacid_lookup['R']])+"\n")
out_file.write("Of which G/X : "+str(total_mutations_from[aminoacid_lookup['G']])+"\n\n\n")
print("Execution terminated succesfully")


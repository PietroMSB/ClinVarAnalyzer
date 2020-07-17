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
check_and_download_files = True  #switch to False if all the structure and interface files were already downloaded, in order to speed up execution 

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

	#method to refine an interface entry, adding information regarding the interacting residue and the interacting_atom
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

#xml tags
tag_root = "pisa_interfaces"
tag_subroot = "pdb_entry"
tag_pdb = "pdb_code"
tag_status = "status"
tag_interface = "interface"
tag_area = "int_area"
tag_molecule = "molecule"
tag_symmetry = "symop"
tag_class = "class"
tag_chain = "residues"
tag_chain_id = "chain_id"
tag_energy = "solv_en"
tag_code = "name"
tag_residue = "residue"
tag_serial = "ser_no"
tag_sequential = "seq_num"
class_protein = "Protein"
class_ligand = "Ligand"
class_DNA = "DNA"
class_RNA = "RNA"
sym_op = "x,y,z"


#functions for parsing XML files

#returns the first child with tag "tag" of node "node"
def getChild(node, tag):
	for i in range(len(node)):
		if(node[i].tag == tag):
			return node[i]
	return None

#returns a list which contains all the children with tag "tag" of node "node"
def getChildren(node, tag):
	res=[]
	for i in range(len(node)):
		if(node[i].tag == tag):
			res.append(node[i])
	return res

#gets a "Molecule" node's class
def getMoleculeClass(molecule):
	class_node = getChild(molecule,tag_class)
	if(class_node == None):
		#tag not found, xml acquisition error
		print("Tag class not found for molecule.")
		return "" 
	else:
		return str(class_node.text)

#checks the symmetry operator applied to this interface (to avoid interfaces which are a side result of the crystallization process)
def checkSymOp(molecule):
	s = str(getChild(molecule,tag_symmetry).text).lower()
	return (s == sym_op)

#returns true if an interface met the parametres and was analyzed, false otherwise
def analyzeInterface(interface_xml, this_pdb_mutations):
	#extract list of residue numbers to check
	res_numbers = list()
	for m in this_pdb_mutations:
		res_numbers.append(int(m.res_number))
	#retrieve interface entries
	interface_entries = getChildren(interface_xml,tag_interface)
	#iterate over all the interface entries
	for k in range(len(interface_entries)):
		area=0
		molecules=[]
		current = interface_entries[k]
		#search for the interface's surface value and for the two molecules involved
		area=float(getChild(current,tag_area).text)
		molecules=getChildren(current,tag_molecule)
		if(len(molecules)!=2):
			print("Error: number of interface molecules not equal to 2.")
			continue
		#symmetry operator filter (both molecules must have symmetry operator equal to "sym_op")
		if( not ((checkSymOp(molecules[0])) and (checkSymOp(molecules[1])))):
			continue
		#determine interface type
		entry_type = None
		#protein-ligand interface identifier
		if(((getMoleculeClass(molecules[0])==class_ligand) and (getMoleculeClass(molecules[1])==class_protein)) or ((getMoleculeClass(molecules[0])==class_protein) and (getMoleculeClass(molecules[1])==class_ligand))):
			entry_type = "Ligand"
		#protein-protein interface identifier
		elif((getMoleculeClass(molecules[0])==class_protein) and (getMoleculeClass(molecules[1])==class_protein)):
			entry_type = "Protein"
		#protein-DNA interface identifier
		elif(((getMoleculeClass(molecules[0])==class_DNA) and (getMoleculeClass(molecules[1])==class_protein)) or ((getMoleculeClass(molecules[0])==class_protein) and (getMoleculeClass(molecules[1])==class_DNA))):
			entry_type = "DNA"
		#protein-RNA interface identifier
		elif(((getMoleculeClass(molecules[0])==class_RNA) and (getMoleculeClass(molecules[1])==class_protein)) or ((getMoleculeClass(molecules[0])==class_protein) and (getMoleculeClass(molecules[1])==class_RNA))):
			entry_type = "RNA"
		#otherwise, this interface entry is not interesting
		else:
			continue
		#scan the two molecules
		for l in range(len(molecules)):
			#avoid non-protein molecules
			if getMoleculeClass(molecules[l])==class_protein:
				#read chain identifier
				chain_id = getChild(molecules[l], tag_chain_id).text
				#retrieve chain node
				chain = getChild(molecules[l], tag_chain)
				#iterate over the chain sequence and search for the mutation residues
				for r in range(len(chain)):
					#skip text nodes between residues
					if chain[r].tag!=tag_residue: 
						continue
					#find the residues involved in at least one mutation
					if int(chain[r].find(tag_sequential).text) in res_numbers:
						marker = None
						#check if the residue is involved in the interface
						if float(chain[r].find(tag_energy).text) != 0:
							marker = entry_type
						else:
							marker = "Out"
						#scan the residue against all the mutations
						for m in this_pdb_mutations:
							if int(m.res_number) == int(chain[r].find(tag_sequential).text) and m.chain == chain_id :
								m.add_interface_type(marker)
								#add the interface entry to the mutation record (avoiding "Out of interfaces" entries)
								if marker != "Out":	
									#retrieve interacting chain's id
									interacting_chain_id = None
									if l==0:
										interacting_chain_id = getChild(molecules[1], tag_chain_id).text
									else:
										interacting_chain_id = getChild(molecules[0], tag_chain_id).text
									#add the entry (interface type, interacting chain) to the mutation record
									m.add_interface_entry(marker, interacting_chain_id)		


#parameters
threshold_percent_sasa = 0.2
path_input = "Data/clinvar_result.txt"
path_plot = "Data/mutation_plot.png"
path_dist_matrix =  "Data/dist_matrix.txt"
path_mutation_matrix = "Data/mutation_matrix.txt"
path_output = "Data/mutation_map.xlsx"
path_varmap_input = "Data/VarMapInput.txt"
path_varmap_output = "Data/VarMap_results.tsv"
path_report = "Data/report.txt"
path_xml = "Data/PISA/Interfaces/"
path_pdb = "Data/PISA/Structures/"
path_pops = "Data/PISA/POPS_output/"
path_list_pisa = "Data/Lists/PISA_unique_pdbs.txt"
pisa_base_url = "http://www.ebi.ac.uk/pdbe/pisa/cgi-bin/interfaces.pisa?"
structure_base_url = "https://files.rcsb.org/download/"
colours = np.array([ [1.0, 1.0, 1.0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0], [1.0, 0.0, 0.0] ])
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
print("Building matrix of minimum codon distances between aminoacids")
dist_matrix = np.zeros((len(aminoacids),len(aminoacids)), dtype = int)
for i in range(len(aminoacids)):
	for j in range(len(aminoacids)):
		#calculate distance as the minimum number of mutations
		min_dist = 4
		for ci in codon_lists[aminoacids[i]]:
			for cj in codon_lists[aminoacids[j]]:
				#calculate distance between ci and cj
				dist = 0
				for k in range(3):
					if ci[k] != cj[k]:
						dist = dist + 1
				#check minimum distance
				if dist < min_dist:
					min_dist = dist
		#update dist_matrix
		dist_matrix[i][j] = min_dist
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
		if len(cells) >=5:			
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
	if len(m.variation_protein)==0:
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
	#correct wrong protein variation entries with protein change entries, if the latter are adequate
	if src not in dict_3L.keys() or dst not in dict_3L.keys():
		#if protein change does not match the pattern, ignore this badly written mutation entry
		first_match = re.match("[A-Z][\d]+[A-Z]", m.protein_change)
		if first_match is None:
			continue
		#otherwise recover protein variation from protein change
		match_str = first_match.group(0)
		src = dict_1L[match_str[0]]
		dst = dict_1L[match_str[-1]]
		num = match_str[1:-1]
		m.variation_protein = src+num+dst
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
page.write_formula(col_tail+row_tail, '=SUM('+col_tail+rows[0]+':'+col_tail+rows[-1]+')', format_colors[0], np.sum(total_mutations_from))
#save .xlsx file
workspace.close()

#produce list of variants for VarMap, in order to retrieve the PDB file of each entry
print("Building list of variants to be retrieved from VarMap")
out_file = open(path_varmap_input, 'w')
#write a \n separated list with one entry for each filtered mutation
for m in filtered_mutations:
	#process nucleic variation to extract the two bases
	src = m.variation_nucleic[-3]
	dst = m.variation_nucleic[-1]
	#write a space separated record reporting: chromosome, genomic coordinates, source base, variant base 
	out_file.write(m.chromosome+" "+m.location+" "+src+" "+dst+"\n")
#close file
out_file.close()

### UPLOAD FILE MANUALLY ON VAR MAP AND WAIT FOR RESULTS

#read VarMap output file
print("Reading VarMap output")
in_file = open(path_varmap_output, "r")
in_text = in_file.read()
in_file.close()
in_lines = in_text.splitlines(in_text.count("\n"))
del in_text
#cycle over input lines (skipping the header)
for i in range(1,len(in_lines)):
	#split line into cells around '\t' characters
	cells = re.split("[\t]", in_lines[i])
	### DEBUG START ###
	'''
	#check parsing quality
	for j in range(len(cells)):
		print(str(j)+"\t|"+cells[j].strip()+"|")
	sys.exit()
	'''
	### DEBUG STOP ###
	#check for correspondence between submitted record (chromosome, location, base nucleotide, mutant nucleotide) and VarMap record
	if filtered_mutations[i-1].chromosome != cells[0].strip() or filtered_mutations[i-1].location != cells[1].strip() or filtered_mutations[i-1].variation_nucleic[-3] != cells[2].strip() or filtered_mutations[i-1].variation_nucleic[-1] != cells[3].strip():
		print(filtered_mutations[i-1].chromosome+" "+filtered_mutations[i-1].location+" "+filtered_mutations[i-1].variation_nucleic[-3]+" "+filtered_mutations[i-1].variation_nucleic[-1])
		print(cells[0].strip()+" "+cells[1].strip()+" "+cells[2].strip()+" "+cells[3].strip())
		sys.exit("ERROR: record "+str(i)+" does not correspond to the submitted mutation")
	#check for existence of a pdb file
	if cells[33].strip() == "TRUE" and len(cells[35].strip()) == 4:
		filtered_mutations[i-1].pdb = cells[35].strip()
		#acquire chain identifier
		filtered_mutations[i-1].chain = cells[36].strip()
		#acquire residue name in pdb file
		filtered_mutations[i-1].res_name = cells[45].strip()
		#acquire residue number in pdb file
		filtered_mutations[i-1].res_number = re.match("[\d]+",cells[46].strip()).group()
	#otherwise mark the "pdb" field as "NO_PDB_FILE", leaving the others as None
	else:
		filtered_mutations[i-1].pdb = "NO PDB FILE"

#count how many filtered mutations have a pdb file, create two lists of mutations with source R and G respectively, build list of pdb files
print("Building list of pdb structure files to be retrieved from PDBe-PISA")
raw_covered_mutations = list()
for m in filtered_mutations:
	#skip mutations without pdb file
	if m.pdb == "NO PDB FILE":
		continue
	#count pdb file occurrence
	raw_covered_mutations.append(m)	

#cut mutations for which the residue at the corresponding location in the PDB structure file does not correspond to "src" nor to "dst"
covered_mutations=list()
R_X_mutations = list()
G_X_mutations = list()
unique_pdbs = list()
for m in raw_covered_mutations:
	#eliminate mutations with residue mismatch
	if (m.res_name.lower() != m.variation_protein[:3].lower()) and (m.res_name.lower() != m.variation_protein[-3:].lower()):
		continue
	#otherwise, append the mutation to the final list
	covered_mutations.append(m)
	#check uniqueness of pdb file
	if m.pdb not in unique_pdbs:
		unique_pdbs.append(m.pdb)
	#check if this mutation has R as a source base
	if m.variation_protein[:3] == 'Arg':
		R_X_mutations.append(m)
	#check if this mutation has K as a source base
	if m.variation_protein[:3] == 'Gly':
		G_X_mutations.append(m)
res_number_conflicts = len(raw_covered_mutations)-len(covered_mutations)

#save list of unique pdbs to file
np.savetxt(path_list_pisa, unique_pdbs, fmt="%s", delimiter=',')

#download the interfaces of each pdb structure in the list of unique pdbs from PDBe-PISA
if check_and_download_files:
	print("Connecting to PDBe-PISA")
	for i in range(len(unique_pdbs)):
		#build path string
		filepath = path_xml+unique_pdbs[i].lower()+".xml"
		#check if the file already exists
		if os.path.exists(filepath):
			print("Checking interface file: "+unique_pdbs[i]+".xml ( "+str(i+1)+" of "+str(len(unique_pdbs))+" )", end='\r')
			#check file integrity
			try:
				dom_check = etree.parse(filepath)
				continue
			#catch keyboard interrupts
			except KeyboardInterrupt:
				print("Execution terminated by keyboard")
				sys.exit()
			#catch integrity check failures and, in case, just download the file again
			except:
				print("Integrity check failed", end='\r')
		#download the file (if it was missing or failed the integrity check)
		print("Retrieving interface file for "+unique_pdbs[i]+" ( "+str(i+1)+" of "+str(len(unique_pdbs))+" )", end='\r')
		response=requests.get(pisa_base_url+unique_pdbs[i].lower())
		if(response.status_code==200):
			outfile = open(filepath,'wb')
			outfile.write(response.content)
			outfile.close()
	print("")

#scan all the interface entries from the XML files downloaded from PDBe-PISA, assigning a category to each PISA-covered mutation
empty_xml_files = 0
for i in range(len(unique_pdbs)):
	print("Analyzing interface file: "+unique_pdbs[i]+".xml ( "+str(i+1)+" of "+str(len(unique_pdbs))+" )", end='\r')
	#build the list of mutations occurring on the same PDB structure
	this_pdb_mutations = list()
	for m in covered_mutations:
		if m.pdb == unique_pdbs[i]:
			this_pdb_mutations.append(m)
	#try parsing the xml file
	try:
		DOM = etree.parse(path_xml+unique_pdbs[i]+".xml")
	#catch exceptions and remove bad xml files 
	except Exception as exc:
		print(exc)
		print("File "+unique_pdbs[i]+".xml was acquired incorrectly")
		os.remove(path_xml+unique_pdbs[i]+".xml")
		print("Removed file: "+unique_pdbs[i]+".xml")
		sys.exit()
	root = DOM.getroot()
	pdb_entry = getChild(root, tag_subroot)
	#if PDB entry status is "Ok", parse the xml tree
	if getChild(pdb_entry,tag_status).text == "Ok":
		#analyze all the entries in this interface file, assigning interface types to each mutation on the same PDB
		analyzeInterface(pdb_entry, this_pdb_mutations)
	#otherwise count the file as empty
	else:
		empty_xml_files += 1
print("")

#count the mutations by interface type
print("Extracting information from mutation data")
pisa_class_counts = {"Protein":0, "Ligand":0, "DNA":0, "RNA":0, "Out":0, "NotFound":0}
#build a list containing each pdb from unique_pdbs having at least one mutation in at least one interface
pdbs_with_interface_mutations = list()
#build a class dictionary for each aminoacid, grouping the mutations by "base aminoacid"
pisa_dicts = dict()
for a in aminoacids:
	pisa_dicts[dict_1L[a]] = {"Protein":0, "Ligand":0, "DNA":0, "RNA":0, "Out":0, "NotFound":0}
#scan all the mutations, filling the dictionaries
for m in covered_mutations:
	#check if the mutation residue was found in at least one interface entry (also if not part of the interface)
	if m.interface_type is None:
		pisa_class_counts["NotFound"] = pisa_class_counts["NotFound"]+1
		m.UpdatePisaDicts(pisa_dicts, "NotFound")
		continue
	#check for classes
	no_class = True
	if "Protein" in m.interface_type:
		pisa_class_counts["Protein"] = pisa_class_counts["Protein"]+1
		m.UpdatePisaDicts(pisa_dicts, "Protein")
		no_class = False
	if "Ligand" in m.interface_type:
		pisa_class_counts["Ligand"] = pisa_class_counts["Ligand"]+1
		m.UpdatePisaDicts(pisa_dicts, "Ligand")
		no_class = False
	if "DNA" in m.interface_type:
		pisa_class_counts["DNA"] = pisa_class_counts["DNA"]+1
		m.UpdatePisaDicts(pisa_dicts, "DNA")
		no_class = False
	if "RNA" in m.interface_type:
		pisa_class_counts["RNA"] = pisa_class_counts["RNA"]+1
		m.UpdatePisaDicts(pisa_dicts, "RNA")
		no_class = False
	#check for errors in the "Out" assignment
	if no_class:
		if "Out" in m.interface_type:
			pisa_class_counts["Out"] = pisa_class_counts["Out"]+1
			m.UpdatePisaDicts(pisa_dicts, "Out")
		else:
			sys.exit("ERROR: mutation not found in xml files was marked as found outside interfaces")
	#otherwise add the mutation's pdb code to the list
	else:
		if m.pdb not in pdbs_with_interface_mutations:
			pdbs_with_interface_mutations.append(m.pdb)

#download the structure files with at least one mutation in at least one interface
if check_and_download_files:
	for i in range(len(unique_pdbs)):
		if os.path.exists(path_pdb+unique_pdbs[i].lower()+".pdb"):
			print("Skipping structure file "+unique_pdbs[i].lower()+".pdb ( "+str(i+1)+" of "+str(len(unique_pdbs))+" )", end='\r')
		else:
			print("Retrieving structure file "+unique_pdbs[i].lower()+".pdb ( "+str(i+1)+" of "+str(len(unique_pdbs))+" )", end='\r')
			response=requests.get(structure_base_url+unique_pdbs[i].upper()+".pdb")
			if(response.status_code==200):
				outfile = open(path_pdb+unique_pdbs[i].lower()+".pdb",'wb')
				outfile.write(response.content)
				outfile.close()
	print("")

#analyze the structure files, looking, for each interface mutation, for the closest atom (and residue) in the other molecule
for i in range(len(pdbs_with_interface_mutations)):
	print("Analyzing contacts in structure file: "+pdbs_with_interface_mutations[i]+".pdb ( "+str(i+1)+" of "+str(len(pdbs_with_interface_mutations))+" )", end='\r')
	#build the list of (interface) mutations pointing to this pdb
	this_pdb_mutations = list()
	for m in covered_mutations:
		#skip mutations which do not appear in interfaces
		if m.get_interface_type_string() in ['-', '----']:
			continue
		if m.pdb == pdbs_with_interface_mutations[i]:
			this_pdb_mutations.append(m)
	#load the structure file
	parser = biopdb.PDBParser(QUIET = True)
	structure = parser.get_structure(pdbs_with_interface_mutations[i].lower(), path_pdb+pdbs_with_interface_mutations[i].lower()+".pdb")
	#create a dict associating each relevant chain's id to its Chain object
	chain_dict = dict()
	for parsed_chain in structure.get_chains():
		chain_dict[parsed_chain.get_id()] = parsed_chain
	#create a list of ligands
	ligands = list()
	for parsed_residue in structure.get_residues():
		if parsed_residue.id[0][0] == 'H': 
			ligands.append(parsed_residue)
	#for each mutation, find the closest residue and atom in the interacting chain
	for m in this_pdb_mutations:
		#get the position of the residue m occurs at
		pos = [0.0, 0.0, 0.0]
		for parsed_residue in chain_dict[m.chain]:
			if str(parsed_residue.get_id()[1]) == m.res_number:
				pos_sum = [0.0, 0.0, 0.0]
				pos_div = 0
				for parsed_atom in parsed_residue:
					c = parsed_atom.get_coord()
					for j in range(3):
						pos_sum[j] = pos_sum[j] + c[j]
					pos_div+=1
				for j in range(3):
					pos[j] = float(pos_sum[j]) / float(pos_div)
		#for each entry, scan the interacting chain and find the closest atom to pos
		for e in range(len(m.entries)):
			min_dist = None
			closest_atom = None
			closest_residue = None
			#if the interacting molecule is a ligand, instead of looking into a chain, ligands must be scanned
			if m.entries[e]['it'] == "Ligand":
				#scan m against the list of ligands
				for parsed_residue in ligands:
					for parsed_atom in parsed_residue:
						coords = parsed_atom.get_coord()
						dist = math.sqrt( ( pos[0]-coords[0])*(pos[0]-coords[0]) + (pos[1]-coords[1])*(pos[1]-coords[1]) + (pos[2]-coords[2])*(pos[2]-coords[2]) )
						#first atom of chain
						if min_dist is None:
							min_dist = dist
							closest_atom = parsed_atom.get_name()
							closest_residue = parsed_residue.get_resname()
						#else compare the distance from this atom to the current minimum distance
						else:
							if dist < min_dist:
								min_dist = dist
								closest_atom = parsed_atom.get_name()
								closest_residue = parsed_residue.get_resname()
				m.refine_interface_entry(e, closest_residue, closest_atom)
				continue
			#otherwise, check each residue in the interacting chain
			for parsed_residue in chain_dict[m.entries[e]['ic']]:
				#avoid waters and other heteroatoms
				if parsed_residue.id[0] in ['H', 'W']:
					continue
				#check each atom in the residue
				for parsed_atom in parsed_residue:
					coords = parsed_atom.get_coord()
					dist = math.sqrt( ( pos[0]-coords[0])*(pos[0]-coords[0]) + (pos[1]-coords[1])*(pos[1]-coords[1]) + (pos[2]-coords[2])*(pos[2]-coords[2]) )
					#first atom of chain
					if min_dist is None:
						min_dist = dist
						closest_atom = parsed_atom.get_name()
						closest_residue = parsed_residue.get_resname()
					#else compare the distance from this atom to the current minimum distance
					else:
						if dist < min_dist:
							min_dist = dist
							closest_atom = parsed_atom.get_name()
							closest_residue = parsed_residue.get_resname()
			#refine the entry with the new data just acquired
			m.refine_interface_entry(e, closest_residue, closest_atom)
print("")

#analyze the structure files with POPS, classifying non-interface mutations according to their location (surface or core)
pops_exclusion_list = list()
pops_skipped_pdbs = 0
pops_processed_pdbs = 0
for i in range(len(unique_pdbs)):
	print("Applying POPS to structure file: "+unique_pdbs[i]+".pdb ( "+str(i+1)+" of "+str(len(unique_pdbs))+" )", end='\r')
	#build the list of out-of-interface mutations pointing to this pdb
	this_pdb_mutations = list()
	for m in covered_mutations:
		#skip mutations which do appear in interfaces or have undefined class
		if m.get_interface_type_string() != '----':
			continue
		if m.pdb.lower() == unique_pdbs[i].lower():
			this_pdb_mutations.append(m)
	#skip files for which no mutations need to be collocated
	if len(this_pdb_mutations) == 0:
		pops_skipped_pdbs += 1
		continue
	#call POPS (Parameter OPtimized Surface of proteins and nucleic acids) on the given protein (skip if running in quick mode)
	if check_and_download_files:
		try:
			os.system("pops --pdb "+path_pdb+unique_pdbs[i]+".pdb --popsOut "+path_pops+unique_pdbs[i]+"_area.txt --multiModel --residueOut>temp.txt 2>&1")
		except:
			pops_exclusion_list.append(unique_pdbs[i])
			continue
	#parse POPS output, assigning a surface/core collocation to each mutation
	try:
		pops_file = open(path_pops+unique_pdbs[i]+"_area.txt","r")
		pops_text = pops_file.read()
		pops_file.close()
		pops_lines = pops_text.splitlines(pops_text.count("\n"))
		del pops_text
	except:
		pops_exclusion_list.append(unique_pdbs[i])
		continue
	#scan the output table
	for pl in pops_lines:
		cells = re.split("[\s]+",pl)
		#skip molecule lines and empty lines
		if len(cells) < 9:
			continue
		#skip headers
		if cells[0] == "ResidNe":
			continue
		#analyze current line
		res_name = cells[1].lower()
		res_chain = cells[2].rstrip()
		res_number = int(cells[3])
		res_qsasa = float(cells[8])
		#skip heteroatoms
		if res_name == "het":
			continue
		#check if any of the mutations on this PDB corresponds to this residue
		for m in this_pdb_mutations:
			if res_chain == m.chain and res_number == int(m.res_number):
				#if the residue has qsasa larger than threshold, assign surface collocation to this mutation
				if res_qsasa >= threshold_percent_sasa:
					m.is_on_surface = True
				#otherwise assign core collocation to this mutation
				else:
					m.is_on_surface = False
	pops_processed_pdbs += 1
print("")

#remove POPS temporary support files
if check_and_download_files:
	os.remove("popsb.out")
	os.remove("sigma.out")

#count surface and core mutations among non-interface mutations
collocation_surface_tot = 0
collocation_core_tot = 0
collocation_undefined_tot = 0
collocation_surface_R_X = 0
collocation_core_R_X = 0
collocation_undefined_R_X = 0
collocation_surface_G_X = 0
collocation_core_G_X = 0
collocation_undefined_G_X = 0
for m in covered_mutations:
	#skip mutations which do appear in interfaces or have undefined class
	if m.get_interface_type_string() != '----':
		continue
	#count undefined collocations
	if m.is_on_surface is None:
		collocation_undefined_tot += 1
		#check R/X
		if m.variation_protein[:3] == 'Arg':
			collocation_undefined_R_X += 1
		#check G/X
		if m.variation_protein[:3] == 'Gly':
			collocation_undefined_G_X += 1
	#count surface collocations
	elif m.is_on_surface:
		collocation_surface_tot += 1
		#check R/X
		if m.variation_protein[:3] == 'Arg':
			collocation_surface_R_X += 1
		#check G/X
		if m.variation_protein[:3] == 'Gly':
			collocation_surface_G_X += 1
	#count core colllocations
	else:
		collocation_core_tot += 1
		#check R/X
		if m.variation_protein[:3] == 'Arg':
			collocation_core_R_X += 1
		#check G/X
		if m.variation_protein[:3] == 'Gly':
			collocation_core_G_X += 1

#write report
print("Writing final report")
out_file = open(path_report, "w")
out_file.write("Number of mutations from ClinVar : "+str(len(mutations))+"\n")
out_file.write("Valid mutations analyzed : "+str(len(filtered_mutations))+"\n")
out_file.write("Of which R/X : "+str(total_mutations_from[aminoacid_lookup['R']])+"\n")
out_file.write("Of which G/X : "+str(total_mutations_from[aminoacid_lookup['G']])+"\n\n\n")
out_file.write("********************************************************************\n\n\n")
out_file.write("Mutations with a PDB file : "+str(len(raw_covered_mutations))+"\n")
out_file.write("Mutations with PDB residue mismatch : "+str(res_number_conflicts)+"\n\n")
out_file.write("Mutations with matching residue : "+str(len(covered_mutations))+"\n")
out_file.write("Of which R/X : "+str(len(R_X_mutations))+"\n")
out_file.write("Of which G/X : "+str(len(G_X_mutations))+"\n")
out_file.write("Unique PDB files : "+str(len(unique_pdbs))+"\n")
out_file.write("Of which empty : "+str(empty_xml_files)+"\n\n\n")
out_file.write("********************************************************************\n\n\n")
out_file.write("Mutations in Protein-DNA interfaces : "+str(pisa_class_counts["DNA"])+"\n")
out_file.write("Of which R/X : "+str(pisa_dicts["Arg"]["DNA"])+"\n")
out_file.write("Of which G/X : "+str(pisa_dicts["Gly"]["DNA"])+"\n\n\n")
out_file.write("********************************************************************\n\n\n")
out_file.write("Mutations in Protein-RNA interfaces : "+str(pisa_class_counts["RNA"])+"\n")
out_file.write("Of which R/X : "+str(pisa_dicts["Arg"]["RNA"])+"\n")
out_file.write("Of which G/X : "+str(pisa_dicts["Gly"]["RNA"])+"\n\n\n")
out_file.write("********************************************************************\n\n\n")
out_file.write("Mutations in Protein-Ligand interfaces : "+str(pisa_class_counts["Ligand"])+"\n")
out_file.write("Of which R/X : "+str(pisa_dicts["Arg"]["Ligand"])+"\n")
out_file.write("Of which G/X : "+str(pisa_dicts["Gly"]["Ligand"])+"\n\n\n")
out_file.write("********************************************************************\n\n\n")
out_file.write("Mutations in Protein-Protein interfaces : "+str(pisa_class_counts["Protein"])+"\n")
out_file.write("Of which R/X : "+str(pisa_dicts["Arg"]["Protein"])+"\n")
out_file.write("Of which G/X : "+str(pisa_dicts["Gly"]["Protein"])+"\n\n\n")
out_file.write("********************************************************************\n\n\n")
out_file.write("Mutations found only out of interfaces : "+str(pisa_class_counts["Out"])+"\n")
out_file.write("Of which R/X : "+str(pisa_dicts["Arg"]["Out"])+"\n")
out_file.write("Of which G/X : "+str(pisa_dicts["Gly"]["Out"])+"\n\n")
out_file.write("Surface collocation : "+str(collocation_surface_tot)+"\n")
out_file.write("Of which R/X : "+str(collocation_surface_R_X)+"\n")
out_file.write("Of which G/X : "+str(collocation_surface_G_X)+"\n\n")
out_file.write("Core collocation : "+str(collocation_core_tot)+"\n")
out_file.write("Of which R/X : "+str(collocation_core_R_X)+"\n")
out_file.write("Of which G/X : "+str(collocation_core_G_X)+"\n\n")
out_file.write("Undefined collocation : "+str(collocation_undefined_tot)+"\n")
out_file.write("Of which R/X : "+str(collocation_undefined_R_X)+"\n")
out_file.write("Of which G/X : "+str(collocation_undefined_G_X)+"\n\n")
out_file.write("Files without non-interface mutations : "+str(pops_skipped_pdbs)+"\n")
out_file.write("Files which could not be processed with POPS due to structural problems : "+str(len(pops_exclusion_list))+"\n")
out_file.write("Files processed with POPS : "+str(pops_processed_pdbs)+"\n\n\n")
out_file.write("********************************************************************\n\n\n")
out_file.write("Mutations not found in the xml file : "+str(pisa_class_counts["NotFound"])+"\n")
out_file.write("Of which R/X : "+str(pisa_dicts["Arg"]["NotFound"])+"\n")
out_file.write("Of which G/X : "+str(pisa_dicts["Gly"]["NotFound"])+"\n\n\n")
out_file.write("********************************************************************\n\n\n")
out_file.write("Legend for tables:\n")
out_file.write("N_Src = Base (source) nucleotide\n")
out_file.write("N_Dst = Mutant (destination) nucleotide\n")
out_file.write("A_Src = Base (source) aminoacid\n")
out_file.write("A_Dst = Mutant (destination) aminoacid\n")
out_file.write("PDB = Protein Data Bank structure identifier\n")
out_file.write("Chain = chain identifier in the PDB structure\n")
out_file.write("RN = Residue name in the PDB structure (might be different from A_src and A_dst)\n")
out_file.write("R# = Residue number in the PDB structure\n")
out_file.write("ITs = (InterfaceTypes) Types of interfaces this mutation was found in\n\n\n")
out_file.write("********************************************************************\n\n\n")
out_file.write("Table of R/X mutations\n\n")
out_file.write("N_Src\tN_Dst\tA_Src\tA_Dst\tPDB\t\tChain\tR#\t\tITs\n\n")
for m in R_X_mutations:
	if m.get_interface_type_string() in ['-', '----']:
		continue
	out_file.write(m.variation_nucleic[-3]+"\t\t")
	out_file.write(m.variation_nucleic[-1]+"\t\t")
	out_file.write(m.variation_protein[:3]+"\t\t")
	out_file.write(m.variation_protein[-3:]+"\t\t")
	out_file.write(m.pdb+"\t")
	out_file.write(m.chain+"\t\t")
	out_file.write(m.get_res_num_string())
	out_file.write(m.get_interface_type_string())
	out_file.write("\n")
out_file.write("\n\n********************************************************************\n\n\n")
out_file.write("Table of G/X mutations\n\n")
out_file.write("N_Src\tN_Dst\tA_Src\tA_Dst\tPDB\t\tChain\tR#\t\tITs\n\n")
for m in G_X_mutations:
	if m.get_interface_type_string() in ['-', '----']:
		continue
	out_file.write(m.variation_nucleic[-3]+"\t\t")
	out_file.write(m.variation_nucleic[-1]+"\t\t")
	out_file.write(m.variation_protein[:3]+"\t\t")
	out_file.write(m.variation_protein[-3:]+"\t\t")
	out_file.write(m.pdb+"\t")
	out_file.write(m.chain+"\t\t")
	out_file.write(m.get_res_num_string())
	out_file.write(m.get_interface_type_string())
	out_file.write("\n")
out_file.write("\n\n********************************************************************\n\n\n")
out_file.write("Detailed table of R/X mutations with list of interface entries in which each mutation occurs\n\n")
out_file.write("N_Src\tN_Dst\tA_Src\tA_Dst\tPDB\t\tChain\tR#\t\tITs\n\n")
for m in R_X_mutations:
	if m.get_interface_type_string() in ['-', '----']:
		continue
	out_file.write(m.variation_nucleic[-3]+"\t\t")
	out_file.write(m.variation_nucleic[-1]+"\t\t")
	out_file.write(m.variation_protein[:3]+"\t\t")
	out_file.write(m.variation_protein[-3:]+"\t\t")
	out_file.write(m.pdb+"\t")
	out_file.write(m.chain+"\t\t")
	out_file.write(m.get_res_num_string())
	out_file.write(m.get_interface_type_string())
	out_file.write("\n")
	for e in m.entries:
		out_file.write("Interface : type = "+str(e['it'])+" , mol = "+str(e['ic'])+" , res = "+str(e['ir'])+" , atom = "+str(e['ia'])+"\n")
	out_file.write("\n")
out_file.write("\n\n********************************************************************\n\n\n")
out_file.write("Detailed table of G/X mutations with list of interface entries in which each mutation occurs\n\n")
out_file.write("N_Src\tN_Dst\tA_Src\tA_Dst\tPDB\t\tChain\tR#\t\tITs\n\n")
for m in G_X_mutations:
	if m.get_interface_type_string() in ['-', '----']:
		continue
	out_file.write(m.variation_nucleic[-3]+"\t\t")
	out_file.write(m.variation_nucleic[-1]+"\t\t")
	out_file.write(m.variation_protein[:3]+"\t\t")
	out_file.write(m.variation_protein[-3:]+"\t\t")
	out_file.write(m.pdb+"\t")
	out_file.write(m.chain+"\t\t")
	out_file.write(m.get_res_num_string())
	out_file.write(m.get_interface_type_string())
	out_file.write("\n")
	for e in m.entries:
		out_file.write("Interface : type = "+str(e['it'])+" , mol = "+str(e['ic'])+" , res = "+str(e['ir'])+" , atom = "+str(e['ia'])+"\n")
	out_file.write("\n")
out_file.close()
print("Execution terminated succesfully")


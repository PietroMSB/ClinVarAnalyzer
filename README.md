# ClinVarAnalyzer
Script for the analysis of sets of mutations from the ClinVar [1] database:
https://www.ncbi.nlm.nih.gov/clinvar/

# Dependencies
The analyzer requires some libraries, programs, and python packages to run.

The script is written in Python:
https://www.python.org/
We used Python 3.6.10, but any Python>=3.5 installation should work.

The POPScomp [4] software can be downloaded from:
https://github.com/Fraternalilab/POPScomp

The following python packages can be installed with pip or Anaconda:
1) numpy : "pip install numpy" or "conda install numpy"
2) requests : "pip install requests" or "conda install requests"
3) matplotlib : "pip install matplotlib" or "conda install matplotlib"
4) xlsxwriter : "pip install XlsxWriter" or "conda install XlsxWriter"
5) biopython : "pip install biopython" or "conda install biopython"
6) lxml : "pip install lxml" or "conda install lxml"

# Input
To run the analyzer or the Blosum62 map builder, download your mutation files from ClinVar, or take one of our sets of mutations from the "InputFiles" directory. "BenignMissenseMutations" contains the list of all benign missense mutations downloaded from ClinVar on June 10, 2020. "PathogenicMissenseMutations" contains the list of all pathogenic missense mutations downloaded from ClinVar on June 10, 2020. The files were downloaded from ClinVar by selecting the "summary" format from the download menu.
The chosen input file (either taken from one of the two sub-directories or downloaded from ClinVar) should be named "clinvar_result.txt" and put in the "Data/" directory before starting. It is possible to change the required filename by changing the parameter "path_input" in both scripts.

# Mutation Maps
By mutation map we mean a matrix, obtained from an input set of ClinVar mutations, in which it is reported how many times aminoacid 1 (row) is substituted by aminoacid 2 (column).

# Blosum62 Mutation Maps
To build a map of the chosen set of ClinVar mutations, based on the Blosum62 matrix, simply run "python build_map_blosum62.py". It is possible to set a number of options by changing the parameters of the script. The most important ones are listed, just at the beginning, in the "execution modifiers" section:
1) normalize : if True, the number of mutations is normalized to a fixed value, to allow comparisons between sets of mutations of different cardinality.
2) normalization_multiplier : if "normalize" is True, sets the fixed value to which the map is normalized.
3) integer_normalization : if True, the normalized numbers of mutations are approximated to the closest integer, after being multiplied by "normalization_multiplier".
The map will be saved in the "Data/" directory, as "mutation_map.xlsx". It is possible to change the default name by changing the "path_output" parameter in the script.

# VarMap Queries
During the analysis, to associate each mutation to the PDB id of the structure it occurs on, the VarMap [2] web tool will be queried:
https://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/DisaStr/GetPage.pl?varmap=TRUE

# PDB Structure Files
The PDB structure files associated to the mutations with VarMap are downloaded from the Protein Data Bank [3]:
https://www.rcsb.org/
The structures are stored in the "Data/PISA/Structures/" sub.directory.

# Interface Files
During the analysis, interface files are downloaded from the PDBe-PISA (Proteins, Interfaces, Structures and Assemblies) website:
https://www.ebi.ac.uk/pdbe/pisa/
The interfaces of each PDB structure are listed in a .xml file, which is stored in the "Data/PISA/Interfaces/" sub-directory.

# POPS Output Files
To classify each mutation as "surface" or "core" mutation, the analyzer resorts to the POPS [4] tool. The output file of each POPS call is stored in the "Data/PISA/POPS_output/" sub-directory.

# ClinVar Analysis
To analyze a set of ClinVar mutations, once the input has been prepared, simply run "python analyze.py". The script will save, in the "Data/" directory, the following files:
1) dist_matrix.txt : distance matrix between aminoacids based on the minimum number of single nucleotide mutations needed to transform aminoacid 1 (row) into aminoacid 2 (column).
2) mutation_map.xlsx : map of the set of ClinVar mutations, based on the distance matrix, in .xlsx format.
3) mutation_matrix.txt : map of the set of ClinVar mutations in .txt format.
4) mutation_plot.png : color visualization of the distance matrix, corresponding to the color base of .xlsx map.
5) report.txt : the output file reporting all the results of the analysis.
6) VarMapInput.txtx : the list of mutation entries to submit to VarMap, in order to associate PDB structures to the mutations.
On a first run, the analyzer will stop after writing "VarMapInput.txt". To proceed, it is necessary to submit this file to the VarMap web tool. Once the results have been produced, they should be inserted in the "Data/" directory and renamed "VarMap_results.tsv" (the default name can be changed with the parameter "path_varmap_output" in the script). Then, the script should be run again as "python analyze.py".
When re-analyzing a set of mutations which was already succesfully analyzed, it is possible to speed up the execution by changing the parameter "check_and_download_files" to False. This requires all the structure files, interface files, and POPS output files to be in the respective sub-directories.

# References
[1] Landrum MJ, Lee JM, Riley GR, et al. ClinVar: public archive of relationships among sequence variation and human phenotype. Nucleic Acids Res. 2014;42(Database issue):D980-D985. doi:10.1093/nar/gkt1113

[2] Stephenson JD, Laskowski RA, Nightingale A, et al. VarMap: a web tool for mapping genomic coordinates to protein sequence and structure and retrieving protein structural annotations, Bioinformatics, Volume 35, Issue 22, 15 November 2019, Pages 4854–4856, https://doi.org/10.1093/bioinformatics/btz482

[3] Berman HM, Westbrook J, Feng Z, et al. The Protein Data Bank. Nucleic Acids Research, 28: 235-242, 2000. doi:10.1093/nar/28.1.235

[4] Cavallo L, Kleinjung J, Fraternali F. POPS: A fast algorithm for solvent accessible surface areas at atomic and residue level. Nucleic Acids Res. 2003;31(13):3364-3366. doi:10.1093/nar/gkg601

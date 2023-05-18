import sys
from Bio.PDB import *
import re
import os
import argparse
from subprocess import call
### File loader ###
description=\
"This script takes a PDB and outputs a tabulated list file of 1) Protein annotation 2) region annotation e.g. MHCa1  \
3) residue number and 4) residue one-letter code \
The script takes the argument chains which flags the MHCpTCRab in the order MHCA,MHCB,peptide,TCRA,TCRB"

### File loader ###

def parse_args():
    """parse arguments for command line"""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--PDB', dest = 'pdb', type = str, required = True, help ='The pdb file to be analysed')
    parser.add_argument('--name', dest = 'name', type = str, required = True, help ='The name of the structure to be analysed')
    parser.add_argument('--MHCclass', dest = 'MHCclass', type = str, required = True, help ='The MHC class of the complex structure: I or II')
    parser.add_argument('--chains', dest = 'chains',   type = str, required = True, help ='Chains of TCR-pMHC complex in order MHCa,MHCb,peptide,TCRa,TCRb')
    parser.add_argument('--output_dir', dest = 'output_dir',   type = str, required = True, help ='The output directory of your analysis.')


    args = parser.parse_args()
    return args

def readFile(FILE, fileType):
    if FILE == None:
        sys.exit('No file to read')
    if FILE.split('.')[-1].lower() != str(fileType):
        sys.exit("\nFile extension must be of type ."+str(fileType)+"\n")
    else:
        print('Reading file: ' + str(FILE))
        return open(FILE, "r")

### Functions ###

def three2one(seq):
    
    d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
    
    if len(seq) %3 == 0:
        upper_seq= seq.upper()
        single_seq=''
        for i in range(len(upper_seq)//3):
            single_seq += d[upper_seq[3*i:3*i+3]]
        return single_seq
	
def lineFormat(line):
    line2=[]
    line2= line.split(" ")	
    line3=line2[1:]
    line4=[] 
    line4.append(line3[2])
    line4.append(line3[3])
    line4.append(line3[1])
    line4.append(three2one(line3[0]))
    return line4
    
def PDB_to_list(PDB_resi, chain_name, chain, MHCclass):
    outTxt=''
    for i in range (0, len(PDB_resi)):
        if is_aa(PDB_resi[i]):
            line = str(PDB_resi[i])
            line = line[8:]
            line = line.replace(' het=  resseq=', ' ')
            if line[-2].isalnum() == False:
                line = line.replace(' icode= >', ' ')
                line = line + chain_name + ' '    
                if chain == 'MHCA' and MHCclass=="I":
                    number = int(re.sub("[^0-9]", "", line))
                    if 50 <= number <= 86:
                        line = line + 'MHCa1'
                    if 140 <= number <= 176:
                        line = line + 'MHCa2'
                if chain == 'MHCA' and MHCclass =='II':            
                    number = int(re.sub("[^0-9]", "", line))
                    if 46 <= number <= 78:
                        line = line + 'MHCa1'            
                if chain == 'MHCB' and MHCclass =='II':            
                    number = int(re.sub("[^0-9]", "", line))
                    if 54 <= number <= 91:        
                        line = line + 'MHCb1'
                lineF=lineFormat(line)
            
                for x in lineF:
                    outTxt+=str(x)+"\t"
                outTxt+="\n"
            else:
                print("Ommitted icode residue!!")
    return outTxt

###### BODY #####
print('     ~  Running PDB_sequence_panClass.py v0.1 BETA  ~')
# Arguments #
args = parse_args()
pdb = args.pdb
structure_name = args.name
MHCclass = args.MHCclass[0:2]
chains = args.chains
output_dir = args.output_dir

origPDB = readFile(pdb, "pdb")
fileName=pdb.rsplit('/')[-1]
PDB_name = fileName.rsplit('.')[0]


print("structure name is ", structure_name)

if not os.path.exists(output_dir):
    print("Creating Directory "+output_dir)
    os.makedirs(output_dir)
    
if not os.path.exists(output_dir+"/sequences"):
    print("Creating Directory "+output_dir+"/sequences")
    os.makedirs(output_dir+"/sequences")
    
filteredFile = readFile(pdb, "pdb")

# Sort chains
MHCa, MHCb, pep = chains[0], chains[1], chains[2]

#create a parser object
parser = PDBParser()

#N.B. PDBParser is set to allow for buggy PDB structures by a defaut PERMISSIVE=1
#can change this to =0 to throw an error outright. I have kept it to allow errors
#as they often don't impact anything.

struct   = parser.get_structure(filteredFile, pdb)

print("The warnings re: discontinous sequences are fine, ignore.")

#subset chains

MHCa_chain = struct[0][MHCa]
MHCb_chain = struct[0][MHCb]
pep_chain  = struct[0][pep]


#get residues
MHCa_residues   = list(Selection.unfold_entities(MHCa_chain, 'R'))
MHCb_residues   = list(Selection.unfold_entities(MHCb_chain, 'R'));
pep_residues    = list(Selection.unfold_entities(pep_chain, 'R'))


#regex out the superflous stuff
print(type(MHCa_residues))
print(MHCa_residues[0])

outTxt=''
outTxt+="Chain \tAnnotation \tResNum \tResCode\n"

MHCaTxt=''
MHCbTxt=''
peptideTxt=''
TCRaTxt=''
TCRbTxt=''


MHCaTxt    = PDB_to_list(MHCa_residues, 'MHCA', 'MHCA', MHCclass)
MHCbTxt    = PDB_to_list(MHCb_residues, 'MHCB', 'MHCB', MHCclass)
peptideTxt = PDB_to_list(pep_residues, 'peptide', 'none',MHCclass)


outTxt+=MHCaTxt
outTxt+=MHCbTxt
outTxt+=peptideTxt





print(outTxt)
outfile = open(output_dir+"/sequences/"+structure_name+'_sequence.txt', 'w')
outfile.write(outTxt)
outfile.close()

print('     ~  End of PDB_sequence_panClass.py v0.1 BETA  ~')



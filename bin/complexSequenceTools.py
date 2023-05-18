from subprocess import call
import os
import sys
import argparse
import time

description=\
"This script takes a PDB file of a TCR-pMHC complex structure and creates two files:\n 1) A fasta file containing information about the \
TCRs gene usage and CDR regions. 2) a tabulated sequence list containing residue, chain and CDR region annotations. It also churns out some intermediate files.\
It uses multiple scripts to do so which must be located in the same directory. These files are created in the directory name/sequences"

### File loader ###

def parse_args():
    """parse arguments for command line"""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--PDB', dest = 'pdb', type = str, required = True, help ='The PDB file to be analysed')
    parser.add_argument('--MHCclass', dest = 'MHCclass', type = str, required = True, help ='The MHC class of the complex structure: I or II')
    parser.add_argument('--Chains', dest = 'chains',   type = str, required = False, help ='Chains of TCR-pMHC complex in order MHCa,MHCb,peptide,TCRa,TCRb', default = "ABCDE")
    
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


##################################################### BODY ##################################################

print('     ~  Running complexSequenceTools.py v1.0 BETA  ~')

# Arguments #
args = parse_args()
pdb = args.pdb
chains = args.chains
MHCclass=args.MHCclass

origPDB = readFile(pdb, "pdb")
fileName=pdb.rsplit('.', 1)[0]

if not os.path.exists(fileName):
    print("Creating Directory "+fileName)
    os.makedirs(fileName)
    
if not os.path.exists(fileName+"/sequences"):
    print("Creating Directory "+fileName+"/sequences")
    os.makedirs(fileName+"/sequences")
    
if not os.path.exists(fileName+"/pdbs"):
    print("Creating Directory "+fileName+"/pdbs")
    os.makedirs(fileName+"/pdbs")
        
########### Make fasta from PDB ############
    
call(["python", "bin/clean_PDB.py", "--PDB", fileName+".pdb"], shell=True)
#Outputted file: 1fyt/pdbs/1fyt_filtered.pdb
    
call(["python", "bin/pdb_to_fasta_v0.3.py", "--PDB", fileName+".pdb", "--Chains", chains], shell=True)
#Outputted file: 1fyt/sequences/1fyt.fasta & 1fyt/sequences/1fyt_TCR.fasta
#time.sleep(5)
########## Find TR gene usage ##############
call(["python", "bin/TCRgeneFinderICBINBlast.py", "--Fasta", fileName+"/sequences/"+fileName+".fasta"], shell=True)
#Outputted file: 1fyt/sequences/1fyt_TCRannot.fasta
#time.sleep(5)
########## Find CDRs #######################
call(["python", "bin/TCRannotate_v0.8.py", "--Fasta", fileName+"/sequences/"+fileName+"_TCRannot.fasta"], shell=True)
#Outputted file: 1fyt/sequences/1fyt_TCRannot_CDRannot.fasta
#time.sleep(5)
###### Make full sequence list from PDB ####
call(["python", "bin/PDB_sequence_panClass.py", "--PDB", fileName+".pdb", "--MHCclass", MHCclass], shell=True)
#Outputted file: 1fyt/sequences/1fyt_sequence.txt
#time.sleep(5)
#### Add CDR annotations to sequence list ##
call(["python", "bin/annotateSequenceFile_v0_2.py", "--SeqTXT", fileName+"/sequences/"+fileName+"_sequence.txt", "--Fasta", fileName+"/sequences/"+fileName+"_TCRannot_CDRannot.fasta"], shell=True)
#Outputted file: 1fyt/sequences/1fyt_sequence_annot.txt

print('     ~  End complexSequenceTools.py v1.0 BETA  ~')


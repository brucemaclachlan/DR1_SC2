#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written using Sublime Text & Spyder

@author:    Bruce J MacLachlan, Division of Infection & Immunity, Cardiff University
            ORCID ID: 0000-0002-2685-2733
            Github: https://github.com/brucemaclachlan
            email: maclachlanb@cardiff.ac.uk
            
            T
"""
import subprocess
import os
import sys
import argparse
import time

description=\
"foobar"

### File loader ###

def parse_args():
    """parse arguments for command line"""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--PDB', dest = 'pdb', type = str, required = True, help ='The PDB file to be analysed')
    parser.add_argument('--MHCclass', dest = 'MHCclass', type = str, required = True, help ='The MHC class of the complex structure: I or II')
    parser.add_argument('--chains', dest = 'chains',   type = str, required = True, help ='Chains of TCR-pMHC complex in order MHCa,MHCb,peptide')
    parser.add_argument('--peptide_colour', dest = 'colour',   type = str, required = False, help ='Hex colour (with #). Colour of the peptide for graphing purposes', default = "#808080")
    parser.add_argument('--output_dir', dest = 'output_dir',   type = str, required = True, help ='The name of the output directory.')


    #parser.add_argument('--SeqAnalysisForce', dest = 'SeqAnalysisForce', type = bool, required = False, help ='Force rerun of complexSequenceTools.py even if sequence information exists. True/False', default = False)
    
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

print('     ~  Running projectContactsFull.py v1.0 BETA  ~')
# Arguments #
args = parse_args()
pdb = args.pdb
chains = args.chains
MHCclass=args.MHCclass
peptide_colour = args.colour
output_dir = args.output_dir
print(peptide_colour)
#SeqAnalysisForce=args.SeqAnalysisForce

print("Running analysis with the following inputs.. ")
print(args)

origPDB = readFile(pdb, "pdb")
fileName=pdb.rsplit('.', 1)[0]

if not os.path.exists(output_dir ):
    print("Creating Directory "+output_dir )
    os.makedirs(output_dir )
    
if not os.path.exists(output_dir +"/contacts"):
    print("Creating Directory "+output_dir +"/contacts")
    os.makedirs(output_dir +"/contacts")

contactPath=str(output_dir+"/contacts/")

pdb_real = os.path.abspath(pdb)

chains_dict = { "MHCA" : chains[0],
                "MHCB" : chains[1],
                "PEPTIDE" : chains[2]
}


#### STRIP ALL HETATMS #####

if not os.path.exists(output_dir +"/pdbs"):
    print("Creating Directory pdbs")
    os.makedirs(output_dir +"/pdbs")

import pymol
pymol.cmd.load(pdb_real, "structure")
pymol.cmd.remove("hetatm")
pymol.cmd.save(output_dir +"/pdbs/"+fileName+"_clean.pdb" ,selection='(all)' , format='pdb')

pdb_clean = output_dir +"/pdbs/"+fileName+"_clean.pdb"
pdb_clean_real = os.path.abspath(pdb_clean)

########################## Check whether we need to run complexSequenceTools.py ############################

### The location of complexSequenceTools.py will need to be changed once we have a /bin"
#print(SeqAnalysisForce)
#if SeqAnalysisForce == True:
#    print("User flagged to rerun sequence analysis tools. Running complexSequenceTools.py.")
#    print("Caught here 1")
#    time.sleep(20)
#    call(["python", "complexSequenceTools.py", "--PDB", pdb, "--MHCclass", MHCclass, "--Chains", chains], shell=True)
#
#    
#
#if not os.path.exists(fileName+"/sequences/"+fileName+"_sequence_annot.txt"):
#    print("Could not find sequence information required. Running complexSequenceTools.py.")
#    print("Caught here 2")
#    time.sleep(20)
#    call(["python", "complexSequenceTools.py", "--PDB", pdb, "--MHCclass", MHCclass, "--Chains", chains], shell=True)  
#    print("Caught here 2")
#    

#################################################### p to MHC ##################################################


if not os.path.exists(contactPath+"pMHC_contacts"):
    print("Creating Directory "+contactPath+"pMHC_contacts")
    os.makedirs(contactPath+"pMHC_contacts")
MHCpath=contactPath+"pMHC_contacts/"
## Run NCONT ##    
pMHC_ncont_parameters=open(MHCpath+fileName+"_ncont_parameters.txt", "w")
pMHC_ncont_parametersT="\
source \""+chains_dict["MHCA"]+","+chains_dict["MHCB"]+"/*\" \n\
target \""+chains_dict["PEPTIDE"]+"/**\"\n\
mindist 0.0\n\
maxdist 4.0\n\
cells OFF\n\
END"
pMHC_ncont_parameters.write(pMHC_ncont_parametersT)
pMHC_ncont_parameters.close()
pMHC_ncont_out=open(MHCpath+fileName+"_ncont_out.txt","w")
 

# problem is passing in the parameters file

# print(pdb_real)
ncont_file = os.path.abspath(MHCpath+fileName+"_ncont_parameters.txt")
# print(ncont_file)

ncont_txt_command = " ".join(["ncont", "XYZIN", pdb_clean_real,"<", ncont_file, ">", MHCpath+fileName+"_ncont_out.txt"])
print(ncont_txt_command)

subprocess.call([ncont_txt_command], shell=True)
pMHC_ncont_out=open(MHCpath+fileName+"_ncont_out.txt","r")

for line in iter(pMHC_ncont_out): print(line)
 
## Annotate Contacts ##
proj_contacts_command = " ".join(["python", "bin/projectContacts_v0_6.py", "--NCONT", MHCpath+fileName+"_ncont_out.txt", "--chains" , chains])
subprocess.call([proj_contacts_command], shell=True)
#Outputted file: 1fyt/TCRpMHC_contacts/1fyt_ncont_out_contacts.txt

## Create per residue data ##

residue_contacts_command = " ".join(["python", "bin/residueContacts_v0.3.py", "--ContactTXT", MHCpath+fileName+"_ncont_out_contacts.txt", "--chains" , chains])
subprocess.call([residue_contacts_command], shell=True)

#Outputted file: 1fyt/TCRpMHC_contacts/1fyt_ncont_out_contacts_residues.txt


###### Make full sequence list from PDB ####

residue_contacts_command = " ".join(["python", "bin/PDB_sequence_panClass.py", "--PDB", pdb_clean_real, "--name" , fileName , "--MHCclass", MHCclass, "--output_dir", output_dir, "--chains" , chains])
subprocess.call([residue_contacts_command], shell=True)
#Outputted file: 1fyt/sequences/1fyt_sequence.txt



## Add per residue data to full PDB sequence ##

full_residue_command = " ".join(["python", "bin/fullResidueContacts_v0.1.py", "--SeqTXT", output_dir+"/sequences/"+fileName+"_sequence.txt", "--ContactTXT", MHCpath+fileName+"_ncont_out_contacts_residues.txt", "--chains", chains])
subprocess.call([full_residue_command], shell=True)



## Create contact maps
contact_map_command = " ".join(["python", "bin/contactMap_pMHC.py", "--Input", MHCpath+fileName+"_ncont_out_contacts_residues_contacts_residues_full.txt", "--MHCclass", MHCclass])
subprocess.call([contact_map_command], shell=True)


print(peptide_colour)
print(type(peptide_colour))
## Create contact plots
contact_plots_command = " ".join(["python", "bin/contacts_graphs.py", "--SeqTXT", output_dir+"/sequences/"+fileName+"_sequence.txt", "--ContactTXT", MHCpath+fileName+"_ncont_out_contacts.txt", '--peptide_colour', peptide_colour[1:], "--output_dir", output_dir])
subprocess.call([contact_plots_command], shell=True)

print('     ~  End projectContactsFull.py v1.0 BETA  ~')

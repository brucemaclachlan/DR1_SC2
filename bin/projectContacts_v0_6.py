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

import os
import sys
import argparse

description=\
"This script takes a CCP4 NCONT contact file and annotates the contact data due to a set of criteria. A tabulated txt file is outputted containing the annotated contact data.\
The NCONT data are filtered to contacts less than 4.0 Angstroms and annotated as van der Waals, hydrogen bond or salt bridge interactions. All interactions involving H2O are ommitted."

### File loader ###

def parse_args():
    """parse arguments for command line"""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--NCONT', dest = 'ncont', type = str, required = True, help ='The ncont output file to be analysed')
    parser.add_argument('--chains', dest = 'chains',   type = str, required = True, help ='Chains of TCR-pMHC complex in order MHCa,MHCb,peptide,TCRa,TCRb')
    parser.add_argument('--Fasta', dest = 'fasta', type = str, required = False, help = 'A fasta file containing annotated TCR CDR loop locations.', default = "")

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

# Parsers and ting #
    
def everythingParser(inFile):
    everythingList=[]
    for lines in inFile.readlines():
        everythingList.append(lines[:-1])
    
    print('\n'+ "Input file contains " + str(len(everythingList)) + ' lines')
    return everythingList
    
def ccp4ContactValidator(inFile):
    inFileYo=inFile
    for lines in inFileYo:
        if "NCONT" in lines:
            print("Input file was validated as a ccp4 contact output file!")
            return True
    else:
        sys.exit("\nThis doesn't look like an NCONT output file. Please try again!")
        
def contactParser(inFile):
	inFileYo=inFile
	contactList=[]
	for lines in inFileYo:
		if "]:" in lines:
			contactList.append(lines)
	print(str(len(contactList)) + " contact lines detected")
	return contactList
    
def doesContactsMatch(everythingList, contactList):
    everythingListYo=everythingList
    contactListYo=contactList
    totalContactLine=''
    reference=0
    analyte=len(contactListYo)
    for lines in everythingListYo:
        if "Total" in lines and "contacts" in lines:
            totalContactLine += lines
    for s in totalContactLine.split(): 
        if s.isdigit():
            reference+=int(s)   
    if analyte == reference:
        print("projectContact recognised the same number of contacts as reported by CCP4!")
    else:
        sys.exit("projectContacts did not find all the contacts in the contact file! Please report this error")
        
def fastaParser(fasta):
    from Bio import SeqIO
    entries=[]
    for seq_record in SeqIO.parse(fasta, "fasta"):
        ID=''
        sequence=''
        entry=[]
        ID= "".join(seq_record.id)
        sequence = "".join(seq_record.seq)
        entry.append(ID)
        entry.append(sequence)
        entries.append(entry)   
    return entries

def depackID(entries):
    newEntries=[]
    for entry in entries:
        ID=entry[0]
        sequence=entry[1]
        newEntry=[]
        id_terms=ID.split("|")
        #print(ID.split("|"))
        for term in id_terms:
            if len(term)!=0:
                newEntry.append(term)
        newEntry.append(sequence)
        
        newEntries.append(newEntry)

    return newEntries

def findLocations(subentries):
    s = "=[]"
    locations = []
    for col in subentries:
        if s[0] in col and s[1] in col and s[2] in col:
            locations.append(col)
    return locations

def depackLocations(subentries):
    locations = []
    for col in subentries:
        location=[]
        location.append(col.rsplit("=", 1)[0])
        locationstring = (col.partition('[')[-1].rpartition(']')[0])
        location += map( int, locationstring.split(',') )
        locations.append(location)
    return locations

def purgeCysLocs(locations):
    output=[]
    for locs in locations:
        if "Cys" not in locs[0]:
            output.append(locs)
    return output
        
#####           MATHS           #######

## Line Filler ##

test=' /1/A/  65(ARG). / CG [ C]:  /1/E/  57(SER). / CB [ C]:   3.98'

def lineNeedsFilled(line):
    if line.count(":") == 1:
        return True
    else:
        return False
        
def sourceAtomsGrabber(line):
    splitLine=line.split(':')
    return str(splitLine[0]+':')
    
def lineFiller(contactLines):
    print('\nFilling source atom lines...')
    newLines=[]
    positionCounter=0
    fillCounter=0
    for lines in contactLines:
        positionCounter+=1
        if lineNeedsFilled(lines) == False:
            newLines.append(lines)
    
        elif lineNeedsFilled(lines) == True:
            fillCounter+=1
            trackBackCounter=2
            newLine=''
            suffix=''
            prefix=''
            while str(sourceAtomsGrabber(contactLines[positionCounter-trackBackCounter])).count(' ') > 20:
                trackBackCounter+=1
            prefix=sourceAtomsGrabber(contactLines[positionCounter-trackBackCounter])
            suffix=lines[28:]
            newLine = prefix + ' ' + suffix
            newLines.append(newLine)
    print(str(fillCounter) + " lines were filled with source atoms...")
    return newLines  
    
def contactMatrixRow(line, TCRAlocations, TCRBlocations):
    
    # newNames = {MHCachain: "MHCa",
    #         MHCbchain: "MHCb",
    #         peptidechain: "peptide",
    #         TCRachain: "TCRa",
    #         TCRbchain: "TCRb",
    #     }
    SnewName=newNames[line[4]]
    TnewName=newNames[line[32]]
    
    annotation1=''
    if SnewName == "TCRa":
        for loop in TCRAlocations:
            for loc in loop[1:]:
                if loc == int(line[6:10]):
                    annotation1 = loop[0]
                    
    if SnewName == "TCRb":
        for loop in TCRBlocations:
            for loc in loop[1:]:
                if loc == int(line[6:10]):
                    annotation1 = loop[0]   
    annotation2=''
    contactRow=[]
    contactRow.append(line[4])
    contactRow.append(SnewName)
    contactRow.append(int(line[6:10]))
    contactRow.append(annotation1)
    contactRow.append(line[11:14])
    contactRow.append(line[19:22])
    contactRow.append(line[24])
    contactRow.append(line[32])
    contactRow.append(TnewName)
    contactRow.append(int(line[34:38]))
    contactRow.append(annotation2)
    contactRow.append(line[39:42])
    contactRow.append(line[47:50])
    contactRow.append(line[52])
    contactRow.append(float(line[-4:]))
    
    return contactRow
 
def contactMatrix(contactLines, TCRAlocations, TCRBlocations):
    omitCounter = 0
    print('\nCreating contact matrix...')
    contactMatrix=[]
    for lines in contactLines:
        if 'HOH' in lines:
            omitCounter +=1 
        else:
            contactMatrix.append(contactMatrixRow(lines, TCRAlocations, TCRBlocations))
    print(str(omitCounter) + ' contacts were omitted due to water "HOH"...' )
    return sorted(contactMatrix, key=lambda items: (items[0], items[2]))
                
def isHBond(contactRow):

    hbDonors=["ARGNE N","ARGNH1N","ARGNH2N","ASNND2N","GLNNE2N","HISND1N","HISNE2N","LYSNZ N",\
    "SEROG O","THROG1O","TRPNE1N","TYROH O","ALAN  N","ARGN  N","ASNN  N","ASPN  N","CYSN  N",\
    "GLUN  N","GLNN  N","GLYN  N","HISN  N","ILEN  N","LEUN  N","LYSN  N","METN  N","PHEN  N",\
    "SERN  N","THRN  N","TRPN  N","TYRN  N","VALN  N"]
    
    hbAcceptors=["ASNOD1O","ASPOD1O","ASPOD2O","GLNOE1O","GLUOE1O","GLUOE2O","HISND1O","HISND2O",\
    "SEROG O","THROG1O","TYROH O","ALAO  O","ARGO  O","ASNO  O","ASPO  O","CYSO  O","GLUO  O",\
    "GLNO  O","GLYO  O","HISO  O","ILEO  O","LEUO  O","LYSO  O","METO  O","PROO  O", "PHEO  O","SERO  O",\
    "THRO  O","TRPO  O","TYRO  O","VALO  O"]

    ID1=contactRow[4]+contactRow[5]+contactRow[6]
    ID2=contactRow[11]+contactRow[12]+contactRow[13]
    if ID1 in hbDonors and ID2 in hbAcceptors:     
        if contactRow[14] <= 3.40:
            return True
    if ID2 in hbDonors and ID1 in hbAcceptors:    
        if contactRow[14] <= 3.40:
            return True
    else:
        return False  

def isSaltBridge(contactRow):
    acidAtoms=["GLUOE1O","GLUOE2O","ASPOD1O","ASPOD2O"]
    baseAtoms=["LYSNZ N","ARGNE N","ARGNH1N","ARGNH2N"]

    ID1=contactRow[4]+contactRow[5]+contactRow[6]
    ID2=contactRow[11]+contactRow[12]+contactRow[13]
    if ID1 in acidAtoms and ID2 in baseAtoms:     
        if contactRow[14] <= 3.40:
            return True
    if ID2 in acidAtoms and ID1 in baseAtoms:    
        if contactRow[14] <= 3.40:
            return True
    else:
        return False  
                      
    
def isVdW(contactRow):
    if contactRow[4] == contactRow[11]:
        if contactRow[5] == contactRow[12]:
            if contactRow[7] == contactRow[13]:
                #print("same atom hit")
                return False
    if contactRow[14] <= 4.0:
        return True
        
def bondAnnotator(contactsRow):
    contactsRow.append('NO') 

    if isVdW(contactsRow)==True:
        contactsRow[15] ='VW'
    if isHBond(contactsRow)==True:
        contactsRow[15] ='HB'
    if isSaltBridge(contactsRow)==True:
        contactsRow[15] ='SB'
        
def annotateAllWrapper(contactMatrix):
    contactMatrixNew =[]
    omitCounter = 0
    print('Anotating contacts...')
    for row in contactMatrix:
        bondAnnotator(row)
        if row[15] == "NO":
            omitCounter +=1
        else:
            contactMatrixNew.append(row)
            
    return contactMatrixNew
    
    


    
    print(str(omitCounter) + ' contacts were omitted due to not meeting annotation criteria "NO"...' )


    
#####       INITIALISER        #######
# Openers #
print('\n''     ~  Running projectContacts.py v0.6  ~')

args = parse_args()
ncont = args.ncont
chains = args.chains
fasta = args.fasta

inFile = readFile(ncont, "txt")

inFileName=ncont.rsplit('.', 1)[0]

if type(inFileName) != str:
    sys.exit('No file was loaded. Please view usage and provide a valid file to be processed')
allLines=everythingParser(inFile)
contactLines=contactParser(allLines)


chains_dict = { "MHCA" : chains[0],
                "MHCB" : chains[1],
                "PEPTIDE" : chains[2]
}

# Validation #
if ccp4ContactValidator(allLines)==False:
    raise IOError("Input file was NOT validated as a ccp4 contact output file! /n Please see usage for information")
if doesContactsMatch(allLines, contactLines) == False:
    raise IOError("projectContacts did not find all the contacts in the contact file!")
    
# Retrieve fasta info

if fasta != "":
    FASTAfile = readFile(fasta, "fasta")
    fastaEntries=fastaParser(fasta)
    fastaEntries=depackID(fastaEntries)
    
    TCRA = []
    for entry in fastaEntries:
        if "TCRA" in entry:
            TCRA = entry
    TCRAlocations = findLocations(TCRA)
    TCRAlocations = depackLocations(TCRAlocations)   
    TCRAlocations = purgeCysLocs(TCRAlocations)
    
    TCRB = []
    for entry in fastaEntries:
        if "TCRB" in entry:
            TCRB = entry
    TCRBlocations = findLocations(TCRB)
    TCRBlocations = depackLocations(TCRBlocations)
    TCRBlocations = purgeCysLocs(TCRBlocations)
    
else:
    TCRAlocations = [""]
    TCRBlocations = [""]

# Sort chains

if len(chains) == 3:
    MHCachain, MHCbchain, peptidechain = chains[0], chains[1], chains[2]

    newNames = {MHCachain: "MHCa",
            MHCbchain: "MHCb",
            peptidechain: "peptide",
            "D": "TCRa",
            "E": "TCRb",
            }

if len(chains) == 5:    
    MHCachain, MHCbchain, peptidechain, TCRachain, TCRbchain = chains[0], chains[1], chains[2], chains[3], chains[4]

    newNames = {MHCachain: "MHCa",
                MHCbchain: "MHCb",
                peptidechain: "peptide",
                TCRachain: "TCRa",
                TCRbchain: "TCRb",
            }
          
                            

#####           BODY           #######
contactMatrixS=[]
contactLines=lineFiller(contactLines)
contactMatrix=contactMatrix(contactLines, TCRAlocations, TCRBlocations)
contactMatrix = annotateAllWrapper(contactMatrix)

for x in contactMatrix:
    print(x)
    
# ## Temp graph generator

# from collections import Counter
# import numpy as np
# import matplotlib as mpl
# import matplotlib.pyplot as plt
# import matplotlib.ticker as ticker

# peptide_contacts = []
# for line in contactMatrix:
#     peptide_contacts.append(line[9])

# print(peptide_contacts)
# peptide_contacts_counted = Counter(peptide_contacts)

# print(peptide_contacts_counted)

# positions = []
# counts = []
# for position in peptide_contacts_counted:
#     positions.append(position)
#     counts.append(peptide_contacts_counted[position])

# total_contacts = float(sum(counts))
# perc_contacts = []
# for p in counts:
#     perc_contacts.append(float(float(p)/total_contacts*100))
    

# a,b = np.asarray(positions),np.asarray(perc_contacts)
# data = np.vstack((a,b)).T

# print(data)
# print(data.shape)


# nrow = 1
# ncol = 1
# fig, axes = plt.subplots(nrow, ncol, figsize=(3.75,1.25), sharex=False, sharey=False, facecolor="None")

# bar_colours = ["grey","grey","grey","k","k","k","k","k","k","k","k","k","grey","grey"]
# print(bar_colours)

# axes.bar(data[:,0], data[:,1], align="center", color=bar_colours)
# axes.set_xlim(0,21)
# axes.set_ylim(0,15)
# axes.set_facecolor("None")
# axes.spines['left'].set_linewidth(1)
# axes.spines['bottom'].set_linewidth(1)
# axes.spines['right'].set_visible(False)
# axes.spines['top'].set_visible(False)

# labels = ["","F","A","R","R","P","P","L","A","E","L","A","A","L","N","L","S","G","S","R","L", ""]
# axes.set_xticks(np.arange(21))
# axes.set_xticklabels(labels)
# labels = axes.get_xticklabels()
# plt.setp(labels, rotation=0, fontsize=12, fontname="Arial", fontweight="bold")

# axes.set_yticks(np.arange(0.0, 15.1, 5.0))
# labels = axes.get_yticklabels()
# plt.setp(labels, rotation=0, fontsize=12, fontname="Arial", fontweight="bold")

# axes.set_ylabel("% total contacts", fontsize=10, fontweight="bold", fontname="Arial")
# axes.set_title("p to MHC contacts", fontsize=12, fontweight="bold", fontname="Arial", y=1.1)

# print(data[1,1])

# axes.text(4, data[0,1]+1.4, "P-3", fontsize=10, rotation=0, fontweight="bold", fontname="Arial",va="center", ha="center", color="grey")
# axes.text(7, data[3,1]+1.4, "P1", fontsize=10, rotation=0, fontweight="bold", fontname="Arial",va="center", ha="center", color="k")
# axes.text(10, data[6,1]+1.4, "P4", fontsize=10, rotation=0, fontweight="bold", fontname="Arial",va="center", ha="center", color="k")
# axes.text(12, data[8,1]+1.4, "P6", fontsize=10, rotation=0, fontweight="bold", fontname="Arial",va="center", ha="center", color="k")
# axes.text(15, data[11,1]+1.4, "P9", fontsize=10, rotation=0, fontweight="bold", fontname="Arial",va="center", ha="center", color="k")



# fig.savefig("perc_contacts.svg", bbox_inches='tight', transparent=True)

#####   OUTPUT GENERATOR       #######   
            
outFile = open(str(inFileName) + '_contacts.txt' , 'w')
output2=''
output2+= 'Donor_Chain_Letter \tDonor_Chain \tDonor_ResNum \tDonor_Annotation \tDonor_ResCode \tDonor_Position \tDonor_Atom \tAcceptor_Chain_Letter \tAcceptor_Chain \tAcceptor_ResNum \tAcceptor_Annotation \tAcceptor_ResCode \tAcceptor_Position \tAcceptor_Atom \tDistance \tType\n'

for x in contactMatrix:
    for y  in x:
        output2+=str(y) + '\t'
    output2+='\n'

outFile.write(output2)
print('\nOutputted file: ' + str(inFileName)+ '_contacts.txt')
print('\n''     ~  End ProjectContacts.py v0.6      ~')

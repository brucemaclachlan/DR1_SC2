import argparse
import sys

description=\
"This script takes a per residue contact table (such as the one outputted by residueContacts.py) and a text based sequence file (such as the one outputted by PDB_sequence_panClass.py) and \
merges both files such that each residue in a sequence has an associated contact annotation - if a contact is made. A single tabulated text file is outputted"

### File loader ###

def parse_args():
    """parse arguments for command line"""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--SeqTXT', dest = 'SeqTXT', type = str, required = True, help ='The tabulated sequence text file to be analysed')
    parser.add_argument('--ContactTXT', dest = 'ContactTXT', type = str, required = True, help ='The contact table text file to be analysed')
    parser.add_argument('--chains', dest = 'chains',   type = str, required = True, help ='Chains of TCR-pMHC complex in order MHCa,MHCb,peptide,TCRa,TCRb')
    
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
    
### Parsers ####
def everythingParser(inFile):
    everythingList=[]
    for lines in inFile.readlines():
        everythingList.append(lines)
    
    print("Input file contains " + str(len(everythingList)) + ' lines')
    return everythingList
    
def matrixParser(allLines):
    matrix=[]
    for lines in allLines:
        splitLine = ''
        splitLine=lines.split('\t')
        matrix.append(splitLine[:-1])
    return matrix
    
    ##### MATHS #######

def addIndexCodeSeq(line):
    index=''
    index+=line[0]+line[2]+line[3]
    line.append(index)    
    return line
    
def addIndexCodeSeqAll(seqMatrix):
    seqMatrixNew=[]
    newLine=[]
    for x in seqMatrix:
        newLine=addIndexCodeSeq(x)
        seqMatrixNew.append(newLine)
    return seqMatrixNew
    
def addIndexCodeContact(line):
    index=''
    index+=line[4]+line[6]+line[7]
    line.append(index)    
    return line
    
def addIndexCodeContactAll(conMatrix):
    conMatrixNew=[]
    newLine=[]
    for x in conMatrix:
        newLine=addIndexCodeSeq(x)
        newLine=addIndexCodeContact(x)
        conMatrixNew.append(newLine)
    return conMatrixNew
        
def pairContacts2parents2(contactMatrix,fullSeqMatrix):
    newFullSeqMatrix=[]
    for x in fullSeqMatrix:
        newSeqLine=[]
        hitCount=0
        for y in contactMatrix: 
            if x[4] == y[11]:
                newSeqLine=[]
#                 print("seq matrix in:")
#                 print(x)
#                 print("contact matrix goes in:")
#                 print(y)
                hitCount+=1
                for terms in x:
                    newSeqLine.append(terms)
                newSeqLine.append(y[4])
                newSeqLine.append(y[5])
                newSeqLine.append(y[6])
                newSeqLine.append(y[7])
                newSeqLine.append(y[8])
                newSeqLine.append(y[9])
                newSeqLine.append(y[10])
                newSeqLine.append(y[11])
                newSeqLine.append(y[12])  
                y[11],y[12] = "Done", "Done"
#                 print("after done:")
#                 print(y)
#                 print("output line:")
#                 print(newSeqLine)
#                 print("\n\n")
                newFullSeqMatrix.append(newSeqLine)
        if hitCount == 0:
            newSeqLine=[]
            for terms in x:
                newSeqLine.append(terms)
            newSeqLine.append('')
            newSeqLine.append('')
            newSeqLine.append('')
            newSeqLine.append('')
            newSeqLine.append('')
            newSeqLine.append('')
            newSeqLine.append('')
            newSeqLine.append('')
            newSeqLine.append('')
            newFullSeqMatrix.append(newSeqLine)           

    return newFullSeqMatrix

def fillAcceptorAnnot(seqMatrix):
    seqMatrixSearch=seqMatrix
    seqMatrixFix=seqMatrix
    for x in seqMatrixFix:
        if len(x[5]) != 0:
            for y in seqMatrixSearch:
                if y[4]==x[13]:
                    x[6]=y[1]    
    return seqMatrixFix
    
def removeIndex(seqMatrix):
    new=[]
    for line in seqMatrix:
        del line[4]
        del line[-2:]
        new.append(line)
    return new

#####       INITIALISER        #######
print('\n''     ~  Running fullResidueContacts.py v0.1  ~')

args = parse_args()
seqFile = args.SeqTXT
contactFile = args.ContactTXT

seqInFile = readFile(seqFile, "txt")
seqInFileName=seqFile.rsplit('.', 1)[0]

contactInFile = readFile(contactFile, "txt")
contactInFileName = contactFile.rsplit('.', 1)[0]
# Openers #
if type(seqInFileName) != str:
    raise IOError('No file was loaded. Please view usage and provide a valid file to be processed')
if type(contactInFileName) != str:
    raise IOError('No file was loaded. Please view usage and provide a valid file to be processed')

#Prep sequence file#
print("\nParsing sequence information..")
allSeqLines=everythingParser(seqInFile)
sequenceList=matrixParser(allSeqLines)
sequenceMatrix=sequenceList[1:]
sequenceMatrix=addIndexCodeSeqAll(sequenceMatrix)
print("Done!\n")

#Prep contact file#
print("Parsing contact information..")
allContactLines=everythingParser(contactInFile)
contactList=matrixParser(allContactLines)
contactMatrix=contactList[1:]
contacMatrix=addIndexCodeContactAll(contactMatrix)
print("Done!\n")

#Find pairs#

pairedSeqMatrix=pairContacts2parents2(contactMatrix, sequenceMatrix)
pairedSeqMatrix=fillAcceptorAnnot(pairedSeqMatrix)
outputContactMatrix=removeIndex(pairedSeqMatrix) 
for x in pairedSeqMatrix:
    print(x)
#####   OUTPUT GENERATOR       #######   

outFile = open(str(contactInFileName) + '_contacts_residues_full.txt' , 'w')

output2=''
output2+= 'Chain \tAnnotation \tResNum \tResCode \tChain \tAnnotation \tResNum \tResCode \tvdW \tHB \tSB \n'

for x in outputContactMatrix:
    for y in x:
        output2+=str(y) + '\t'
    output2+='\n'

outFile.write(output2)
print('\nOutputted file: ' + str(contactInFileName)+ '_contacts_residues_full.txt')
print('\n''     ~  End fullResidueContacts.py v0.1      ~')


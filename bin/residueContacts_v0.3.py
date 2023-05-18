import argparse
import sys

description=\
"This script takes a contact table such as the one outputted by projectContacts.py and reduces all contacts down to a per residue basis i.e. Residue x contacts y via a,b & c contact types. \
Example: TCRa        28    V    MHCb        81    H    True    False    False. This becomes important for downsteam analysis of contacts between two macromolecules. The script outputs one \
tabulated text file."

### File loader ###

def parse_args():
    """parse arguments for command line"""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--ContactTXT', dest = 'inputFile', type = str, required = True, help ='The contact table text file to be analysed')
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
        everythingList.append(lines[:-1])
    
    print('\n'+ "Input file contains " + str(len(everythingList)) + ' lines')
    return everythingList
    
def matrixParser(allLines):
    matrix=[]
    for lines in allLines:
        splitLine = ''
        splitLine=lines.split('\t')
        matrix.append(splitLine)
    return matrix
    
#####           MATHS           #######

def dataStripper(contactMatrix):
    newMatrix=[]
    for x in contactMatrix:
        y=[]
        z=''
        y.append(x[0])
        z+=x[0]
        y.append(x[2])
        z+=x[2]
        y.append(x[4])
        z+=x[4]
        y.append(x[7])
        z+=x[7]
        y.append(x[9])
        z+=x[9]
        y.append(x[11])
        z+=x[11]
        y.append(x[15])
        y.append(z)
            
        newMatrix.append(y)
    return newMatrix

def contactGrouper(contactMatrix):
    """
    
    Returns a three layer list: grouped contacts -> contacts -> contact parameters
    
    """
    print("\nReducing contacts to single residue level..")
    from itertools import groupby
    
    result=[]
    for key,group in groupby(contactMatrix, lambda x: x[7]):
        result.append(list(group))
    #print(result)
    return result    
    
    
def vDwHBSB(group):
    """
    
    Returns three bools: VW, HB, SB for each group
    
    """
  
    vdw=False
    HB=False
    SB=False
    for x in group:
        if 'VW' in x[6]:
            vdw=True
        if 'HB' in x[6]:
            HB=True
        if 'SB' in x[6]:
            SB=True 
    return vdw, HB, SB
    
def newGroup(group):
    
    x=[]
    y=[]
    y=group[0]
    vdw, HB, SB = vDwHBSB(group)
    x.append(y[0])
    x.append('')
    x.append(y[1])
    x.append(y[2])
    x.append(y[3])
    x.append('')
    x.append(y[4])
    x.append(y[5])
    x.append(vdw)
    x.append(HB)
    x.append(SB)
    return x
    
def simpleContactMatrix(groupedContactMatrix):
    
    print("\nGenerating single residue contact matrix..")
    
    y=[]
    for x in groupedContactMatrix:
        #print(x)
        y.append(newGroup(x))
    return y  
     
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
        
def replaceThree2one(aa):
    old=[]
    old=aa
    new=[]
    for x in old:
        x[3]=three2one(x[3])
        x[7]=three2one(x[7])
        
        new.append(x)
        
    return new

def nameChains(line,where):
    oldCode=line[where]
    newLine=line
    newCode=''
    if oldCode == chains_dict["MHCA"]:
        newCode ='MHCA'
    if oldCode == chains_dict["MHCB"]:
        newCode ='MHCB'
    if oldCode == chains_dict["PEPTIDE"]:
        newCode ='peptide'
    newLine[where]=newCode              
    return newLine
    

#####       INITIALISER        #######
# Openers #
print('\n''     ~  Running projectContactsSimple.py v0.2  ~')
args = parse_args()
fileIn = args.inputFile
chains = args.chains

inFile = readFile(fileIn, "txt")
inFileName=fileIn.rsplit('.', 1)[0]

if type(inFileName) != str:
    raise IOError('No file was loaded. Please view usage and provide a valid file to be processed')
allLines=everythingParser(inFile)
contactMatrix=matrixParser(allLines)
contactMatrix=contactMatrix[1:]

chains_dict = { "MHCA" : chains[0],
                "MHCB" : chains[1],
                "PEPTIDE" : chains[2]
}

#####           BODY           #######
contactMatrix=dataStripper(contactMatrix)

groupedContactMatrix=contactGrouper(contactMatrix)
outputContactMatrix=simpleContactMatrix(groupedContactMatrix)
outputContactMatrix2=outputContactMatrix[:]
outputContactMatrix2=replaceThree2one(outputContactMatrix2)

for x in outputContactMatrix:
    nameChains(x,0)
    nameChains(x,4)
for x in outputContactMatrix:
    print(x)



#####   OUTPUT GENERATOR       #######   

outFile = open(str(inFileName) + '_residues.txt' , 'w')

output2=''
output2+= 'Chain \tAnnotation \tResNum \tResCode \tChain \tAnnotation \tResNum \tResCode \tvdW \tHB \tSB \n'

for x in outputContactMatrix:
    for y in x:
        output2+=str(y) + '\t'
    output2+='\n'

outFile.write(output2)
print('\nOutputted file: ' + str(inFileName)+ '_residues.txt')
print('\n''     ~  End ProjectContactsSimple.py v0.1      ~')
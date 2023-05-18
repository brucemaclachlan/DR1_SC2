import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys
import itertools
import os
import re

description=\
"This script takes an ouput text from projectContacts which describes pMHC contact interactions (currently __ncont_out_contacts_residues_contacts_residues_full.txt)\
The script parses the annotated contact information and generates a networkx map where residues between the donor (MHC) and acceptor (p) act as nodes.\
Contacts between residues (nodes) are generated as edges. The resulting nodes and edges are drawn on a matplotlib plot. The script takes twp arguments which provide the \
input file and the MHC class (I/II)"

### File loader ###

def parse_args():
    """parse arguments for command line"""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--Input', dest = 'TXT', type = str, required = True, help ='The contact table text file to be analysed')
    parser.add_argument('--MHCclass', dest = 'MHCclass', type = str, required = True, help ='The MHC class of the complex structure: I or II')

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

##### Make index codes #######

def addIndexCodeDonor(line):
    index=''
    index+='%-7s' % (line[0],)
    index+='%-5s' % (line[1],)
    index+='%5s' % (line[2],)
    index+='%1s'% line[3]   
    return index
        
def addIndexCodeAcceptor(line):
    index=''
    index+='%-7s' % (line[4],)
    index+='%-5s' % (line[5],)
    index+='%5s' % (line[6],)
    index+='%1s'% line[7]   
    return index
    
def simplifiedMatrix(conMatrix):
    output=[]
    for x in conMatrix:
        newLine=[]
        newLine.append(addIndexCodeDonor(x))
        newLine.append(addIndexCodeAcceptor(x))
        newLine.append(x[8])
        newLine.append(x[9])
        newLine.append(x[10])
        output.append(newLine)
    return output

## Node functions ##
    
def removeEdgelessDonors(workingMatrix, donorList, acceptorList):
    print(donorList)
    passList = []
    for group in donorList:
        found = 0
        for line in workingMatrix:
            if group in line[0]:
                for find in acceptorList:
                    if re.search(find, line[1]) != None:
                        found += 1
        if found > 0:
            passList.append(group)
                    
                
    print(passList)
    return passList
        
def removeEdgelessAcceptors(workingMatrix, acceptorList, donorList):
    print(acceptorList)
    passList = []
    for group in acceptorList:
        found = 0
        for line in workingMatrix:
            if group in line[1]:
                for find in donorList:
                    if re.search(find, line[0]) != None:
                        found += 1
        if found > 0:
            passList.append(group)
                    
                
    print(passList)
    return passList    
    
   
def createNodeList(workingMatrix, groupList, f):
    output = []
    gapCount=1
    for group in groupList:
        out = []
        names = []
        for line in workingMatrix:
            if group in line[0]:
                if line[0] not in names:
                    names.append(line[0])
        print("Number of "+group+" nodes is " + str(len(names)))
        out.append(group)
        out.append(len(names))
        out.append(names)
        output.append(out)
        
        out = []
        names = []
        out.append("gap")
        out.append(f)
        for i in range(1,f+1):
            names.append("gap"+str(gapCount))
            gapCount+=1
        out.append(names)
        output.append(out)

    return output[:-1]


def reverseNodes(nodeList):
    output = []
    for group in nodeList:
        newGroup = []
        newGroup.append(group[0])
        newGroup.append(group[1])
        newGroup.append(list(reversed(group[2])))
        output.append(newGroup)
    return output

    
def dishOutLocations(nodeList, i, xcentre, y): 
    
    output = []
    numberOfNodes = 0.0
    for element in nodeList:
        numberOfNodes+=element[1]
    print("\nThere are a total of ", numberOfNodes, " nodes to dish out (including gap nodes)")
    xcoord=xcentre
    
    print("The middle of these nodes will be at "+ str(xcoord))
    print("i is", i)
    leftnright = 0.0
    leftnright=(numberOfNodes-1)/2
    xcoord+=leftnright*(i)*-1
    
    print("The far left of these nodes will be at "+ str(xcoord) + "\n") 
    print("Dishing out node locations..")
    
    for group in nodeList:
        counter = 1
        groupNodes = group[1]
        line = []
        line.append(group[0])
        line.append(groupNodes)
        line.append(group[2])
        positions = []
        while counter < groupNodes+1:
            positions.append(np.array([xcoord,y]))
            xcoord+=i
            counter+=1
        line.append(positions)
        output.append(line)

    print("Node locations established!\n")
    counter = 1
    print("Num", "\t", "Annot", "\t", "Res", "\t", "xcoord")
    for element in output:
        for node,coord in zip(element[2],element[3]):
            print(counter, "\t", element[0][:4], "\t", node[-1],"\t", coord[0])
            counter+=1
        
    return output

def generateNodes(groupList, nodedictionary, nodelabeldictionary, yOffset):
    for group in groupList:
        tracker=0
        print(group[0]+" nodes..")
        for node in group[2]:
            nodedictionary[node] = group[3][tracker]+np.array([0,yOffset])
            if group[0] != "gap":
                nodelabeldictionary[node] = group[2][tracker][-1]
            tracker+=1
    print("Done!")
    return None      

## Edge functions ##

def generateEdges(workingMatrix, donorsList, acceptorsList, bondType):
    print("\nGenerating edges between "+' '.join(donorsList)+" and "+' '.join(acceptorsList)+" of type "+str(bondType)+"..")
    output = []
    if bondType == None:
        for line in workingMatrix:
            if any(x in line[0] for x in donorsList):
                if any(x in line[1] for x in acceptorsList):
                    e = (line[0], line[1])
                    output.append(e)
    if bondType == "vdW":
        for line in workingMatrix:
            if any(x in line[0] for x in donorsList):
                if any(x in line[1] for x in acceptorsList):
                    if line[2] == "True":
                        e = (line[0], line[1])
                        output.append(e)    
    if bondType == "HB":
        for line in workingMatrix:
            if any(x in line[0] for x in donorsList):
                if any(x in line[1] for x in acceptorsList):
                    if line[3] == "True":
                        e = (line[0], line[1])
                        output.append(e)  
    if bondType == "SB":
        for line in workingMatrix:
            if any(x in line[0] for x in donorsList):
                if any(x in line[1] for x in acceptorsList):
                    if line[4] == "True":
                        e = (line[0], line[1])
                        output.append(e)
    print("Number of edges found is "+str(len(output)))
    return output

## Annotation functions
    
def generateAnnotations(groupList, offset):  
    output = []
    
    for group in groupList:
        if group[0] != "gap":
            left = int(group[2][0][13:-1])
            right = int(group[2][-1][13:-1])
            middle = group[0]
            leftLoc = group[3][0]+np.array([0,offset])
            rightLoc = group[3][-1]+np.array([0,offset])
            middleLoc = (leftLoc+rightLoc)/2
            
            if "a" in middle:
                middle = middle.replace("a", r"$\alpha$")
                print(middle)
            if "b" in middle:
                middle = middle.replace("b", r"$\beta$")  
                print(middle)
                
            line=[]
            line.append(middle)
            line.append(middleLoc)
            output.append(line)              
            line=[]
            line.append(left)
            line.append(leftLoc)
            output.append(line)
            
            line=[]
            line.append(right)
            line.append(rightLoc)
            output.append(line)
            
    return output

########## BODY #############################
    
#####       INITIALISER        #######
print('\n''     ~  Running contactMap_pMHC.py v1.0  ~')

args = parse_args()
seqFile = args.TXT
MHCclass = args.MHCclass
seqInFile = readFile(seqFile, "txt")
fileName=seqFile.rsplit('.', 1)[0]
filePath = fileName.rsplit('/', 1)[0]

mapPath = filePath+"/maps"

# Openers #
if type(fileName) != str:
    raise IOError('No file was loaded. Please view usage and provide a valid file to be processed')

# Make directories #
if not os.path.exists(mapPath):
    print("Creating Directory "+mapPath)
    os.makedirs(mapPath)
    
if not os.path.exists(mapPath+"/svg"):
    print("Creating Directory "+mapPath+"/svg")
    os.makedirs(mapPath+"/svg")
    
if not os.path.exists(mapPath+"/png"):
    print("Creating Directory "+mapPath+"/png")
    os.makedirs(mapPath+"/png")

if not os.path.exists(mapPath+"/pdf"):
    print("Creating Directory "+mapPath+"/pdf")
    os.makedirs(mapPath+"/pdf")

#Prep sequence file#

print("\nParsing sequence information..")
allSeqLines=everythingParser(seqInFile)
sequenceList=matrixParser(allSeqLines)
sequenceMatrix=sequenceList[1:]
workingMatrix=simplifiedMatrix(sequenceMatrix)
print("Done!\n")

# Is donors or acceptors longer?

######################### These parameters modify the size, scale and axes of the plot #######################

#  Set the boundaries of the plot #
xboundaries = (-3.0,3.0)
yboundaries = (0.0,2.0)

xsize = xboundaries[1]-xboundaries[0]
ysize = yboundaries[1]-yboundaries[0]

# x Increment between nodes
i = 0.09
# Increment between groups of nodes - number of fake nodes to create spaces
f = 2

font_size = 9
edge_width = 0.75

# x center location
xcentre = np.mean(xboundaries)

# Node offset to make nodes sit just below node labels (i.e. resdiue letter)
# This will need to be changed if y axis height is changed (OR FONT SIZE)

MHC1LabelOffset = -0.08
MHC2LabelOffset = +0.10
peptideOffsetUp = +0.10
peptideOffsetDown = -0.08

# Annotation (i.e. peptide/CDR label and residue number) offset to allow annotations to sit above/below residue letters
# This will need to be changed if y axis height is changed

MHC1Offset = +0.12
MHC2Offset = -0.16


# y axis height of the donor and acceptor lines

MHC1y = 1.6
MHC2y = 0.4
peptidey = 1.0

# Draw axis? on/off

axisYorN = "off"
#axisYorN = "on"

############################################################################################################
if MHCclass == "I":
    MHC1donorsList = ["MHCa1"]
    MHC2donorsList = ["MHCa2"]
    
    MHC1direction = "forward"
    MHC2direction = "forward"

if MHCclass == "II":
    MHC1donorsList = ["MHCa1"]
    MHC2donorsList = ["MHCb1"]
    
    MHC1direction = "forward"
    MHC2direction = "reverse"
    
acceptorsList = ["peptide"]
peptidedirection = "forward"


# determine the number of donor (MHC) nodes for each group
print("\nDetermining MHC node locations..\n")
MHC1donorsNodeList= createNodeList(workingMatrix, MHC1donorsList, f)
MHC2donorsNodeList= createNodeList(workingMatrix, MHC2donorsList, f)

if MHC2direction == "reverse":
    MHC2donorsNodeList = reverseNodes(MHC2donorsNodeList)


# Determine the donor graph locations
MHC1donors = dishOutLocations(MHC1donorsNodeList, i, xcentre, MHC1y)
MHC2donors = dishOutLocations(MHC2donorsNodeList, i, xcentre, MHC2y)

# determine the number of acceptor (middle line) nodes for each group
print("\nDetermining acceptor node locations..\n")
acceptorsNodeList= createNodeList(workingMatrix, acceptorsList, f)

# Determine the acceptor graph locations
acceptors = dishOutLocations(acceptorsNodeList, i, xcentre, peptidey)

# Generate Nodes
 
nodesToAdd = {}
nodesToAddLabels = {}
print("\nGenerating donor nodes..")
generateNodes(MHC1donors, nodesToAdd, nodesToAddLabels, 0)
generateNodes(MHC2donors, nodesToAdd, nodesToAddLabels, 0)

print("\nGenerating acceptor nodes..")
generateNodes(acceptors, nodesToAdd, nodesToAddLabels, 0)

# Generate Edges

#edgestoadd = generateEdges(workingMatrix, donors, acceptors, None)

MHC1vdWedges = generateEdges(workingMatrix, MHC1donorsList, acceptorsList, "vdW")
MHC1HBedges = generateEdges(workingMatrix, MHC1donorsList, acceptorsList, "HB")
MHC1SBedges = generateEdges(workingMatrix, MHC1donorsList, acceptorsList, "SB")

MHC2vdWedges = generateEdges(workingMatrix, MHC2donorsList, acceptorsList, "vdW")
MHC2HBedges = generateEdges(workingMatrix, MHC2donorsList, acceptorsList, "HB")
MHC2SBedges = generateEdges(workingMatrix, MHC2donorsList, acceptorsList, "SB")


# Generate annotations
MHC1donorAnnotations = generateAnnotations(MHC1donors, MHC1Offset)
MHC2donorAnnotations = generateAnnotations(MHC2donors, MHC2Offset)

#acceptorAnnotations = generateAnnotations(acceptors, acceptorAnnotationOffset)

# Generate plot and axes
fig = plt.figure(1, figsize=(xsize,ysize), dpi=300, frameon=False)
ax = plt.axes()

# plot and axes parameters

ax.set_xlim(xboundaries)
ax.set_ylim(yboundaries)

ax.axis(axisYorN)

# Generate network
G = nx.DiGraph()
G.add_nodes_from(nodesToAdd)
#G.add_edges_from(edgestoadd)

# Draw network
#nx.draw_networkx_nodes(G, nodesToAdd, ax=ax, node_size = 10, node_shape="s", node_color="r")

nx.draw_networkx_labels(G, nodesToAdd, ax=ax, labels = nodesToAddLabels, font_family="monospace", font_size=font_size)

# offset nodes to account for text height problem

# REgenerate nodes with an offset on the y direction

nodesToAdd = {}
nodesToAddLabels = {}

print("\nRegenerating donor nodes with an offset for text labels..")
generateNodes(MHC1donors, nodesToAdd, nodesToAddLabels, MHC1LabelOffset)

print("\nRegenerating acceptor nodes with an offset for text labels..")
generateNodes(acceptors, nodesToAdd, nodesToAddLabels, peptideOffsetUp)

nx.draw_networkx_edges(G, nodesToAdd, ax=ax, edgelist=MHC1vdWedges, edge_color='gray', style="solid", arrows=False, width=edge_width)
nx.draw_networkx_edges(G, nodesToAdd, ax=ax, edgelist=MHC1HBedges, edge_color='blue', style="dashed", arrows=False, width=edge_width)
nx.draw_networkx_edges(G, nodesToAdd, ax=ax, edgelist=MHC1SBedges, edge_color='red', style="dotted", arrows=False, width=edge_width)

nodesToAdd = {}
nodesToAddLabels = {}

print("\nRegenerating donor nodes with an offset for text labels..")
generateNodes(MHC2donors, nodesToAdd, nodesToAddLabels, MHC2LabelOffset)

print("\nRegenerating acceptor nodes with an offset for text labels..")
generateNodes(acceptors, nodesToAdd, nodesToAddLabels, peptideOffsetDown)



nx.draw_networkx_edges(G, nodesToAdd, ax=ax, edgelist=MHC2vdWedges, edge_color='gray', style="solid", arrows=False, width=edge_width)
nx.draw_networkx_edges(G, nodesToAdd, ax=ax, edgelist=MHC2HBedges, edge_color='blue', style="dashed", arrows=False, width=edge_width)
nx.draw_networkx_edges(G, nodesToAdd, ax=ax, edgelist=MHC2SBedges, edge_color='red', style="dotted", arrows=False, width=edge_width)
# Annotations
for annot in MHC1donorAnnotations:
    ax.annotate(annot[0], (annot[1][0], annot[1][1]), ha="center", size=font_size-3, family="monospace")
for annot in MHC2donorAnnotations:
    ax.annotate(annot[0], (annot[1][0], annot[1][1]), ha="center", size=font_size-3, family="monospace")

#for annot in acceptorAnnotations:    
#    ax.annotate(annot[0], (annot[1][0], annot[1][1]), ha="center", size=6)

#Show plot
#plt.axvline(x=0)
#plt.show()

outputName = "pMHC_contactMap"

    

plt.savefig(mapPath+"/svg/"+outputName+".svg", format="svg", transparent=True)
plt.savefig(mapPath+"/png/"+outputName+".png", format="png", transparent=True)
plt.savefig(mapPath+"/pdf/"+outputName+".pdf", format="pdf", transparent=True)
print('\n''     ~  End of contactMap_pMHC.py v1.0  ~')
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

### File loader ###

description=""


import argparse
import sys
import os

# use latex text rendering
from matplotlib import rc

from collections import Counter
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

import pandas as pd

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

mpl.rcParams['text.latex.preamble'] = r'\usepackage{siunitx} \sisetup{detect-all} \usepackage{helvet} \usepackage{sansmath} \sansmath'
# mpl.rcParams['axes.unicode_minus'] = False


def parse_args():
    """parse arguments for command line"""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--ContactTXT', dest = 'inputFile', type = str, required = True, help ='The contact table text file to be analysed')
    parser.add_argument('--SeqTXT', dest = 'SeqTXT', type = str, required = True, help ='The tabulated sequence text file to be analysed')
    parser.add_argument('--peptide_colour', dest = 'colour',   type = str, required = False, help ='Colour of the peptide to be used for plotting colour. Do not include the hash.. argsparse does not like', default = "808080")
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

def force_aspect(ax,aspect,xmin,xmax,ymin,ymax):
    
    ax.set_aspect(np.absolute((xmax-xmin)/(ymax-ymin))*(aspect[1]/aspect[0]))
    ax.set_xlim([xmin, xmax])
    ax.set_ylim([ymin, ymax])  
    return None

def bruce_magic_axis(ax, aspect, xmin, xmax, xticks, ymin, ymax, yticks):
    
    force_aspect(ax, aspect, xmin, xmax, ymin, ymax)
    ax.xaxis.set_ticks(xticks)
    ax.yaxis.set_ticks(yticks)            
    
    return None

def boxed_label(ax, text, x_left, x_right, y_center, pad, height, facecolor, alpha, fontsize):

    '''
    Creates a text box which spans a certain distance with text centered in the middle.

    ax: the matplotlib axis to add the text to
    text: the text string you want added
    x_left: the value you want the left hand side of the box to start at
    x_right: the value you want the right hand side of the box to end
    y_center: the value you want the middle of the box to be (and text) on the y axis
    pad: an offset value to remove from the edges of the box. Use this if you have two boxes next to each other
    height: the height of the box


    if you want to work in axes values delete transform=ax.transAxes (hacky for now)

    e.g.

    boxed_label(ax=ax, text='Spike', x_left=ax.get_xlim()[0], x_right=15.0, y_center=-5.75, pad=0.25, height=1.0)
    boxed_label(ax=ax, text='non-Spike', x_left=15.0, x_right=ax.get_xlim()[1], y_center=-5.75, pad=0.25, height=1.0)

    '''
    import matplotlib.patches as patches

    ax.add_patch(patches.Rectangle(xy=(x_left+pad , y_center-(height/2.0) ), width=x_right-x_left-pad, height=height, linewidth=0, clip_on=False, facecolor=facecolor, alpha=alpha, transform=ax.transAxes))
    ax.text(x=(x_left+x_right)/2.0, y=y_center, s=text, fontsize=fontsize, horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)

#####       INITIALISER        #######
# Openers #
print('\n''     ~  Running contacts_graphs.py  ~')
args = parse_args()
fileIn = args.inputFile
seqFile = args.SeqTXT
output_dir = args.output_dir

# args parse doesn't seem to like the hash in the hex code, so 
peptide_colour = "#"+args.colour

inFile = readFile(fileIn, "txt")
inFileName=fileIn.rsplit('/')[0]

seqInFile = readFile(seqFile, "txt")
seqInFileName=seqFile.rsplit('.', 1)[0]


if type(inFileName) != str:
    raise IOError('No file was loaded. Please view usage and provide a valid file to be processed')
allLines=everythingParser(inFile)
contactMatrix=matrixParser(allLines)
contactMatrix=contactMatrix[1:]


#Prep sequence file#
print("\nParsing sequence information..")
allSeqLines=everythingParser(seqInFile)
sequenceList=matrixParser(allSeqLines)
sequenceMatrix=sequenceList[1:]

print("Done!\n")


#########################
# MATHS                 #
#########################

res_nums = []
res_codes = []

for line in sequenceMatrix:
	if line[0] == 'peptide':
		res_nums.append(int(line[2]))
		res_codes.append(line[3])

print(res_nums)
print(res_codes)

peptide_contacts = []
for line in contactMatrix:
    peptide_contacts.append(int(line[9]))


print("Peitide contacts...")
print(peptide_contacts)
peptide_contacts_counted = Counter(peptide_contacts)
print("Peitides counted...")
print(peptide_contacts_counted)


counted_contacts_full = {}

for pos in res_nums:
	if pos in peptide_contacts_counted:
		print(pos)
		print(peptide_contacts_counted[pos])
		counted_contacts_full[int(pos)] = peptide_contacts_counted[pos]

	else:
		counted_contacts_full[int(pos)] = 0

print("counted_contacts_full...")
print(counted_contacts_full)

counted_contacts_full_list = []
for pos in counted_contacts_full:
	counted_contacts_full_list.append(counted_contacts_full[pos])

total_contacts = 0

for position in counted_contacts_full:
    total_contacts += peptide_contacts_counted[position]

print("Total contacts is.." , total_contacts)

perc_contacts_full = {}

for position in counted_contacts_full:
    perc_contacts_full[int(position)] = (float(float(peptide_contacts_counted[position])/total_contacts*100))
    
print("perc_contacts_full...")
print(perc_contacts_full)


perc_contacts_full_list = []
for pos in perc_contacts_full:
	perc_contacts_full_list.append(perc_contacts_full[pos])



###########################################
# build a pandas dataframe                #
###########################################

df = pd.DataFrame.from_dict({
							'position' : res_nums,
							'residue'  : res_codes,
							'contacts' : counted_contacts_full_list,
							'perc_contacts' : perc_contacts_full_list
							})

print(df)

graphPath = output_dir+"/"+"contacts/pMHC_contacts/plots"

print(inFileName)
print(graphPath)



# Make directories #
if not os.path.exists(graphPath):
    print("Creating Directory "+graphPath)
    os.makedirs(graphPath)

#########################
# Raw contacts          #
#########################

# Initialise figure and subplots
nrow = 1
ncol = 1
fig, ax = plt.subplots(nrow, ncol, figsize=(2.625 , 1.75), sharex=False, sharey=False, facecolor="None")

label_size = 8

labels = df['residue'].tolist()
labels.insert(0,"")
labels.append("")



markerline, stemlines, baseline = ax.stem(df['position'], df['contacts'], basefmt = " ")

plt.setp(stemlines, color=peptide_colour, linewidth=2)
plt.setp(markerline, color=peptide_colour, linewidth=2, markersize=9.0)

for i in range(len(df['contacts'])):
	print(df['contacts'][i])
	ax.text(x=df['position'][i], y=df['contacts'][i], s=df['contacts'][i], size=label_size-1, va="center", ha="center", color="k", weight="bold", fontfamily="monospaced")


ax.text(x=df['position'].min()-1, y=-12.0, s="", fontsize=label_size, va="center", ha="center", color="k")
ax.text(x=df['position'].max()+1, y=-12.0, s="", fontsize=label_size, va="center", ha="center", color="k")

for index, row in df.iterrows():
    ax.text(x=row['position'], y=-8.5, s=row['residue'], fontsize=label_size, va="center", ha="center", color="k") 

boxed_label(ax, "p-HLA contacts (Total = "+str(total_contacts) + ")", x_left=0.0, x_right=1.0, y_center=1.080, pad=0.005, height=0.15, facecolor=peptide_colour, alpha=0.8, fontsize=label_size)
ax.text(x=-0.135, y=0.5, s="No. contacts", rotation=90.0, transform=ax.transAxes, fontsize=label_size, va="center", ha="center")


axes_vitals_contacts = {
        "aspect" :  (1.5,1.0),
        "xmin"  :  df['position'].min()-1,
        "xmax"  :  df['position'].max()+1,
        "xticks":  np.arange(df['position'].min(), df['position'].max()+1 , 1),
        "ymin"  :  -3.5,
        "ymax"  :  +45.0+5.0,
        "yticks":  np.arange(0.0, 45.0+1.0, 5.0)    
        }
bruce_magic_axis(ax=ax, **axes_vitals_contacts)
ax.spines['bottom'].set_position(("data", axes_vitals_contacts["ymin"]))








    
#stylise the frame
ax.set_facecolor("None")
ax.spines['left'].set_linewidth(0.75)
ax.spines['bottom'].set_linewidth(0.75)
ax.spines['right'].set_visible(0.75)
ax.spines['top'].set_visible(0.75)

print(labels)
ax.set_xticklabels([])

y_labels = ax.get_yticklabels()
plt.setp(y_labels, rotation=0, fontsize=label_size)


fig.savefig(graphPath+"/p_MHC_contacts.pdf", bbox_inches="tight", transparent=True)

###################################
#### Perentage total contacts #####
###################################

# Initialise figure and subplots
nrow = 1
ncol = 1
fig, ax = plt.subplots(nrow, ncol, figsize=(2.625, 1.75), sharex=False, sharey=False, facecolor="None")

label_size = 8

labels = df['residue'].tolist()
labels.insert(0,"")
labels.append("")

markerline, stemlines, baseline = ax.stem(df['position'], df['perc_contacts'], basefmt = " ")
plt.setp(stemlines, color=peptide_colour, linewidth=2)
plt.setp(markerline, color=peptide_colour, linewidth=2, markersize=9.0)

for i in range(len(df['perc_contacts'])):
    print(df['perc_contacts'][i])
    ax.text(x=df['position'][i], y=df['perc_contacts'][i], s='%d' % df['perc_contacts'][i], size=label_size-1, va="center", ha="center", color="k", weight="bold", fontfamily="monospaced")


ax.text(x=df['position'].min()-1, y=-3.5, s="", fontsize=label_size, va="center", ha="center", color="k")
ax.text(x=df['position'].max()+1, y=-3.5, s="", fontsize=label_size, va="center", ha="center", color="k")

for index, row in df.iterrows():
    ax.text(x=row['position'], y=-5.5, s=row['residue'], fontsize=label_size, va="center", ha="center", color="k") 

boxed_label(ax, "p-HLA contacts (\% Total)", x_left=0.0, x_right=1.0, y_center=1.080, pad=0.005, height=0.15, facecolor=peptide_colour, alpha=0.8, fontsize=label_size)
ax.text(x=-0.135, y=0.5, s="\% Total contacts", rotation=90.0, transform=ax.transAxes, fontsize=label_size, va="center", ha="center")


# this is a dump from peptide binding assays. Trying to get the plots to look similar style.



axes_vitals_perc = {
        "aspect" :  (1.5,1.0),
        "xmin"  :  df['position'].min()-1,
        "xmax"  :  df['position'].max()+1,
        "xticks":  np.arange(df['position'].min(), df['position'].max()+1, 1),
        "ymin"  :  -2.5,
        "ymax"  :  +20.0+5.0,
        "yticks":  np.arange(0.0, 20.0+1.0, 5.0)    
        }
bruce_magic_axis(ax=ax, **axes_vitals_perc)
ax.spines['bottom'].set_position(("data", axes_vitals_perc["ymin"]))


#stylise the frame
ax.set_facecolor("None")
ax.spines['left'].set_linewidth(0.75)
ax.spines['bottom'].set_linewidth(0.75)
ax.spines['right'].set_visible(0.75)
ax.spines['top'].set_visible(0.75)

print(labels)
ax.set_xticklabels([])

y_labels = ax.get_yticklabels()
plt.setp(y_labels, rotation=0, fontsize=label_size)

fig.savefig(graphPath+"/p_MHC_perc_contacts.pdf", bbox_inches="tight", transparent=True)



##################################
# HYDROGEN BONDS                 #
##################################

# For now I'm just going to copy and paste the above and modify it to pick out the h bonds. This is a little clunky but I will
# streamline the code later if interesting.

print("§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§ DID I MAKE HERE §§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§")

#########################
# MATHS                 #
#########################

res_nums = []
res_codes = []

for line in sequenceMatrix:
    if line[0] == 'peptide':
        res_nums.append(int(line[2]))
        res_codes.append(line[3])

print(res_nums)
print(res_codes)

peptide_contacts = []

######################### This is the key code that makes it HB specific Note for future when streamlining code ######################### 

for line in contactMatrix:
    if line[15] == "HB": # <----HERE
        peptide_contacts.append(int(line[9]))

######################################################################################################################################### 

print("Peitide contacts...")
print(peptide_contacts)
peptide_contacts_counted = Counter(peptide_contacts)
print("Peitides counted...")
print(peptide_contacts_counted)


counted_contacts_full = {}

for pos in res_nums:
    if pos in peptide_contacts_counted:
        print(pos)
        print(peptide_contacts_counted[pos])
        counted_contacts_full[int(pos)] = peptide_contacts_counted[pos]

    else:
        counted_contacts_full[int(pos)] = 0

print("counted_contacts_full...")
print(counted_contacts_full)

counted_contacts_full_list = []
for pos in counted_contacts_full:
    counted_contacts_full_list.append(counted_contacts_full[pos])

total_contacts = 0

for position in counted_contacts_full:
    total_contacts += peptide_contacts_counted[position]

print("Total contacts is.." , total_contacts)

perc_contacts_full = {}

for position in counted_contacts_full:
    perc_contacts_full[int(position)] = (float(float(peptide_contacts_counted[position])/total_contacts*100))
    
print("perc_contacts_full...")
print(perc_contacts_full)


perc_contacts_full_list = []
for pos in perc_contacts_full:
    perc_contacts_full_list.append(perc_contacts_full[pos])



###########################################
# build a pandas dataframe                #
###########################################

df = pd.DataFrame.from_dict({
                            'position' : res_nums,
                            'residue'  : res_codes,
                            'contacts' : counted_contacts_full_list,
                            'perc_contacts' : perc_contacts_full_list
                            })

print(df)

graphPath = output_dir+"/"+"contacts/pMHC_contacts/plots"

print(inFileName)
print(graphPath)



# Make directories #
if not os.path.exists(graphPath):
    print("Creating Directory "+graphPath)
    os.makedirs(graphPath)


#########################
# Raw contacts          #
#########################

# Initialise figure and subplots
nrow = 1
ncol = 1
fig, ax = plt.subplots(nrow, ncol, figsize=(2.625 , 1.75), sharex=False, sharey=False, facecolor="None")

label_size = 8

labels = df['residue'].tolist()
labels.insert(0,"")
labels.append("")



markerline, stemlines, baseline = ax.stem(df['position'], df['contacts'], basefmt = " ")

plt.setp(stemlines, color=peptide_colour, linewidth=2)
plt.setp(markerline, color=peptide_colour, linewidth=2, markersize=9.0)

for i in range(len(df['contacts'])):
    print(df['contacts'][i])
    ax.text(x=df['position'][i], y=df['contacts'][i], s=df['contacts'][i], size=label_size-1, va="center", ha="center", color="k", weight="bold", fontfamily="monospaced")


ax.text(x=df['position'].min()-1, y=-12.0, s="", fontsize=label_size, va="center", ha="center", color="k")
ax.text(x=df['position'].max()+1, y=-12.0, s="", fontsize=label_size, va="center", ha="center", color="k")

for index, row in df.iterrows():
    ax.text(x=row['position'], y=-8.5, s=row['residue'], fontsize=label_size, va="center", ha="center", color="k") 

boxed_label(ax, "p-HLA HBs (Total = "+str(total_contacts) + ")", x_left=0.0, x_right=1.0, y_center=1.080, pad=0.005, height=0.15, facecolor=peptide_colour, alpha=0.8, fontsize=label_size)
ax.text(x=-0.135, y=0.5, s="No. contacts", rotation=90.0, transform=ax.transAxes, fontsize=label_size, va="center", ha="center")


axes_vitals_contacts = {
        "aspect" :  (1.5,1.0),
        "xmin"  :  df['position'].min()-1,
        "xmax"  :  df['position'].max()+1,
        "xticks":  np.arange(df['position'].min(), df['position'].max()+1 , 1),
        "ymin"  :  -3.5,
        "ymax"  :  +45.0+5.0,
        "yticks":  np.arange(0.0, 45.0+1.0, 5.0)    
        }
bruce_magic_axis(ax=ax, **axes_vitals_contacts)
ax.spines['bottom'].set_position(("data", axes_vitals_contacts["ymin"]))
    
#stylise the frame
ax.set_facecolor("None")
ax.spines['left'].set_linewidth(0.75)
ax.spines['bottom'].set_linewidth(0.75)
ax.spines['right'].set_visible(0.75)
ax.spines['top'].set_visible(0.75)

print(labels)
ax.set_xticklabels([])

y_labels = ax.get_yticklabels()
plt.setp(y_labels, rotation=0, fontsize=label_size)


fig.savefig(graphPath+"/p_MHC_HBs.pdf", bbox_inches="tight", transparent=True)

###################################
#### Perentage total contacts #####
###################################

# Initialise figure and subplots
nrow = 1
ncol = 1
fig, ax = plt.subplots(nrow, ncol, figsize=(2.625, 1.75), sharex=False, sharey=False, facecolor="None")

label_size = 8

labels = df['residue'].tolist()
labels.insert(0,"")
labels.append("")

markerline, stemlines, baseline = ax.stem(df['position'], df['perc_contacts'], basefmt = " ")
plt.setp(stemlines, color=peptide_colour, linewidth=2)
plt.setp(markerline, color=peptide_colour, linewidth=2, markersize=9.0)

for i in range(len(df['perc_contacts'])):
    print(df['perc_contacts'][i])
    ax.text(x=df['position'][i], y=df['perc_contacts'][i], s='%d' % df['perc_contacts'][i], size=label_size-1, va="center", ha="center", color="k", weight="bold", fontfamily="monospaced")


ax.text(x=df['position'].min()-1, y=-3.5, s="", fontsize=label_size, va="center", ha="center", color="k")
ax.text(x=df['position'].max()+1, y=-3.5, s="", fontsize=label_size, va="center", ha="center", color="k")

for index, row in df.iterrows():
    ax.text(x=row['position'], y=-5.5, s=row['residue'], fontsize=label_size, va="center", ha="center", color="k") 

boxed_label(ax, "p-HLA contacts (\% Total)", x_left=0.0, x_right=1.0, y_center=1.080, pad=0.005, height=0.15, facecolor=peptide_colour, alpha=0.8, fontsize=label_size)
ax.text(x=-0.135, y=0.5, s="\% Total HBs", rotation=90.0, transform=ax.transAxes, fontsize=label_size, va="center", ha="center")


# this is a dump from peptide binding assays. Trying to get the plots to look similar style.



axes_vitals_perc = {
        "aspect" :  (1.5,1.0),
        "xmin"  :  df['position'].min()-1,
        "xmax"  :  df['position'].max()+1,
        "xticks":  np.arange(df['position'].min(), df['position'].max()+1, 1),
        "ymin"  :  -2.5,
        "ymax"  :  +20.0+5.0,
        "yticks":  np.arange(0.0, 20.0+1.0, 5.0)    
        }
bruce_magic_axis(ax=ax, **axes_vitals_perc)
ax.spines['bottom'].set_position(("data", axes_vitals_perc["ymin"]))


#stylise the frame
ax.set_facecolor("None")
ax.spines['left'].set_linewidth(0.75)
ax.spines['bottom'].set_linewidth(0.75)
ax.spines['right'].set_visible(0.75)
ax.spines['top'].set_visible(0.75)

print(labels)
ax.set_xticklabels([])

y_labels = ax.get_yticklabels()
plt.setp(y_labels, rotation=0, fontsize=label_size)

fig.savefig(graphPath+"/p_MHC_perc_HBs.pdf", bbox_inches="tight", transparent=True)
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 10:31:55 2019
Python 2.7 Anaconda Distribution recommended
Written using Spyder

@author:    Bruce J MacLachlan, Monash Biomedicine Discovery Institute, Clayton, VIC, 3800
            ORCID ID: 0000-0002-2685-2733
            Github: https://github.com/brucemaclachlan
            email: bruce.maclachlan@monash.edu
            
            BSA functions
            
"""
def force_aspect(ax,aspect,xmin,xmax,ymin,ymax):
    import numpy as np
    
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

def is_float_or_int(string):
    try:
        int(string)
        return int(string)
    except ValueError:
        try:
            float(string)
            return float(string)
        except ValueError:
            return string

def parsePISA(pisa_output):
    '''
    Parses the output of a PISA request for interface information. PISA, why you no output as csv.
    '''
    print("Parsing PISA interface output to make a dict of interfaces..")
    output = {}
    headerline = ''
    dataList = []
    dataList2 = []
    for line in pisa_output:
        if "## Id" in line:
            headerline = line
        if any(i.isdigit() for i in line[0:4]):
            dataLine = line[:-1].split(" ")
            while "" in dataLine: 
                dataLine.remove("")
            while "|" in dataLine: 
                dataLine.remove("|")
            dataList.append(dataLine)
        
        for data in dataList:
            data2 = []
            for item in data:
                if "|" in item:
                    data2.append(item.replace("|", ""))
                else:
                    data2.append(item)
            dataList2.append(data2)
        

    headers = headerline[:-1].split(" ")
    while "" in headers: 
        headers.remove("")
    while "|" in headers: 
        headers.remove("|")
    headers = [head.replace('Symmetry', 'Symmetry.Op') for head in headers]
    while "operation" in headers: 
        headers.remove("operation")
    # parse these lists into a dictionary of interfaces (##) containing a dictionary of properties #
    
    index = 0
    while index < len(dataList2):
    
        interface_dict = {}
        
        indey = 0
        while indey < 10:
            
            if type(is_float_or_int(dataList2[index][indey])) == float:
                interface_dict[headers[indey]] = float(dataList2[index][indey])
            if type(is_float_or_int(dataList2[index][indey])) == int:
                interface_dict[headers[indey]] = int(dataList2[index][indey]) 
            if type(is_float_or_int(dataList2[index][indey])) == str:
                interface_dict[headers[indey]] = str(dataList2[index][indey])            
                  
            indey = indey + 1
        
        output[dataList2[index][0]] = interface_dict
        index = index + 1
    print("Done!")
    return output


def parsePISA_detail(detail_interface):

    d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
 'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
 'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
 'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

    detail_file = open(detail_interface , "r")
    Lines = detail_file.readlines()

    index = 0
    for line in Lines:
        # print(line)
        if "4. Interfacing Residues" in line:
            header = index+2
            first = index+4
            break
        index+=1

    index = first
    for line in Lines[index:]:
        if "----------------------" in line:
            last = index
            break
        index+=1

    data = Lines[first:last]
    header=Lines[header]
    print(header)

    for line in data:
        print(line)

    output = {}

    for line in data:
        line_dict = {}
        line_dict["##"] = int(line[0:5])
        line_dict["I"] = line[7]
        line_dict["Chain"] = line[10]
        line_dict["ResName"] =line[12:15]
        line_dict["ResCode"] =d[line[12:15]]

        line_dict["ResNum"] = int(line[16:20])


        line_dict["HS"] = line[22:24]
        line_dict["ASA"] = float(line[25:32])
        line_dict["BSA"] = float(line[32:39])
        line_dict["DeltaG"] =float(line[39:46])

        output[int(line[0:5])] = line_dict 


    print(output)
    return output

def remove_symm_interfaces(interface_dict):
    new_dict = interface_dict.copy()
    for interface in interface_dict:
        if interface_dict[interface]["Sym.Id"] != 1555:
            del new_dict[interface]
          
    return new_dict
        
def find_interface(interface_tuple, interface_dict):
    chain1 = interface_tuple[0]
    chain2 = interface_tuple[1]
    for interface in interface_dict:
        if interface_dict[interface]["Monomer1"] == chain1:
            if interface_dict[interface]["Monomer2"] == chain2:
                return interface
        elif interface_dict[interface]["Monomer2"] == chain1:
            if interface_dict[interface]["Monomer1"] == chain2:
                return interface  

# Call PISA

def call_pisa(sesh_name, pdb, output_dir):

    import subprocess

    bsa_dir = output_dir+"/bsa"


    ## Create pisa session -analyse
    command = " ".join(['pisa', sesh_name, '-analyse', pdb+".pdb"])
    print(command)
    subprocess.call([command], shell=True)


    # Create interface list
    command = " ".join(['pisa', sesh_name, '-list interfaces', '>' , bsa_dir+"/"+pdb+'_pisa_interfaces.txt'])
    print(command)
    subprocess.call([command], shell=True)

    return None

# Parse and clean the BSA interface information from PISA

def clean_pisa(pdb, output_dir):

    bsa_dir = output_dir+"/bsa"

    interfaces=open(bsa_dir+"/"+pdb+'_pisa_interfaces.txt', 'r')

    interface_dict = parsePISA(interfaces)

    interfaces_dict = remove_symm_interfaces(interface_dict)

    return interfaces_dict

def generate_bsa_areas(interfaces_of_interest, interfaces_dict, colours, labels):
    '''
    e.g. 
    interfaces_of_interest = [("A","J"),("P","J"),("A","K"),("P","K")]
    '''

    # output a list of dictionaries
    output = []
    for ioi, label, colour in zip(interfaces_of_interest, labels, colours):
        ioi_dict ={}
        # print(interfaces_dict[find_interface(ioi, interfaces_dict)])
        ioi_dict["name"] = ioi[0]+"_to_"+ioi[1]
        ioi_dict["label"] = label
        ioi_dict["colour"] = colour
        ioi_dict["area"] = interfaces_dict[find_interface(ioi, interfaces_dict)]["Area"]*2

        ioi_dict["num_label"] = ('%.0f' % ioi_dict["area"] + " "+ r"$\AA^2$")


        output.append(ioi_dict)

    return output


def detail_interface(sesh_name, pdbnn, interface, interfaces_dict, output_dir):
    '''
    e.g. 
    interfaces_of_interest = ("A","J")
    '''

    interface_index = find_interface(interface, interfaces_dict)

    print(interface_index)
    print(type(interface_index))
    print(interfaces_dict[interface_index])
    interface_id = interfaces_dict[interface_index]["##"]


    print("Calling PISA to detail the interface..")
    # Call pisa to detail that interface
    command = " ".join(['pisa', sesh_name, '-detail interfaces', str(interface_id), '>' , output_dir+"/bsa/"+pdbnn+'_pisa_detail_interface_'+str(interface_id)+'.txt'])
    print(command)

    import subprocess

    subprocess.call([command], shell=True)

    print("Parsing PISA detail interface file..")
    detail_interface_dict = parsePISA_detail(output_dir+"/bsa/"+pdbnn+'_pisa_detail_interface_'+str(interface_id)+'.txt')

    return detail_interface_dict


def calculate_bsa(sesh_name, pdb, interfaces_of_interest, labels, colours, do_pisa, output_dir):

    print("Calculating BSA using PISA..")


    pdbn = pdb.split("/")[-1]
    pdbnn = pdbn.split(".")[0]

    import os
    if not os.path.exists(pdbn):
        print("Copying pdb file to working directory..")
        import shutil
        shutil.copy2(pdb, pdbn)

    print(pdbn)
    analysis_dir = output_dir
    bsa_dir = output_dir+"/bsa"

    if not os.path.exists(analysis_dir):
        print("Creating Directory "+analysis_dir)
        os.makedirs(analysis_dir)

    if not os.path.exists(bsa_dir):
        print("Creating Directory "+bsa_dir)
        os.makedirs(bsa_dir)


    if do_pisa == True:
        print("Calling PISA..")
        call_pisa(sesh_name, pdbnn, output_dir)
        print("..DONE!")

    print("..DONE!")
    print("Cleaning PISA file..")
    print(bsa_dir)
    interfaces_dict = clean_pisa(pdbnn, output_dir)
    print("..DONE!")

    print("Generating dictionary of areas of interfaces of interest..")
    area_dict = generate_bsa_areas(interfaces_of_interest = interfaces_of_interest, interfaces_dict = interfaces_dict, labels = labels, colours = colours)
    print("..DONE!")

    print("Generating per residue data for the areas of interest")


    asa_bsa_dG_dicts = {}

    for interface in interfaces_of_interest:
        print(interface)
        asa_bsa_dG_dict = detail_interface(sesh_name, pdbnn, interface, interfaces_dict, output_dir)
        asa_bsa_dG_dicts["_".join(interface)] = asa_bsa_dG_dict



    # peptide_bsa_asa = pd.Dataframe(asa_bsa_dG_dicts)

    return area_dict , asa_bsa_dG_dicts




def plot_donut(bsa_plot_df, saveas):

    # use latex text rendering
    from matplotlib import rc
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib as mpl

    import pandas as pd
    import numpy as np

    rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    rc('text', usetex=True)

    mpl.rcParams['text.latex.preamble'] = r'\usepackage{siunitx} \sisetup{detect-all} \usepackage{helvet} \usepackage{sansmath} \sansmath'
    # mpl.rcParams['axes.unicode_minus'] = False

    # Begin plotting a donut

    total = sum(bsa_plot_df.area)
    label_size = 8

    nrow = 1
    ncol = 1
    fig, ax = plt.subplots(nrow, ncol, figsize=(2.5,2.5), dpi = 400, sharex=True, sharey=True, facecolor="w")

    size = 0.4
    outline_size = 0.5

    patches0i, texts0i = ax.pie(bsa_plot_df.area, radius=1-size, labeldistance=0.65, colors = bsa_plot_df.colour, wedgeprops=dict(width=size, edgecolor='w', linewidth=outline_size), labels = bsa_plot_df.num_label)
    
    #replace_small_labels(ax[0], texts0i)
    plt.setp(texts0i, size=label_size, va='center', ha='center')

    ax.set_ylim((-1,1))
    ax.set(aspect="equal")

    ax.text(0.5, 0.125, 'Total BSA = '+'%.0f' % total + " "+ r"$(\mathrm{\AA}^2)$", horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize = label_size*1.2)

    fig.savefig(saveas, bbox_inches="tight", transparent=True)

def plot_bsa_asa(df, colour, saveas):

    # use latex text rendering
    from matplotlib import rc

    from collections import Counter
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker

    import pandas as pd
    import numpy as np

    rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    rc('text', usetex=True)

    mpl.rcParams['text.latex.preamble'] = r'\usepackage{siunitx} \sisetup{detect-all} \usepackage{helvet} \usepackage{sansmath} \sansmath'
    # mpl.rcParams['axes.unicode_minus'] = False


    # Initialise figure and subplots
    nrow = 1
    ncol = 1
    fig, ax = plt.subplots(nrow, ncol, figsize=(2.625, 1.75), sharex=True, sharey=False, facecolor="None")

    label_size = 8

    # ASA

    bars = ax.bar(df['position'], df['ASA'], color=colour)
    bars = ax.bar(df['position'], df['BSA']*-1, color=colour)
    ax.axhline(color="k", linewidth=0.75)

    ticks =  ax.get_xticks()
    ax.set_xticklabels("")

    # for i in range(len(df['ASA'])):
    #     print(df['ASA'][i])
    #     ax.text(x=df['position'][i], y=df['ASA'][i], s=df['ASA'][i], size=label_size-1, va="center", ha="center", color="k", weight="bold", fontfamily="monospaced")


    ax.text(x=df['position'].min()-1, y=-12.0, s="", fontsize=label_size, va="center", ha="center", color="k")
    ax.text(x=df['position'].max()+1, y=-12.0, s="", fontsize=label_size, va="center", ha="center", color="k")

    for index, row in df.iterrows():
        ax.text(x=row['position'], y=-355.0, s=row['residue'], fontsize=label_size, va="center", ha="center", color="k") 

    boxed_label(ax, "Peptide ASA / BSA", x_left=0.0, x_right=1.0, y_center=1.080, pad=0.005, height=0.15, facecolor=colour, alpha=0.8, fontsize=label_size)
    ax.text(x=-0.19, y=0.75, s="ASA "+r"$(\mathrm{\AA}^2)$", rotation=90.0, transform=ax.transAxes, fontsize=label_size, va="center", ha="center")
    ax.text(x=-0.19, y=0.25, s="BSA "+r"$(\mathrm{\AA}^2)$", rotation=90.0, transform=ax.transAxes, fontsize=label_size, va="center", ha="center")



    axes_vitals = {
        "aspect" :  (1.5,1.0),
        "xmin"  :  df['position'].min()-1,
        "xmax"  :  df['position'].max()+1,
        "xticks":  np.arange(df['position'].min(), df['position'].max()+1 , 1),
        "ymin"  :  -300.0,
        "ymax"  :  +300.0,
        "yticks":  np.arange(-250.0, 250.0+1.0, 50.0)    
        }

    bruce_magic_axis(ax=ax, **axes_vitals)
    ax.spines['bottom'].set_position(("data", axes_vitals["ymin"]))

    ticks =  ax.get_yticks()
    ax.set_yticklabels([int(abs(tick)) for tick in ticks])


    #stylise the frame
    ax.set_facecolor("None")
    ax.spines['left'].set_linewidth(0.75)
    ax.spines['bottom'].set_linewidth(0.75)
    ax.spines['right'].set_visible(0.75)
    ax.spines['top'].set_visible(0.75)

    y_labels = ax.get_yticklabels()
    plt.setp(y_labels, rotation=0, fontsize=label_size)


    fig.savefig(saveas, bbox_inches="tight", transparent=True)



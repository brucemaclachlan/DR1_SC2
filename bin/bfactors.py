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
            
            B factor functions
            
"""

def extract_b(selection):
    import pandas as pd
    import pymol

    myspace = {
                'position' : [],
                'residue'     : [],
                'atom'    : [],
                'bfactors' : []

                }

    tempfactor = 0
    atomnumber = 0
    pymol.cmd.iterate(selection, "bfactors.append(b)", space=myspace)
    pymol.cmd.iterate(selection, "position.append(resi)", space=myspace)
    pymol.cmd.iterate(selection, "atom.append(name)", space=myspace)
    pymol.cmd.iterate(selection, "residue.append(oneletter)", space=myspace)

    output = {}

    output['position'] = list(map(int, myspace['position']))
    output['residue'] = myspace['residue']
    output['atom'] = myspace['atom']
    output['bfactors'] = list(map(int, myspace['bfactors']))

    return pd.DataFrame.from_dict(output)

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


def plot_b_factors(df, colour, saveas):

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

    # B factors

    scatters = ax.scatter(x=df['position'], y=df['bfactors'], s=20.0, c=colour, alpha=0.5, edgecolors='none')

    ticks =  ax.get_xticks()
    ax.set_xticklabels("")

    # for i in range(len(df['ASA'])):
    #     print(df['ASA'][i])
    #     ax.text(x=df['position'][i], y=df['ASA'][i], s=df['ASA'][i], size=label_size-1, va="center", ha="center", color="k", weight="bold", fontfamily="monospaced")


    ax.text(x=df['position'].min()-1, y=-12.0, s="", fontsize=label_size, va="center", ha="center", color="k")
    ax.text(x=df['position'].max()+1, y=-12.0, s="", fontsize=label_size, va="center", ha="center", color="k")

    for index, row in df.iterrows():
        ax.text(x=row['position'], y=-30.0, s=row['residue'], fontsize=label_size, va="center", ha="center", color="k") 

    boxed_label(ax, "Peptide B-factors (backbone atoms)", x_left=0.0, x_right=1.0, y_center=1.15, pad=0.005, height=0.30, facecolor=colour, alpha=0.8, fontsize=label_size)
    ax.text(x=-0.19, y=0.5, s="B factor "+r"$(\mathrm{\AA}^2)$", rotation=90.0, transform=ax.transAxes, fontsize=label_size, va="center", ha="center")


    axes_vitals = {
        "aspect" :  (3.0,1.0),
        "xmin"  :  df['position'].min()-1,
        "xmax"  :  df['position'].max()+1,
        "xticks":  np.arange(df['position'].min(), df['position'].max()+1 , 1),
        "ymin"  :  0.0,
        "ymax"  :  180.0,
        "yticks":  np.arange(0.0, 180.0+1.0, 30.0)    
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
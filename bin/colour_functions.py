#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written using Sublime Text & Spyder

@author:    Bruce J MacLachlan, Division of Infection & Immunity, Cardiff University
            ORCID ID: 0000-0002-2685-2733
            Github: https://github.com/brucemaclachlan
            email: maclachlanb@cardiff.ac.uk
            
            Colour functions
"""

def hex_to_fraction(hex_str):
    '''
    pymol deals with colours as rgb values in fractions from 0-1.0
    Use to convert hex values  to fraction format.
    i.e.
 
    "#FFFFFF" to [1.0, 1.0, 1.0] e.g. white
    '''
    
    if hex_str.startswith('#'):
        hex_str = hex_str[1:]
    rgb255list = tuple([int(hex_str[i:i + 2], 16) for i in range(0, len(hex_str), 2)])
    output = []
    for v in rgb255list:
        output.append(v/255.0)
    print(output)
    return output
 
 
def set_new_colour(name, rgblist):
    '''
    Adds a new colour to the pymol colourset list.
    '''
    import pymol
    pymol.cmd.set_color(name, rgblist, mode=0, quiet=1)
    return None

def show_palette(palette):
    '''
    Makes a small plot in matplotlib and outputs png file of the colour palette requested from palettable. It pops out in the root directory.
    Select your colours from bottom [0] to top [-1]
    '''
    import matplotlib.pyplot as plt
    print("Selecting from palettable palette:"  )
    print(palette.name)
    print("Number of colours in palette:")
    print(palette.number)
    fig = plt.figure(frameon=True)
 
    ax = fig.add_axes([0,0, palette.number, 1])
    ax.axis('on')
 
    i = 0
    while i < palette.number:
        ax.axhspan(i, i+1.0, facecolor=palette.mpl_colors[i], zorder=1)
        ax.annotate(str(i), xy=(0.5, float(i)+0.5), zorder=10)
        i+=1
 
    print("Saving a colour swatch of colours in palette. See /"+palette.name+".png ..")
    fig.savefig(palette.name+".png")
 
 
    return None


def colour_hsv(base_colour, delta_value):
    import colorsys
    output = []
    if type(base_colour) == str:
        
        return output
        
    if type(base_colour) == list:
        
        modulate = (0.0, 0.0, delta_value)
        
        colour_in = colorsys.rgb_to_hsv(*base_colour)
        
        print(colour_in)
        
        colour_out = [x + y for x, y in zip(colour_in, modulate)]
        print("new value in HSV", colour_out)
        print("new value in rgb", colorsys.hsv_to_rgb(*colour_out))
        return list(colorsys.hsv_to_rgb(*colour_out))
    
def rgb255_to_fraction(rgb255list):  
    '''
    pymol deals with colours as rgb values in fractions from 0-1.0
    Use to convert rgb values in 255 format to fraction format.
    i.e.
    
    [255,255,255] to [1.0, 1.0, 1.0] e.g. white
    '''
    output = []
    for v in rgb255list:
        output.append(v/255.0)
    return output
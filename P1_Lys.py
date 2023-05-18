# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written using Sublime Text & Spyder

@author:    Bruce J MacLachlan, Division of Infection & Immunity, Cardiff University
            ORCID ID: 0000-0002-2685-2733
            Github: https://github.com/brucemaclachlan
            email: maclachlanb@cardiff.ac.uk
            
            T
"""

description = "\
This script uses PyMOL to make figures relating to a peptide-HLA structure\
I will use this to build and collate all analyses of a peptide-HLA structure to try and automate \
It will be very messy as there will be lots of cut, paste and force, for now.. \
"

peptide_colours = { 
                    "PKY_15mer"      :   "#DDDEDF"   ,
                    "SC2_1"          :   "#65C6BF"   ,
                    "SC2_2"          :   "#65C6BF"   ,
                    "SC2_3"          :   "#CBE4B4"   ,
                    "SC2_4"          :   "#65C6BF"   ,
                    "SC2_5"          :   "#65C6BF"   ,
                    "SC2_6"          :   "#65C6BF"   ,
                    "SC2_7"          :   "#65C6BF"   ,
                    "SC2_8"          :   "#65C6BF"   ,
                    "SC2_9"          :   "#2A9FC1"   ,
                    "SC2_10"         :   "#65C6BF"   ,
                    "SC2_11"         :   "#65C6BF"   ,
                    "SC2_12"         :   "#65C6BF"   ,
                    "SC2_13"         :   "#65C6BF"   ,
                    "SC2_14"         :   "#65C6BF"   ,
                    "SC2_15"         :   "#1884BD"   ,
                    "SC2_16"         :   "#1884BD"   ,
                    "SC2_17"         :   "#1884BD"   ,
                    "SC2_18"         :   "#214C9C"   ,
                    "SC2_19"         :   "#1884BD"   ,
                    "SC2_20"         :   "#1884BD"   ,
                    "SC2_21"         :   "#37B5C3"   ,
                    "SC2_22"         :   "#EEF1B7"   ,
                    "SC2_23"         :   "#52BDC2"   ,
                    "SC2_24"         :   "#32A9C2"   ,
                    "SC2_25"         :   "#65C6BF"   ,
                    "SC2_26"         :   "#1C66AC"   ,
                    "SC2_27"         :   "#214C9C"   ,
                    "X2_A2"          :   "#65C6BF"   ,
                    "X2_A7"          :   "#65C6BF"   ,
                    "X2_A2_omicron"  :   "#B73979"   ,
                    "SC2_6_omicron"  :   "#B73979"   ,
                    "PKY"            :   "#DDDEDF"   ,
                    "PKY_11R"        :   "#DDDEDF"

                   }

import os
import sys
import time
import subprocess
import pymol
import itertools
import argparse
import colorsys
import shutil

# import bin.peptide_MHC_visualise as pMHC_omit

import bin.data.viewSet as viewSet
import bin.data.colourSet as colourSet

from PIL import Image, ImageDraw

def parse_args():
    """parse arguments for command line"""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--ray', dest = 'do_ray',   action='store_true', required = False, help ='Do you want to render and save images? Add flag to generate images (time consuming).')
    parser.add_argument('--no_ray', dest = 'do_ray',   action='store_false', required = False, help ='Do you want to render and save images? Add flag to generate no images (quicker).')
    parser.add_argument('--view', dest = 'do_view',   action='store_true', required = False, help ='Do you want to view the pymol session on end? Add flag to leave open at end.')
    parser.add_argument('--mhc_class', dest = 'mhc_class',   type=str, required = True, help ='What is the MHC class of the structure (options: I or II')
    parser.add_argument('--peptide', dest = 'peptide',   type=str, required = True, help ='What is the name of the peptide')
    parser.add_argument('--chains', dest = 'chains',   type=str, required = False, help ='What are the chains of the MHC complex, i.e. A = DR1a, B = DR1b, C = peptide. ABC as here is default.')

    parser.add_argument('--pdb', dest = 'pdb',   type=str, required = True, help ='Input .pdb file')
    parser.add_argument('--mtz', dest = 'mtz',   type=str, required = False, help ='Input .mtz file')



    parser.set_defaults(do_apbs=False)
    parser.set_defaults(do_ray=False)
    parser.set_defaults(do_view=False)
    parser.set_defaults(do_bsa=False)
    parser.set_defaults(chains="ABC")
    
    args = parser.parse_args()
    return args

def initialisePymol():
    '''
    This function asks python to start a new pymol session and apply a set of parameters related pymol renders the molecules.
    i.e. I don't like shadows, so they are turned off.
    This helps to keep all figures consistent.
    '''
    print("\nInitialising pymol...\n")
    pymol.finish_launching(['pymol', '-c'])
    pymol.cmd.reinitialize()
    # set PyMOL parameters
    # these are my favourite parameters for making figures
    pymol.cmd.set("ray_shadows","0")
    pymol.cmd.set("specular", "off")
    pymol.cmd.set("orthoscopic", "on")
    pymol.cmd.bg_color("white")
    pymol.cmd.set("valence", 0)
    pymol.cmd.set("ray_opaque_background", "0")
    pymol.cmd.set("ray_trace_mode", 1)
    pymol.cmd.set("transparency_mode", 2)
    pymol.cmd.set("dash_round_ends", 0)
    return None

#### Some colour functions #####
 
 
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
    Bit obsolete because it ended up being a one-liner.
    '''
    pymol.cmd.set_color(name, rgblist, mode=0, quiet=1)
    return None
 
 
def show_palette(palette):
    '''
    Makes a small plot in matplotlib and outputs png file of the colour palette requested from palettable. It pops out in the root directory.
    Select your colours from bottom [0] to top [-1]
    '''
 
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

def set_new_colour(name, rgblist):
    pymol.cmd.set_color(name, rgblist, mode=0, quiet=0)
    return None

def wait4ray(query):  
    counter = 0
    spinner = itertools.cycle(['-', '/', '|', '\\'])
    while not os.path.exists(query):
        toWrite=next(spinner) + " Time elapsed: "+str(counter)+" seconds"
        sys.stdout.write(toWrite)  # write the next character
        sys.stdout.flush()                # flush stdout buffer (actual character display)
        sys.stdout.write(len(toWrite)*"\b")           # erase the last written chars
        time.sleep(1)
        counter += 1
    return None



def rayTime(saveas, do_ray):
    if do_ray == 0:
        return None
    else:
        print("Outputting image.. This may take a few seconds..")
        if os.path.exists(saveas):
            print("Removing "+saveas+" as it already exists!")
            os.remove(saveas)
        time.sleep(2)
        pymol.cmd.png(saveas,ray=do_ray,width=3000,height=3000, dpi=200)
        print("after cmd.png")
        wait4ray(saveas) 
        print("Done! "+str(saveas)+ " was outputted" )
        return None


def magickray(layers, saveas, view, direc, do_ray):
    
    if do_ray == 1:
        
        pymol_dir = direc
        
        print("~~~~~~~~~~~~~~~~~~~~~~~~ Generating "+saveas+" using magickray ~~~~~~~~~~~~~~~~~~~~~~~~")
        pymol.cmd.hide("everything", "all")
        new_dir = pymol_dir+"/"+saveas
        if os.path.exists(new_dir):
            print("Removing directory "+saveas+" and it's contents as it already exists!")
            shutil.rmtree(new_dir)
        os.mkdir(new_dir)
    
        i = 1
        for layer in layers:
            pymol.cmd.hide("everything", "all")
            for ob in layer:
                for representation in layer[ob]:
                    pymol.cmd.show(representation, ob)
            
            pymol.cmd.set_view(view)
            time.sleep(5)
            print("Generating layer "+ str(i)+ " as /"+pymol_dir+"/"+saveas + "/"+saveas+str(i)+".png")
            rayTime(pymol_dir+"/"+saveas+"/"+saveas+str(i)+".png", do_ray)
            i+=1
        print("Merging sub-images of "+saveas+" to "+saveas+" using subprocess.call to ImageMagick")
        command = " ".join(["convert", pymol_dir+"/"+saveas+"/"+saveas+"*", "-background", "none", "-flatten", pymol_dir+"/"+saveas+".png"])
        print(command)
        subprocess.call([command], shell=True)
        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Done! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        
    pymol.cmd.hide("everything", "all")
    
    for layer in layers:
        for ob in layer:
            for representation in layer[ob]:
                pymol.cmd.show(representation, ob)
    
    pymol.cmd.set_view(view)
    pymol.cmd.scene(key=saveas, action="store")
    return None   


def magick_outline(layer, view, image_name, subimage_name, direc, do_ray):
    '''
    Hacky way to generate an outline of an object then overlay it on the image stack.
    subimage_name should be image_name+i (can only do layers on top).
    This code isn't very flexible right now.
    '''
    if do_ray == True:
        
        pymol_dir = direc
        
        for lay in layer:
            selection = lay
            representation = layer[lay][0]
        
            pymol.cmd.hide("everything", "all")

            print(selection)
            pymol.cmd.create(selection.replace(" ", "_")+"o",selection)
            pymol.cmd.show(representation,selection+"o")

        pymol.cmd.set_view(view)
        
        rayTime(pymol_dir+"/"+image_name+"/"+subimage_name+".png", do_ray)
        time.sleep(2)

        command = " ".join(["convert", pymol_dir+"/"+image_name+"/"+subimage_name+".png", "-alpha", "extract", "-threshold", "0", "-negate", "-transparent", "white", pymol_dir+"/"+image_name+"/"+subimage_name+"_1.png"])
        print(command)
        subprocess.call([command], shell=True)
        command = " ".join(["convert", pymol_dir+"/"+image_name+"/"+subimage_name+"_1.png", "-negate", "-threshold", "1", "-edge", "10", pymol_dir+"/"+image_name+"/"+subimage_name+"_2.png"])
        print(command)
        subprocess.call([command], shell=True)
        command = " ".join(["convert", pymol_dir+"/"+image_name+"/"+subimage_name+"_2.png", "-negate", pymol_dir+"/"+image_name+"/"+subimage_name+"_3.png"])
        print(command)
        subprocess.call([command], shell=True)
        command = " ".join(["convert", pymol_dir+"/"+image_name+"/"+subimage_name+"_3.png", "-fuzz", "20%", "-transparent", "white", pymol_dir+"/"+image_name+"/"+subimage_name+"_4.png"])
        print(command)
        subprocess.call([command], shell=True)

        os.remove(pymol_dir+"/"+image_name+"/"+subimage_name+".png")
        os.remove(pymol_dir+"/"+image_name+"/"+subimage_name+"_1.png")
        os.remove(pymol_dir+"/"+image_name+"/"+subimage_name+"_2.png")
        os.remove(pymol_dir+"/"+image_name+"/"+subimage_name+"_3.png")
        
        command = " ".join(["convert", pymol_dir+"/"+image_name+"/*.png", "-background", "none", "-flatten", pymol_dir+"/"+image_name+".png"])
        print(command)
        subprocess.call([command], shell=True)
        
    return None
##########################################################
# Some functions that are nice to remove imagemagick     #
def overlay_two_images(layer_1, layer_2):
    """
    A wrapper to open an image at a given path, load the files as PIL.Image and perform a composite of each
    To be used in conjunction with recursive_alpha_composite

    Args:
        layer_1: Path to image of the first layer (lower)
        layer_2: Path to image of the second layer (lower)

    Saves:
        None
    Returns:
        A PIL Image of the output composite two images.

    """
    from PIL import Image

    im1 = Image.open(layer_1)
    im2 = Image.open(layer_2)

    return Image.alpha_composite(im1 , im2)


def recursive_alpha_composite(layers, im_size, image_output_name):
    """
    A recursive function that allows the alpha composition stack of multiple image layers. Utilises PIL Image.alpha_composite

    Requires:
        def overlay_two_images(layer_1, layer_2)
    Args:
        layers: A list of paths to images which are to be layered atop each other. First layer in list at bottom. Last layer in list at top.
        im_size: A touple containing the image size of the canvas which all images are pasted on. At present all images must be the same size i.e. (3000,3000)
        image_output_name: A file name path where the output image will be saved.

    Saves:
        A composite image of all layers pasted together. The first image in the layers list will be rendered at the bottom, the last layer rendered on the top.
    Returns:
        A PIL Image of the output composite image.

    """
    from PIL import Image

    print("\n", "Using PIL Image.alpha_composite to merge the images..")
    for layer in layers:
        print(layer)

    canvas = Image.new('RGBA', im_size).save(image_output_name)

    for layer in layers:
        canvas = overlay_two_images(image_output_name, layer).save(image_output_name)

    print("Done!")
    print("Composite image saved as" , image_output_name, "\n")
    return canvas

def image_single_colour_filter(image, input_colour, output_colour):

    """
    A single colour rgb pass filter. Takes an image and for every pixel which has that EXACT rgb value, the colour is passed. 
    All other pixels are changed to white fully transparent (255,255,255,0)
    There is also the ability to change the pass colour to a different output colour

    Args:
        image: a path to the image to be analysed
        input_colour: tuple of the rgb colour to be let through in the filter i.e. (200,200,200)
        output_colour: tuple of the rgb colour that all passed pixels will be converted to. Make this the same as input_colour if you want no change

    Saves:
        None
    Returns:
        A PIL Image of the output colour filtered image.

    """
    im = Image.open(image)
    im = im.convert("RGBA")
    datas = im.getdata()

    newData = []
    for pixel in datas:
        if pixel[0] == input_colour[0] and pixel[1] == input_colour[1] and pixel[2] == input_colour[2]:
            alpha = (pixel[3],)
            ca = output_colour + alpha
            newData.append(ca)
        else:
            newData.append((255,255,255,0))

    im.putdata(newData)
    return im


def image_not_transparent_filter(image, output_colour):

    im = Image.open(image)
    im = im.convert("RGBA")
    datas = im.getdata()

    newData = []
    for pixel in datas:

        if pixel[3] != 0.0:
            alpha = (pixel[3],)
            ca = output_colour + alpha
            newData.append(ca)
        else:
            newData.append((255,255,255,0))

    im.putdata(newData)
    return im

##########################################################


def pMHC_load(structure, file_path):
    pymol.cmd.load(file_path, structure)
    # Let's remove hydrogens from all structures (if they have them) otherwise they can mess with alignments even if we don't see them
    pymol.cmd.remove("hydro")
    return None

def pMHC_align(mobile, target, mhc_class, chains_dict):
    # This aligns the a1 and a2 helices of our moving (mobile) and reference (target) MHC molecules

    #class I
    if mhc_class == "I":
        pymol.cmd.align(mobile+" and chain "+chains_dict["MHCA"]+" and resi 137-181,49-85", target+" and chain A and resi 137-181,49-85")

    #class II
    if mhc_class == "II":
        pymol.cmd.align(mobile+" and chain "+chains_dict["MHCA"]+" and resi 46-78 or "+mobile+" and chain "+chains_dict["MHCB"]+" and resi 54-91", target+" and chain A and resi 46-78 or "+target+" and chain B and resi 54-91")

    return None
        
def pMHC_create_objects(structure, mhc_class, chains_dict):


    if mhc_class == 'I':

        pymol.cmd.select(structure+"_MHCa1a2", selection=structure+" and chain "+chains_dict["MHCA"]+" and resi 137-181,49-85")
        pymol.cmd.create(structure+"_MHCa1a2_obj", selection=structure+"_MHCa1a2")

        pymol.cmd.select(structure+"_MHCa1", selection=structure+" and chain "+chains_dict["MHCA"]+" and resi 49-85")
        pymol.cmd.create(structure+"_MHCa1_obj", selection=structure+"_MHCa1")

        pymol.cmd.select(structure+"_MHCa2", selection=structure+" and chain "+chains_dict["MHCA"]+" and resi 137-181")
        pymol.cmd.create(structure+"_MHCa2_obj", selection=structure+"_MHCa2")
        
        pymol.cmd.select(structure+"_MHCgroove", selection=structure+" and chain "+chains_dict["MHCA"]+" and resi 1-136")
        pymol.cmd.create(structure+"_MHCgroove_obj", selection=structure+"_MHCgroove")

        pymol.cmd.select(structure+"_p", selection=structure+" and chain "+chains_dict["PEPTIDE"]+"")
        pymol.cmd.create(structure+"_p_obj", selection=structure+"_p")

        pymol.cmd.select(structure+"_pMHCgroove", selection=structure+" and chain "+chains_dict["MHCA"]+" and resi 1-136 or "+structure+" and chain "+chains_dict["PEPTIDE"]+"")
        pymol.cmd.create(structure+"_pMHCgroove_obj", selection=structure+"_pMHCgroove")    
      

    if mhc_class == 'II':

        pymol.cmd.select(structure+"_MHCa1a2", selection=structure+" and chain "+chains_dict["MHCA"]+" and resi 46-78 or "+structure+" and chain "+chains_dict["MHCB"]+" and resi 54-91")
        pymol.cmd.create(structure+"_MHCa1a2_obj", selection=structure+"_MHCa1a2")

        pymol.cmd.select(structure+"_MHCa1", selection=structure+" and chain "+chains_dict["MHCA"]+" and resi 46-78")
        pymol.cmd.create(structure+"_MHCa1_obj", selection=structure+"_MHCa1")

        pymol.cmd.select(structure+"_MHCa2", selection=structure+" and chain "+chains_dict["MHCB"]+" and resi 54-91")
        pymol.cmd.create(structure+"_MHCa2_obj", selection=structure+"_MHCa2")
        
        pymol.cmd.select(structure+"_MHCgroove", selection=structure+" and chain "+chains_dict["MHCA"]+" and resi 1-78 or "+structure+" and chain "+chains_dict["MHCB"]+" and resi 6-91")
        pymol.cmd.create(structure+"_MHCgroove_obj", selection=structure+"_MHCgroove")

        pymol.cmd.select(structure+"_p", selection=structure+" and chain "+chains_dict["PEPTIDE"]+"")
        pymol.cmd.create(structure+"_p_obj", selection=structure+"_p")
      
        pymol.cmd.select(structure+"_pMHCgroove", selection=structure+" and chain "+chains_dict["MHCA"]+" and resi 1-78 or "+structure+" and chain "+chains_dict["MHCB"]+" and resi 6-91 or "+structure+" and chain "+chains_dict["PEPTIDE"]+"")
        pymol.cmd.create(structure+"_pMHCgroove_obj", selection=structure+"_pMHCgroove")        



    pymol.cmd.delete("hetatms")

    return None


    ################################ BODY ########################################
# This is the "start" of the script. This should be main.


################################
# Housekeeping & Startup       #
################################

# Handle commoand line arguments
args = parse_args()
print("Running analysis with the following inputs.. ")
print(args)
ray = args.do_ray
view = args.do_view
do_apbs = args.do_apbs
peptide_name = args.peptide
do_bsa = args.do_bsa
chains = args.chains

mhc_class = args.mhc_class

pdb = args.pdb
pdb_name = pdb.split(".")[0]
pdb = os.path.realpath(pdb)

mtz = args.mtz
mtz_name = mtz.split(".")[0]
mtz = os.path.realpath(mtz)

bin_dir = os.getcwd()+"/bin"
home_dir = os.getcwd()

output_dir = os.getcwd()+"/"+pdb_name+"_"+chains

if os.path.exists(output_dir) == False:
    os.mkdir(output_dir)

if os.path.exists(output_dir+"/visualisation") == False:
    os.mkdir(output_dir+"/visualisation")

# need to add this line as pymol.png wants 0 or 1 not bool
if ray == True:
    do_ray = 1
if ray == False:
    do_ray = 0

if mhc_class == 'I':
    template = 'bin/data/I_cdr_template.pdb'
if mhc_class == 'II':
    template = 'bin/data/II_cdr_template.pdb'

chains_dict = { "MHCA" : chains[0],
                "MHCB" : chains[1],
                "PEPTIDE" : chains[2]
}

########################
# Visualisations       #
########################

# This creates a new session of PyMOL and runs my favourite viewing parameters for PyMOL
# These parameters can be seen in def intialisePymol() function above
initialisePymol()

structure = pdb_name

pMHC_load(structure, pdb)

pymol.cmd.load(template, "align_target")
pMHC_align(structure, "align_target", mhc_class, chains_dict)
# pymol.cmd.delete('align_target')

chains_dict["MHCA"]
pymol.cmd.remove(structure+ " and NOT chain "+chains_dict["MHCA"]+","+chains_dict["MHCB"]+","+chains_dict["PEPTIDE"])
pymol.cmd.remove("hetatm and NOT resn DHL")

pMHC_create_objects(structure, mhc_class, chains_dict)


peptide_colour = hex_to_fraction(peptide_colours[peptide_name])
set_new_colour(peptide_name, peptide_colour)

    
# Set colours
pymol.cmd.color("gray90",structure+" and chain "+chains_dict["MHCA"])

pymol.cmd.color("gray90",structure+"_MHCa1_obj")
pymol.cmd.color("gray90",structure+"_MHCa2_obj")
pymol.cmd.color("gray90",structure+"_MHCgroove_obj")

pymol.cmd.color("gray60",structure+" and chain "+chains_dict["MHCB"])

print(peptide_colour , structure+" and chain "+chains_dict["PEPTIDE"])
pymol.cmd.color(peptide_name,structure+" and chain "+chains_dict["PEPTIDE"])


pymol.cmd.util.cnc(structure+" and chain "+chains_dict["PEPTIDE"])   

#### VIEWS ####

#### IMAGE 1 ####
# overview of all structure

pymol.cmd.hide("everything", "all")

pymol.cmd.show("cartoon",structure+" and chain "+chains_dict["MHCA"]+" or "+structure+" and chain "+chains_dict["MHCB"])
pymol.cmd.show("sticks",structure+" and chain "+chains_dict["PEPTIDE"]) 


if peptide_name in ["X2_A2", "X2_A2_omicron"]:

    
    pymol.cmd.unbond(structure+" and resn DHL", structure+" and chain "+chains_dict["PEPTIDE"]+" and resi 3")
    pymol.cmd.bond(structure+" and resn DHL and name SG", structure+" and chain "+chains_dict["PEPTIDE"]+" and resi 3 and name SG" )

    pymol.cmd.set("stick_transparency" , 0.5 , structure+" and resn DHL")

    pymol.cmd.select("SGs", structure+" and chain C and name SG")
    pymol.cmd.set("stick_transparency" , 0.5 , "SGs")

pymol.cmd.set_view(viewSet.pMHC_1)
pymol.cmd.scene(structure+"_pMHC", "store")


########################
# Clip cut-away        #
########################

pymol.cmd.hide("everything" , "all")

# must disable depth cue and shadows
pymol.cmd.set("depth_cue" , 0)

# this controls the z depth of the slice plane
# (sets it halfway between the clipping planes)
fraction = 0.39

# pymol.cmd.center(structure+"_p_obj")
# pymol.cmd.turn("x", -90)

pymol.cmd.set_view("\
    -0.782565236,    0.267936289,   -0.561967373,\
     0.592841208,    0.596321344,   -0.541243196,\
     0.190093905,   -0.756709278,   -0.625500619,\
    -0.000302987,    0.000438260, -193.702819824,\
    64.388168335,   28.303731918,   36.209182739,\
   160.413558960,  259.850311279,   20.000000000" )


# view = pymol.cmd.get_view()
# near_dist = fraction*(view[16]-view[15])
# far_dist = (view[16]-view[15]) - near_dist
# pymol.cmd.clip("near", -near_dist)
# pymol.cmd.center(structure+"_p_obj")

# pymol.cmd.show("sticks" , structure+"_p_obj")
pymol.cmd.show("surface" , structure+"_MHCgroove_obj")

# pymol.cmd.center(structure+"_p_obj")

pymol.cmd.set("ray_interior_color", "grey80")
pymol.cmd.set("surface_color", "white")
pymol.cmd.scene(structure+"_cutaway", "store")
rayTime(output_dir+"/visualisation/"+ structure+"_cutaway_back.png", do_ray)
# pymol.cmd.hide("sticks" , "all")

# Make a coverer to hide that annoying hole in the clipped surface
image = Image.new('RGBA', (3000, 3000))
draw = ImageDraw.Draw(image)
draw.rectangle((1100, 1300, 1300, 1600), fill = (155,155,155), outline =None)
image.save(output_dir+"/visualisation/"+ structure+'_coverer.png')


layer_1 = output_dir+"/visualisation/"+ structure+"_cutaway_back.png"
layer_2 = output_dir+"/visualisation/"+ structure+"_coverer.png"
layers = [layer_1, layer_2]
overlay_output_name = output_dir+"/visualisation/"+ structure+"_cutaway_back_edit.png"

cutaway_back_edit = recursive_alpha_composite(layers, (3000,3000), overlay_output_name)


######################################################################################################################
# Make a version of the above but with a greenscreen cutaway colour.. this is for the flat/goodsell image used later #
######################################################################################################################

greenscreen_color = "0fa110"
greenscreen_color_srgb = 'sRGB(15,161,16)'
greenscreen_color_rgb = (15, 161, 16)

set_new_colour("greenscreen", list((15, 161, 16)))
pymol.cmd.set("ray_interior_color" , "greenscreen")
rayTime(output_dir+"/visualisation/"+ structure+"_cutaway_back_green.png", do_ray)

# The greenscreen_color I chose isn't coming out exactly that rgb value when rayed through pymol
# I can't work out why so for now, let's just use colour picker to re-assign the variable to what comes out in reality
# This doesn't matter anyway as the green colour is just arbitrary


greenscreen_color = "0b790c"
greenscreen_color_srgb = 'sRGB(11,121,12)'
greenscreen_color_rgb = (11,121,12)


image = Image.new('RGBA', (3000, 3000))
draw = ImageDraw.Draw(image)
draw.rectangle((1100, 1300, 1300, 1600), fill = greenscreen_color_rgb, outline =None)
image.save(output_dir+"/visualisation/"+ structure+'_coverer_green.png')


layer_1 = output_dir+"/visualisation/"+ structure+"_cutaway_back_green.png"
layer_2 = output_dir+"/visualisation/"+ structure+"_coverer_green.png"
layers = [layer_1, layer_2]
overlay_output_name = output_dir+"/visualisation/"+ structure+"_cutaway_back_green_edit.png"

cutaway_back_green_edit = recursive_alpha_composite(layers, (3000,3000), overlay_output_name)

######################################################################################################################


pymol.cmd.hide("everything" , "all")

pymol.cmd.set_view("\
    -0.782565236,    0.267936289,   -0.561967373,\
     0.592841208,    0.596321344,   -0.541243196,\
     0.190093905,   -0.756709278,   -0.625500619,\
    -0.000302098,    0.000437962, -193.702819824,\
    64.444366455,   28.357858658,   36.271736145,\
  -5579.197753906, 5999.861328125,   20.000000000" )

pymol.cmd.show("sticks" , structure+"_p_obj")
pymol.cmd.color(peptide_name, structure+"_p_obj")
pymol.cmd.util.cnc(structure+"_p_obj")

rayTime(output_dir+"/visualisation/"+ structure+"_cutaway_front.png", do_ray)


layer_1 = output_dir+"/visualisation/"+ structure+"_cutaway_back_edit.png"
layer_2 = output_dir+"/visualisation/"+ structure+"_cutaway_front.png"
layers = [layer_1, layer_2]
overlay_output_name = output_dir+"/visualisation/"+ structure+"_cutaway.png"

cutaway = recursive_alpha_composite(layers, (3000,3000), overlay_output_name)

# pymol.cmd.select("near_copy1_P1_Lys", "byres chain A within 8.0 of (chain C and resi 3) or byres chain B within 8.0 of (chain C and resi 3) ")
# pymol.cmd.show("sticks", "near_copy1_P1_Lys")
# pymol.cmd.util.cnc("near_copy1_P1_Lys")


#################################
# Clip cut-away P1 focus        #
#################################

pymol.cmd.set("ray_interior_color", "grey80")

pymol.cmd.hide("sticks" , "all")
pymol.cmd.show("surface" , structure+"_MHCgroove_obj")
# pymol.cmd.set("transparency", 0.5, structure+"_MHCgroove_obj")

pymol.cmd.set_view("\
    -0.782565236,    0.267936289,   -0.561967373,\
     0.592841208,    0.596321344,   -0.541243196,\
     0.190093905,   -0.756709278,   -0.625500619,\
    -0.000251163,    0.000273665,  -88.982124329,\
    74.280212402,   22.479646683,   32.407936096,\
    55.663780212,  155.100418091,   20.000000000")
### cut above here and paste into script ###)

pymol.cmd.scene(structure+"_pMHC_P1_Lys_focus_surface", "store")

rayTime(output_dir+"/visualisation/"+ structure+"_P1_focus_cutaway_surface.png", do_ray)


image = Image.new('RGBA', (3000, 3000))
draw = ImageDraw.Draw(image)
draw.rectangle((1750, 1400, 2250, 1850), fill = (155,155,155), outline =None)
image.save(output_dir+"/visualisation/"+ structure+'zoom_coverer.png')


layer_1 = output_dir+"/visualisation/"+ structure+"_P1_focus_cutaway_surface.png"
layer_2 = output_dir+"/visualisation/"+ structure+"zoom_coverer.png"
layers = [layer_1, layer_2]
overlay_output_name = output_dir+"/visualisation/"+ structure+"_P1_focus_cutaway_surface_edit.png"

cutaway_zoom_surface_edit = recursive_alpha_composite(layers, (3000,3000), overlay_output_name)


pymol.cmd.set_view("\
    -0.782565236,    0.267936289,   -0.561967373,\
     0.592841208,    0.596321344,   -0.541243196,\
     0.190093905,   -0.756709278,   -0.625500619,\
    -0.000273535,    0.000260269,  -88.982124329,\
    71.976158142,   20.260526657,   29.843364716,\
   -97.736183167,  155.500396729,   20.000000000" )
### cut above here and paste into script ###)

pymol.cmd.hide("surface" , "all")
pymol.cmd.select("custom_near_P1", structure+" and chain "+chains_dict["MHCA"]+" and resi 7,24,25,26,31,32,43,48 or "+structure+" and chain "+chains_dict["MHCB"]+" and resi 85,89,90,153")

# MHCB residue 82 Asn is in the front of the pocket. I think it makes sense to not show it to make viewing the open pocket easier
# pymol.cmd.select("custom_near_P1", structure+" and chain "+chains_dict["MHCA"]+" and resi 7,24,25,26,31,32,43,48 or "+structure+" and chain "+chains_dict["MHCB"]+" and resi 82,85,89,90,153")

pymol.cmd.show("sticks", "custom_near_P1")
pymol.cmd.util.cnc("custom_near_P1")

pymol.cmd.set("stick_transparency" , 0.5 , structure+" and chain "+chains_dict["MHCA"]+" and resi 31,32,43")

# MHCB residue 82 Asn is in the front of the pocket. I think it makes sense to not show it to make viewing the open pocket easier
# pymol.cmd.set("stick_transparency" , 0.5 , structure+" and chain "+chains_dict["MHCB"]+" and resi 82")

# pymol.cmd.show("cartoon" , structure+"_MHCgroove_obj")
# pymol.cmd.set("cartoon_transparency" , 0.9 , structure+"_MHCgroove_obj")

pymol.cmd.scene(structure+"_pMHC_P1_Lys_focus_pocket", "store")
rayTime(output_dir+"/visualisation/"+ structure+"_P1_focus_cutaway_pocket.png", do_ray)

pymol.cmd.hide("sticks" , "all")

pymol.cmd.show("sticks" , structure+"_p_obj")

pymol.cmd.scene(structure+"_pMHC_P1_Lys_focus_peptide", "store")
rayTime(output_dir+"/visualisation/"+ structure+"_P1_focus_cutaway_peptide.png", do_ray)


layer_1 = output_dir+"/visualisation/"+ structure+"_P1_focus_cutaway_surface_edit.png"
layer_2 = output_dir+"/visualisation/"+ structure+"_P1_focus_cutaway_pocket.png"
layers = [layer_1, layer_2]
overlay_output_name = output_dir+"/visualisation/"+ structure+"_P1_focus_cutaway.png"

cutaway_zoom = recursive_alpha_composite(layers, (3000,3000), overlay_output_name)

########### Goodselly peptide ###########


cell_press_colours = {"cell_1" : { "hex" : '0F283D',
                                   "rgb" : (15,40,61)
                                 },
                      "cell_2" : { "hex" : '6C9EC8',
                                   "rgb" : (108,158,200)
                                 },
                      "cell_3" : { "hex" : '2E7DBC',
                                   "rgb" : (45,125,188)
                                 },
                      "cell_4" : { "hex" : '21313D',
                                   "rgb" : (33,49,61)
                                 },
                      "cell_5" : { "hex" : '215B8A',
                                   "rgb" : (33,91,138)
                                 }
                     }
    
pymol.cmd.hide("everything", "all")
# pymol.cmd.show("spheres" ,  structure+"_p_obj")
pymol.cmd.show("cartoon" ,  structure+"_p_obj")
pymol.cmd.show("sticks" ,  structure+"_p_obj and sc."+" or "+ structure+"_p_obj and name ca")
pymol.cmd.hide("everything", "hetatm")


pymol.cmd.set("cartoon_loop_radius",0.4)
pymol.cmd.set("stick_radius" , 0.4)


# qutemol_like()

pymol.cmd.set_view("\
    -0.782565236,    0.267936289,   -0.561967373,\
     0.592841208,    0.596321344,   -0.541243196,\
     0.190093905,   -0.756709278,   -0.625500619,\
    -0.000302098,    0.000437962, -193.702819824,\
    64.444366455,   28.357858658,   36.271736145,\
  -5579.197753906, 5999.861328125,   20.000000000" )


set_new_colour("cell_3", list(cell_press_colours["cell_3"]["rgb"]))
pymol.cmd.color("cell_3", structure+"_p_obj")
rayTime(output_dir+"/visualisation/"+ structure+"_peptide_spheres.png", do_ray)

groove_front = image_single_colour_filter(output_dir+"/visualisation/"+ structure+"_cutaway_back_green_edit.png" , greenscreen_color_rgb, cell_press_colours["cell_5"]["rgb"]).save(output_dir+"/visualisation/"+ structure+"_PILtoon1.png")
groove_all = image_not_transparent_filter(output_dir+"/visualisation/"+ structure+"_cutaway_back_green_edit.png" , cell_press_colours["cell_2"]["rgb"]).save(output_dir+"/visualisation/"+ structure+"_PILtoon2.png")
peptide_blob = image_not_transparent_filter(output_dir+"/visualisation/"+ structure+"_peptide_spheres.png", cell_press_colours["cell_3"]["rgb"]).save(output_dir+"/visualisation/"+ structure+"_PILtoon3.png")

layer_1 = output_dir+"/visualisation/"+ structure+"_PILtoon2.png"
layer_2 = output_dir+"/visualisation/"+ structure+"_PILtoon1.png"
layer_3 = output_dir+"/visualisation/"+ structure+"_PILtoon3.png"

overlay_output_name = output_dir+"/visualisation/"+ structure+"_PILtoon.png"

layers = [layer_1, layer_2, layer_3]

piltoon = recursive_alpha_composite(layers, (3000,3000), overlay_output_name)


if peptide_name in ["X2_A2_omicron", "SC2_6_omicron"]:
    pymol.cmd.hide("everything", "all")
    if peptide_name == "X2_A2_omicron":

        pymol.cmd.show("sticks" ,  structure+"_p_obj and resi 8,11,13 and sc."+" or "+ structure+"_p_obj and resi 8,11,13 and name ca")
        pymol.cmd.color(peptide_name, structure+"_p_obj")

    if peptide_name == "SC2_6_omicron":

        pymol.cmd.show("sticks" ,  structure+"_p_obj and sc."+" or "+ structure+"_p_obj and name ca")
        pymol.cmd.show("cartoon" , structure+"_p_obj")
        pymol.cmd.color(peptide_name, structure+"_p_obj")

    rayTime(output_dir+"/visualisation/"+ structure+"_peptide_spheres_omicron.png", do_ray)
    peptide_residue_blob = image_not_transparent_filter(output_dir+"/visualisation/"+ structure+"_peptide_spheres_omicron.png", (198,99,143)).save(output_dir+"/visualisation/"+ structure+"_PILtoon4.png")



    layer_4 = output_dir+"/visualisation/"+ structure+"_PILtoon4.png"
    layers = [layer_1, layer_2, layer_3, layer_4]
    overlay_output_name = output_dir+"/visualisation/"+ structure+"_PILtoon_omicron.png"

    piltoon_omicron = recursive_alpha_composite(layers, (3000,3000), overlay_output_name)

###########
# This is the end. We finally want to save the PyMOL session file and quit pymol. 
# If we have asked to view the results in PyMOL, we will finish by opening up the session file.
###########

pymol.cmd.save(output_dir+"/"+pdb_name+"_P1_lys.pse")


pymol.cmd.quit()
if view == True:
    subprocess.call(["pymol", output_dir+"/"+pdb_name+"_P1_lys.pse"])
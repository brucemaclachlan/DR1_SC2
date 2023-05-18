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

description = "\
This script uses PyMOL to make figures relating to a peptide-HLA structure\
I will use this to build and collate all analyses of a peptide-HLA structure to try and automate \
It will be very messy as there will be lots of cut, paste and force, for now.. \
"

peptide_colours = { 
                    "PKY_15mer"      :   "#DDDEDF"   ,
                    "SC2_1"          :   "#65C6BF"   ,
                    "SC2_2"          :   "#65C6BF"   ,
                    "SC2_3"          :   "#CBE484"   ,
                    "SC2_4"          :   "#65C6BF"   ,
                    "SC2_5"          :   "#65C6BF"   ,
                    "SC2_6"          :   "#65C6BF"   ,
                    "SC2_7"          :   "#65C6BF"   ,
                    "SC2_8"          :   "#65C6BF"   ,
                    "SC2_9"          :   "#B4D9B9"   ,
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
                    "PKY"            :   "#FFFF00"   ,
                    "PKY_11R"        :   "#FFFF00"

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

def parse_args():
    """parse arguments for command line"""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--ray', dest = 'do_ray',   action='store_true', required = False, help ='Do you want to render and save images? Add flag to generate images (time consuming).')
    parser.add_argument('--no_ray', dest = 'do_ray',   action='store_false', required = False, help ='Do you want to render and save images? Add flag to generate no images (quicker).')
    parser.add_argument('--view', dest = 'do_view',   action='store_true', required = False, help ='Do you want to view the pymol session on end? Add flag to leave open at end.')
    parser.add_argument('--do_apbs', dest = 'do_apbs',   action='store_true', required = False, help ='Do you want to run apbs? Add flag to do so if you need to generate electrostatic maps (time consuming).')
    parser.add_argument('--mhc_class', dest = 'mhc_class',   type=str, required = True, help ='What is the MHC class of the structure (options: I or II')
    parser.add_argument('--peptide', dest = 'peptide',   type=str, required = True, help ='What is the name of the peptide')
    parser.add_argument('--do_bsa', dest = 'do_bsa',   action='store_true', required = False, help ='Do you want to generate bsa analysis and plots (slower).')
    parser.add_argument('--chains', dest = 'chains',   type=str, required = False, help ='What are the chains of the MHC complex, i.e. A = DR1a, B = DR1b, C = peptide. ABC as here is default.')

    parser.add_argument('--pdb', dest = 'pdb',   type=str, required = True, help ='Input .pdb file')
    parser.add_argument('--mtz', dest = 'mtz',   type=str, required = True, help ='Input .mtz file')



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
    pymol.cmd.set_color(name, rgblist, mode=0, quiet=1)
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

def pMHC_omit_maps(pdb, mtz, mhc_class, chains, ray, output_dir, file_name, peptide_name):

    # Sort chains
    MHCachain, MHCbchain, peptidechain = chains_dict["MHCA"] , chains_dict["MHCB"], chains_dict["PEPTIDE"]
    pdb_name = file_name.lower()

    # Make output folder #
    print(file_name)

    if not os.path.exists(output_dir):
        print("Creating Directory " + output_dir)
        os.makedirs(output_dir)

    if not os.path.exists(output_dir + "/visualisation"):
        print("Creating Directory " + output_dir + "/visualisation")
        os.makedirs(output_dir + "/visualisation")

    if not os.path.exists(output_dir+ "/maps"):
        print("Creating Directory " + output_dir + "/maps")
        os.makedirs(output_dir + "/maps")
   
    if not os.path.exists(output_dir + "/electrostatics"):
        print("Creating Directory " + output_dir + "/electrostatics")
        os.makedirs(output_dir + "/electrostatics")

    if not os.path.exists(output_dir + "/pdbs"):
        print("Creating Directory " + output_dir + "/pdbs")
        os.makedirs(output_dir + "/pdbs")

    if not os.path.exists(output_dir + "/bsa"):
        print("Creating Directory " + output_dir + "/bsa")
        os.makedirs(output_dir + "/bsa")

    if mtz != "ebi":
        print("A map.mtz file was provided!", mtz, "will be moved to", output_dir + "/maps/" + file_name + ".mtz")

        shutil.copy(mtz, output_dir + "/maps/" + file_name + ".mtz")

    if mtz == "ebi":
        if not os.path.exists(output_dir + "/maps/" + file_name + ".mtz"):
            print("Downloading map.mtz for entry", pdb_name, "from the PDBe (EBI)")

            try:
                urllib.urlretrieve("http://www.ebi.ac.uk/pdbe/coordinates/files/" + pdb_name + "_map.mtz",
                               output_dir + "/maps/" + file_name + ".mtz")
            except:
                print("Could not retrieve url. Please try again, making sure you either supply a file," \
                      " or your file shares its name with one on PDB")
                # quit early, get rid of pymol
                pymol.cmd.quit()
                sys.exit()


        else:
            "Did not need to download from ebi as map.mtz already exists"

        if os.path.exists(output_dir + "/maps/" + file_name + ".mtz"):
            print("Download successful!")

    mtz = output_dir + "/maps/" + file_name + ".mtz"

    # COMMENT OUT FROM HERE TO SPEED UP #


    tempTxt = "rmchain " + peptidechain + "\n" + "END"
    temp = open(output_dir + "/pdbs/" + file_name + "_pdbcurPARAM.tmp", "w")
    temp.write(tempTxt)
    temp.close()

    # delete chain
    os.system("pdbcur XYZIN " + pdb + " XYZOUT " + output_dir + "/pdbs/" + file_name +
              "_nopeptide.pdb" + " < " + output_dir + "/pdbs/" + file_name + "_pdbcurPARAM.tmp")
    #os.remove(file_name + "/pdbs/" + file_name + "_pdbcurPARAM.tmp")


    #Run one cycle of refmac without peptide in the groove
    print("~"*40+"Running refmac without peptide"+"~"*40)
    os.system("refmac5 XYZIN " + output_dir + "/pdbs/" + file_name + "_nopeptide.pdb" + " XYZOUT " +
              output_dir + "/pdbs/" + file_name + "_refmac5_omitmap.pdb" + " HKLIN " + output_dir + "/maps/"
              + file_name + ".mtz" + " HKLOUT " + output_dir + "/maps/" + file_name + "_refmac5_omitmap.mtz" + " LIBOUT "
              + output_dir + "/maps/" + file_name + "_refmac_omitmap.cif" + " < " + "bin/data/refmacOMITparams.tmp")

    # Use CCP4 to generate a 2fo-1fc map of the refined mtz (F1=F-obs-filtered PHI=SigF-obs-filtered

    print("~"*40+"EDMparam1 below"+"~"*40)
    os.system("fft HKLIN " + mtz + " MAPOUT " + output_dir + "/maps/" +
              file_name + ".map1.tmp" + " < " + "bin/data/EDMparam1.tmp")
    print("~"*40+"EDMparam2 below"+"~"*40)
    os.system("mapmask MAPIN " + output_dir + "/maps/" + file_name + ".map1.tmp" +
              " MAPOUT " + output_dir + "/maps/" + file_name + ".map.ccp4" + " XYZIN " + pdb + " < " + "bin/data/EDMparam2.tmp")



    # Use CCP4 to generate fo-fc (difference) map of the omit refined mtz (i.e. from the above refmac round)

    print("~"*40+"EDMparam3 below"+"~"*40)
    os.system("fft HKLIN " + output_dir + "/maps/" + file_name + "_refmac5_omitmap.mtz"+ " MAPOUT " + output_dir + "/maps/" +
              file_name + ".map3.tmp" + " < " + "bin/data/EDMparam3.tmp")
    print("~"*40+"EDMparam4 below"+"~"*40)
    os.system("mapmask MAPIN " + output_dir + "/maps/" + file_name + ".map3.tmp" +
              " MAPOUT " + output_dir + "/maps/" + file_name + ".difference_map.ccp4" +
              " XYZIN " + pdb + " < " + "bin/data/EDMparam4.tmp")

    os.remove(output_dir + "/maps/" + file_name + ".map1.tmp")
    os.remove(output_dir + "/maps/" + file_name + ".map3.tmp")



    # COMMENT OUT UP TO HERE TO SPEED UP #

    edmap = output_dir + "/maps/" + file_name + ".map.ccp4"

    diffmap = output_dir + "/maps/" + file_name + ".difference_map.ccp4"

    # Find the MHC helices
    if mhc_class == "I":
        a1locs = list(range(50, 86))
        MHCa1 = ["MHCa"]
        MHCa1 = MHCa1 + a1locs
        a2locs = list(range(140, 176))
        MHCa2 = ["MHCa"]
        MHCa2 = MHCa2 + a2locs

    if mhc_class == "II":
        a1locs = list(range(46, 78))
        MHCa1 = ["MHCa"]
        MHCa1 = MHCa1 + a1locs
        a2locs = list(range(54, 91))
        MHCa2 = ["MHCb"]
        MHCa2 = MHCa2 + a2locs

    # Let's get started

    initialisePymol()
    pymol.cmd.reinitialize()
    initialisePymol()

    pymol.cmd.load(pdb, "structure")
    pymol.cmd.load(edmap, file_name + "_map")
    pymol.cmd.load(diffmap, file_name + "_dmap")

    pymol.cmd.remove("structure and NOT chain "+chains_dict["MHCA"]+","+chains_dict["MHCB"]+","+chains_dict["PEPTIDE"])
    pymol.cmd.save("whatsonhere.pse")
    pymol.cmd.remove("hetatm and NOT resn DHL")

    # align to template
    print("\nAligning file to template...\n")



    if mhc_class == 'I':
        template = 'bin/data/I_cdr_template.pdb'
    if mhc_class == 'II':
        template = 'bin/data/II_cdr_template.pdb'

    pymol.cmd.load(template, mhc_class+"_template")
    pymol.cmd.align("structure", mhc_class+"_template")
    pymol.cmd.matrix_copy("structure", file_name + "_map")
    pymol.cmd.matrix_copy("structure", file_name + "_dmap")
    pymol.cmd.delete(mhc_class + "_template")
    print("\nAlignment to " + mhc_class + "_cdr_template.pdb complete!\n")

    # Make chains objects
    pymol.cmd.select("MHCas", selection="chain " + MHCachain)
    pymol.cmd.select("MHCbs", selection="chain " + MHCbchain)
    pymol.cmd.select("ps", selection="chain " + peptidechain)


    peptide_colour = hex_to_fraction(peptide_colours[peptide_name])
    set_new_colour(peptide_name, peptide_colour)

    pymol.cmd.hide("all")
    pymol.cmd.show("sticks", "ps")
    pymol.cmd.color(peptide_name, "ps")
    pymol.cmd.util.cnc("ps")

    # Select the MHCa helices
    # MHCa1, MHCa2 = None, None

    locs = '+'.join(str(x) for x in MHCa1[1:])
    pymol.cmd.select("MHCa1", selection=MHCa1[0] + "s and resi " + locs)
    locs = '+'.join(str(x) for x in MHCa2[1:])
    pymol.cmd.select("MHCa2", selection=MHCa2[0] + "s and resi " + locs)

    # Make electron density map
    # pymol.cmd.map_double(file_name + "_map", -1)
    pymol.cmd.isomesh("p_map_1sigma", file_name + "_map", 1.0, "ps", carve=1.6)
    pymol.cmd.isomesh("p_map_05sigma", file_name + "_map", 0.5, "ps", carve=1.6)
    pymol.cmd.set("mesh_width", 0.5)
    pymol.cmd.set("mesh_radius" , 0.02)

    pymol.cmd.hide("mesh", "all")

    pymol.cmd.set_view(viewSet.pMHC_2)

    pymol.cmd.show("mesh", "p_map_1sigma")
    pymol.cmd.color("grey50", "p_map_1sigma")

    # we need to go to ray trace mode 0 as the mesh doesn't work with a ray_trace_gain (outline)
    pymol.cmd.set("ray_trace_mode", 0)

    # Photo op here
    pymol.cmd.scene(key="PeptideEdm1sig", action="store")
    if ray:
        PeptideEdm1sig = output_dir + "/visualisation/" + "PeptideEdm1sig.png"
        rayTime(PeptideEdm1sig, do_ray=ray)

    pymol.cmd.hide("mesh", "all")

    pymol.cmd.show("mesh", "p_map_05sigma")
    pymol.cmd.color("grey50", "p_map_05sigma")


    pymol.cmd.set_view(viewSet.pMHC_2)
    # Photo op here
    pymol.cmd.scene(key="PeptideEdm05sig", action="store")
    if ray:
        PeptideEdm05sig = output_dir + "/visualisation/" + "PeptideEdm05sig.png"
        rayTime(PeptideEdm05sig, do_ray=ray)

    # Make difference map

    pymol.cmd.hide("mesh", "all")
    pymol.cmd.isomesh("posdiffmesh", file_name + "_dmap", 3.0, "ps", carve=1.6)
    pymol.cmd.color("green", "posdiffmesh")
    pymol.cmd.show("mesh", "posdiffmesh")
    pymol.cmd.isomesh("negdiffmesh", file_name + "_dmap", -3.0, "ps", carve=1.6)
    pymol.cmd.color("red", "negdiffmesh")
    pymol.cmd.show("mesh", "negdiffmesh")

    # Photo op here
    pymol.cmd.set_view(viewSet.pMHC_2)
    pymol.cmd.scene(key="differencemap", action="store")
    if ray:
        differencemap = output_dir + "/visualisation/" + "differencemap.png"
        rayTime(differencemap, do_ray=ray)

    # pMHC helices
    pymol.cmd.hide("mesh", "all")
    pymol.cmd.show("cartoon", "MHCa1")
    pymol.cmd.color(colourSet.generalColourSet["MHCa"], "MHCa1")
    pymol.cmd.set("cartoon_transparency", 0.5)

    # Photo op here
    pymol.cmd.set_view(viewSet.pMHC_2)
    pymol.cmd.scene(key="MHChelixPeptide1", action="store")
    if ray:
        MHChelixPeptide1 = output_dir + "/visualisation/" + "MHChelixPeptid1e.png"
        rayTime(MHChelixPeptide1, do_ray=ray)

    # Copy above scene but with a pink peptide
    pymol.cmd.color("magenta", "ps")
    pymol.cmd.util.cnc("ps")

    # Photo op here
    pymol.cmd.set_view(viewSet.pMHC_2)

    pymol.cmd.scene(key="MHChelixPeptide2", action="store")
    if ray:
        MHChelixPeptide2 = output_dir + "/visualisation/" + "MHChelixPeptide2.png"
        rayTime(MHChelixPeptide2, do_ray=ray)


    # Copy above scene but with a pink peptide
    pymol.cmd.color("magenta", "ps")
    pymol.cmd.util.cnc("ps")

    # Photo op here
    pymol.cmd.set_view(viewSet.pMHC_2)

    pymol.cmd.scene(key="MHChelixPeptide3", action="store")
    if ray:
        MHChelixPeptide3 = output_dir + "/visualisation/" + "MHChelixPeptide3.png"
        rayTime(MHChelixPeptide3, do_ray=ray)

    # Save the session
    pymol.cmd.save(output_dir + "/visualisation/" + file_name + "_peptideMHCvis.pse")

    # Quit pymol
    #pymol.cmd.quit()
    print('     ~  End peptideMHCvisualisation.py v0.1 BETA  ~')

    return None

def apbs_electrostatics(pdb_name, pdb, do_apbs, do_ray, template, mhc_class, neg_colour, pos_colour):

    pymol.cmd.reinitialize()
    initialisePymol()

    print(pdb)
    print(pdb_name)
    
    if os.path.exists(pdb_name) == False:
        os.mkdir(pdb_name)
    
    if os.path.exists(pdb_name+"/apbs") == False:
        os.mkdir(pdb_name+"/apbs")
        
    apbs_dir = pdb_name+"/apbs"
    
    logger = open(pdb_name+"/apbs/apbs_log.txt", "w")

    pymol.cmd.load(pdb, "structure")
    pymol.cmd.load(template, "align_target")
    pMHC_align("structure", "align_target", mhc_class)
    # pymol.cmd.delete('align_target')
    pymol.cmd.remove("hetatm")

    pymol.cmd.select("pMHC", "structure and chain A+B+C")
    pymol.cmd.create("pMHC"+"o", "pMHC") 
    pMHC_create_objects("structure", mhc_class)
    
    pymol.cmd.save(pdb_name+"/apbs/"+pdb_name+"_electrostatics.pse")

    obj_name = "pMHCo"
    pqr_name = pdb_name+"_"+obj_name+".pqr"
    pqr_name_full = apbs_dir+"/"+pdb_name+"_"+obj_name+".pqr"
    apbs_output_name = pdb_name+"_"+obj_name+"_map"
    apbs_output_name_full = apbs_dir+"/"+pdb_name+"_"+obj_name+"_map"
    apbs_input_name = pdb_name+"_"+obj_name+"_apbs.in"
    prepped_pdb = apbs_dir+"/"+pdb_name+"_"+obj_name+".pdb"

    pymol.cmd.save(prepped_pdb, obj_name)
    
    if do_apbs == True:
        
        print("Converting pdb file of object to .pqr using subprocess.call to p2b2pqr")
        command = " ".join(["pdb2pqr30", prepped_pdb, pqr_name_full,"--ff=AMBER"])
        print(command)
        subprocess.call([command], shell=True, stdout=logger, stderr=logger)

        print("######## END pdb2pqr30 ##########")

        with open("bin/data/apbs_template.txt", 'r') as file :
            filedata = file.read()
        
        filedata = filedata.replace('PQR_IN', pqr_name)
        filedata = filedata.replace('APBS_OUT', apbs_output_name)
        
        with open(apbs_dir+"/"+apbs_input_name, 'w') as file:
            file.write(filedata)

        print('apbs seems to have issues where .in file is not the working directory.')
        print('To work around, changing the working directory to the apbs directory using os.chdir - promise to change it straight back after.')

        print(os.getcwd())
        revert = os.getcwd()
        os.chdir(apbs_dir)
        print(os.getcwd())

        print("Calculating electrostatic map using subprocess.call to apbs...")

        command = " ".join(["apbs", apbs_input_name])
        print(command)
        subprocess.call([command], shell=True, stdout=logger, stderr=logger)

        print("######## END apbs ##########")
        os.chdir(revert)
        print(os.getcwd())

        print("Done!")
    

    pymol.cmd.load(apbs_output_name_full+".dx", "map")
    pymol.cmd.ramp_new("map_esramp", "map", [-5,0,5])
    pymol.cmd.set("surface_color", "map_esramp", "structure_pMHCgroove_obj")    
    pymol.cmd.set("surface_ramp_above_mode", 1)

    pymol.cmd.show("surface" , "structure_pMHCgroove_obj")
    pymol.cmd.colour("grey90" , "structure_MHCa1a2_obj")
    pymol.cmd.colour("grey90" , "structure_p_obj")

    pymol.cmd.set("stick_transparency", 0.7, "structure_p_obj")
    pymol.cmd.set("cartoon_transparency", 0.7, "structure_MHCa1a2_obj")

    layer1 = {"structure_pMHCgroove_obj": ["surface"]
             }
    layer2 = {"structure_MHCa1a2_obj": ["cartoon"],
              "structure_p_obj" : ["sticks"]
             }

    layers = [layer1, layer2]

    # MHC birdseye view of electrostatics
    pymol.cmd.disable("*_esramp")
    magickray(layers, "electrostatic_pMHC", viewSet.birdsEyeView, pdb_name+"/apbs", do_ray)

    spike1 = {"structure_p_obj": ["surface"]}
    magick_outline(spike1, viewSet.birdsEyeView, "electrostatic_pMHC", "electrostatic_pMHC3", pdb_name+"/apbs", ray)
    spike2 = {"structure_MHCgroove and chain A": ["surface"]}
    magick_outline(spike2, viewSet.birdsEyeView, "electrostatic_pMHC", "electrostatic_pMHC4", pdb_name+"/apbs", ray)
    spike3 = {"structure_MHCgroove and chain A": ["surface"]}
    magick_outline(spike3, viewSet.birdsEyeView, "electrostatic_pMHC", "electrostatic_pMHC5", pdb_name+"/apbs", ray)

    pymol.cmd.save(pdb_name+"/apbs/"+pdb_name+"_electrostatics.pse")

    return None



###### script specific functions ########
        
def figure1(structure, do_ray, peptide_colour):
    
    # Set colours
    pymol.cmd.color("gray90",structure+" and chain "+chains_dict["MHCA"])
    
    pymol.cmd.color("gray90",structure+"_MHCa1_obj")
    pymol.cmd.color("gray90",structure+"_MHCa2_obj")
    pymol.cmd.color("gray90",structure+"_MHCgroove_obj")
   
    pymol.cmd.color("gray60",structure+" and chain "+chains_dict["MHCB"])
    
    pymol.cmd.color(peptide_colour,structure+" and chain "+chains_dict["PEPTIDE"])
    pymol.cmd.util.cnc(structure+" and chain "+chains_dict["PEPTIDE"])   
    
    #### VIEWS ####
    
    #### IMAGE 1 ####
    # overview of all structure
    
    pymol.cmd.hide("everything", "all")
    
    pymol.cmd.show("cartoon",structure+" and chain "+chains_dict["MHCA"]+" or "+structure+" and chain "+chains_dict["MHCB"])
    pymol.cmd.show("sticks",structure+" and chain "+chains_dict["PEPTIDE"]) 


    if peptide_name in ["X2_A2", "X2_A2_omicron"]:

        pymol.cmd.unbond(structure+" and chain "+chains_dict["PEPTIDE"]+" and resi 3 and name SG and altloc B", structure+" and chain "+chains_dict["PEPTIDE"]+" and resi 3 and name SG and altloc A")
        pymol.cmd.unbond(structure+" and resn DHL", structure+" and chain "+chains_dict["PEPTIDE"]+" and resi 3")
        pymol.cmd.bond(structure+" and resn DHL and name SG and altloc A", structure+" and chain "+chains_dict["PEPTIDE"]+" and resi 3 and name SG and altloc A" )
        pymol.cmd.bond(structure+" and resn DHL and name SG and altloc B", structure+" and chain "+chains_dict["PEPTIDE"]+" and resi 3 and name SG and altloc B" )


        pymol.cmd.set("stick_transparency" , 0.5 , structure+" and resn DHL")

        pymol.cmd.select("SGs", structure+" and chain C and name SG")
        pymol.cmd.set("stick_transparency" , 0.5 , "SGs")
    
    pymol.cmd.set_view(viewSet.pMHC_1)
    pymol.cmd.scene(structure+"_pMHC", "store")
    rayTime(output_dir+"/visualisation/"+ structure+"_pMHC.png", do_ray)
    
    #### IMAGE 2 ####
    # peptide sitting in a1 and a2 helices
    
    pymol.cmd.hide("everything", "all")
    
    pymol.cmd.set("cartoon_transparency", 0.80, structure+"_MHCa1_obj")
    pymol.cmd.set("cartoon_transparency", 0.95, structure+"_MHCa2_obj")
    
    pymol.cmd.show("cartoon",structure+"_MHCa1_obj")
    pymol.cmd.show("cartoon",structure+"_MHCa2_obj")
    pymol.cmd.show("sticks",structure+" and chain "+chains_dict["PEPTIDE"]) 
    
    pymol.cmd.set_view(viewSet.pMHC_2)
    pymol.cmd.scene(structure+"_pMHCa1", "store")
    rayTime(output_dir+"/visualisation/"+ structure+"_pMHCa1.png", do_ray)
    
    #### IMAGE 3 ####
    
    pymol.cmd.set("cartoon_transparency", 0.80, structure+"_MHCgroove_obj")
    pymol.cmd.set("cartoon_transparency", 0.95, structure+"_MHCa2_obj")
    
    pymol.cmd.hide("everything", "all")
    
    pymol.cmd.show("cartoon",structure+"_MHCgroove_obj")
    pymol.cmd.show("cartoon",structure+"_MHCa2_obj")
    pymol.cmd.show("sticks",structure+" and chain "+chains_dict["PEPTIDE"]) 
    
    pymol.cmd.set_view(viewSet.pMHC_2)
    pymol.cmd.scene(structure+"_pMHCgroove", "store")

    rayTime(output_dir+"/visualisation/"+ structure+"_pMHCgroove.png", do_ray)

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


figure1(structure, do_ray, peptide_name)

pymol.cmd.save(output_dir+"/"+pdb_name+".pse")


#######################
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

# pymol.cmd.hide("sticks" , "all")


rayTime(output_dir+"/visualisation/"+ structure+"_cutaway_back.png", do_ray)
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


print("Merging sub-images using subprocess.call to ImageMagick")
command = " ".join(["convert", output_dir+"/visualisation/"+ structure+"_cutaway_back.png", output_dir+"/visualisation/"+ structure+"_cutaway_front.png","-background", "none", "-flatten", output_dir+"/visualisation/"+ structure+"_cutaway.png"])
subprocess.call([command], shell=True)


# pymol.cmd.select("near_copy1_P1_Lys", "byres chain A within 8.0 of (chain C and resi 3) or byres chain B within 8.0 of (chain C and resi 3) ")
# pymol.cmd.show("sticks", "near_copy1_P1_Lys")
# pymol.cmd.util.cnc("near_copy1_P1_Lys")


#################################
# Clip cut-away P1 focus        #
#################################

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

print("Merging sub-images using subprocess.call to ImageMagick")
command = " ".join(["convert", output_dir+"/visualisation/"+ structure+"_P1_focus_cutaway_surface.png", output_dir+"/visualisation/"+ structure+"_P1_focus_cutaway_pocket.png", output_dir+"/visualisation/"+ structure+"_P1_focus_cutaway_peptide.png","-background", "none", "-flatten", output_dir+"/visualisation/"+ structure+"_P1_focus_cutaway.png"])
subprocess.call([command], shell=True)

################
# B-factor     #
################

if os.path.exists(output_dir+"/bfactors") == False:
    os.mkdir(output_dir+"/bfactors")

import bin.bfactors as bfactors
selection = structure+ " and chain "+chains_dict["PEPTIDE"]+" and bb."
print("Using pymol.cmd.iterate to extract backbone atom b factors..")
backbone_b_df = bfactors.extract_b(selection)

print("Min b factor is: " , backbone_b_df["bfactors"].min())
print("Max b factor is: " , backbone_b_df["bfactors"].max())

print("Plotting scatter plot of backbone atom b factors from dataframe..")
bfactors.plot_b_factors(backbone_b_df, peptide_colours[peptide_name], output_dir+"/bfactors/"+pdb_name+"_backbone_bfactors.pdf")


########################
# Electrostatics       #
########################


#### electrostatics.pMHC_electrostatics(pdb, do_apbs)

# apbs_electrostatics(pdb_name=pdb_name, pdb=pdb, do_apbs=do_apbs, do_ray=do_ray, template=template, mhc_class=mhc_class, neg_colour="blue", pos_colour="red")



################
# Maps         #
################

pymol.cmd.reinitialize()


import shutil

pMHC_omit_maps(pdb, mtz, mhc_class="II", chains=chains_dict, ray=do_ray, file_name=pdb_name, peptide_name=peptide_name, output_dir=output_dir)
pymol.cmd.save(output_dir+"/"+pdb_name+"_omit.pse")

###########
# This is the end. We finally want to save the PyMOL session file and quit pymol. 
# If we have asked to view the results in PyMOL, we will finish by opening up the session file.
###########

# pymol.cmd.quit()
# if view == True:
#     subprocess.call(["pymol", output_dir+"/"+pdb_name+".pse"])




################
# Contacts     #
################

peptide_colour = peptide_colours[peptide_name]
print(peptide_colour)
print(type(peptide_colour))
subprocess.call(["python", "bin/projectContactsFullpMHC.py", "--PDB",args.pdb, "--MHCclass", mhc_class, "--chains" , chains, "--peptide_colour" , peptide_colour, "--output_dir" , output_dir])

### everything above here works 28/03/2022
print("~"* 30 , "\n" , "everything above here works 28/03/2022" , "\n" , "~" * 30)


################
# BSA          #
################

import bin.bsa as bsa
import pandas as pd

if do_bsa == True:

    interfaces = [[chains_dict["PEPTIDE"],chains_dict["MHCA"]], [chains_dict["PEPTIDE"], chains_dict["MHCB"]]]
    labels = ["DRA", "DRB"]
    colours = [(0.5,0.5,0.5), (0.7,0.7,0.7)]

    area_dict , asa_bsa_dG_dicts = bsa.calculate_bsa(sesh_name = pdb_name, pdb = pdb, interfaces_of_interest = interfaces, labels = labels, colours = colours, do_pisa = True, output_dir=output_dir)

    # This plots a BSA donut - not really useful for this but it was leftover from another analysis
    area_dict_df = pd.DataFrame(area_dict)
    bsa.plot_donut(area_dict_df, saveas=output_dir+"/bsa/"+"bsa_donut.pdf")

    # This plots a bar chart with BSA going down and ASA up (hopefully)
    peptide_MHCA_df = pd.DataFrame(asa_bsa_dG_dicts[chains_dict["PEPTIDE"]+"_"+chains_dict["MHCA"]]).T
    peptide_MHCB_df = pd.DataFrame(asa_bsa_dG_dicts[chains_dict["PEPTIDE"]+"_"+chains_dict["MHCB"]]).T
    
    print(peptide_MHCA_df)
    print(peptide_MHCB_df)

    # merge the BSA of C_A and C_B. We are not interested in per residue per chain data atm

    peptide_MHCA_bsa = peptide_MHCA_df["BSA"].tolist()
    peptide_MHCB_bsa = peptide_MHCB_df["BSA"].tolist()


    peptide_MHCA_peptide_MHCB_bsa = []
    for x, y in zip(peptide_MHCA_bsa, peptide_MHCB_bsa):
        peptide_MHCA_peptide_MHCB_bsa.append(x+y)

    for line in peptide_MHCA_peptide_MHCB_bsa:
        print(line)


    df = pd.DataFrame.from_dict({
                            'position' : peptide_MHCA_df['ResNum'],
                            'residue'  : peptide_MHCB_df['ResCode'],
                            'BSA' : peptide_MHCA_peptide_MHCB_bsa,
                            'ASA' : peptide_MHCA_df['ASA']
                            })

    print(df)

    bsa.plot_bsa_asa(df=df, colour=peptide_colour, saveas=output_dir+"/bsa/"+"asa_bsa_updown.pdf")


#
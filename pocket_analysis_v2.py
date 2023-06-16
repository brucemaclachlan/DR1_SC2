# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written using Sublime Text & Spyder

@author:    Bruce J MacLachlan, Division of Infection & Immunity, Cardiff University
            ORCID ID: 0000-0002-2685-2733
            Github: https://github.com/brucemaclachlan
            email: maclachlanb@cardiff.ac.uk

            pocket_analysis.py
            
            
"""

description = "\
This script uses PyMOL to make figures which relate to a peptide-HLA-II structure. \
The output is a series of images which describe the HLA-II binding pockets and the peptide that occupies them.\
"

import os
import sys
import time
import subprocess
import pymol
import itertools
import argparse
import colorsys
import shutil


import bin.data.viewSet as viewSet
import bin.data.colourSet as colourSet
import bin.data.peptide_colours as peptide_colours

import bin.colour_functions as colour_functions
import bin.ray_functions as ray_functions
import bin.image_functions as image_functions
import bin.pymol_functions as pymol_functions

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

    parser.set_defaults(do_ray=False)
    parser.set_defaults(do_view=False)
    parser.set_defaults(chains="ABC")
    
    args = parser.parse_args()
    return args

    ################################ def main() ########################################
def main():

    ################################
    # Housekeeping & Startup       #
    ################################

    # Handle command line arguments
    args = parse_args()
    print("Running analysis with the following inputs.. ")
    print(args)
    ray = args.do_ray
    view = args.do_view
    peptide_name = args.peptide
    chains = args.chains

    mhc_class = args.mhc_class

    pdb = args.pdb
    pdb_name = pdb.split(".")[0]
    pdb = os.path.realpath(pdb)

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
    pymol_functions.initialisePymol()

    structure = pdb_name

    pymol_functions.pMHC_load(structure, pdb)

    pymol.cmd.load(template, "align_target")
    pymol_functions.pMHC_align(structure, "align_target", mhc_class, chains_dict)
    # pymol.cmd.delete('align_target')

    chains_dict["MHCA"]
    pymol.cmd.remove(structure+ " and NOT chain "+chains_dict["MHCA"]+","+chains_dict["MHCB"]+","+chains_dict["PEPTIDE"])
    pymol.cmd.remove("hetatm and NOT resn DHL")

    pymol_functions.pMHC_create_objects(structure, mhc_class, chains_dict)

    peptide_colour = colour_functions.hex_to_fraction(peptide_colours.peptide_colours[peptide_name])
    colour_functions.set_new_colour(peptide_name, peptide_colour)

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

    ########################################################################################################################
    # Clip cut-away                                                                                                        #  
    ########################################################################################################################

    pymol.cmd.hide("everything" , "all")

    # must disable depth cue and shadows
    pymol.cmd.set("depth_cue" , 0)

    # this controls the z depth of the slice plane
    # (sets it halfway between the clipping planes)
    fraction = 0.39

    pymol.cmd.set_view("\
        -0.782565236,    0.267936289,   -0.561967373,\
         0.592841208,    0.596321344,   -0.541243196,\
         0.190093905,   -0.756709278,   -0.625500619,\
        -0.000302987,    0.000438260, -193.702819824,\
        64.388168335,   28.303731918,   36.209182739,\
       160.413558960,  259.850311279,   20.000000000" )

    pymol.cmd.show("surface" , structure+"_MHCgroove_obj")

    pymol.cmd.set("ray_interior_color", "grey80")
    pymol.cmd.set("surface_color", "white")

    pymol.cmd.scene(structure+"_cutaway", "store")
    ray_functions.rayTime(output_dir+"/visualisation/"+ structure+"_cutaway_back.png", do_ray)

    # Make a coverer to hide that hole in the clipped surface
    image = Image.new('RGBA', (3000, 3000))
    draw = ImageDraw.Draw(image)
    draw.rectangle((1100, 1300, 1300, 1600), fill = (155,155,155), outline =None)
    image.save(output_dir+"/visualisation/"+ structure+'_coverer.png')

    layer_1 = output_dir+"/visualisation/"+ structure+"_cutaway_back.png"
    layer_2 = output_dir+"/visualisation/"+ structure+"_coverer.png"
    layers = [layer_1, layer_2]
    overlay_output_name = output_dir+"/visualisation/"+ structure+"_cutaway_back_edit.png"

    cutaway_back_edit = image_functions.recursive_alpha_composite(layers, (3000,3000), overlay_output_name)

    # Make a version of the above but with a greenscreen cutaway colour.. this is for the graphical abstract used later  #
    ######################################################################################################################

    greenscreen_color = "0fa110"
    greenscreen_color_srgb = 'sRGB(15,161,16)'
    greenscreen_color_rgb = (15, 161, 16)

    colour_functions.set_new_colour("greenscreen", list((15, 161, 16)))
    pymol.cmd.set("ray_interior_color" , "greenscreen")
    ray_functions.rayTime(output_dir+"/visualisation/"+ structure+"_cutaway_back_green.png", do_ray)

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

    cutaway_back_green_edit = image_functions.recursive_alpha_composite(layers, (3000,3000), overlay_output_name)

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

    ray_functions.rayTime(output_dir+"/visualisation/"+ structure+"_cutaway_front.png", do_ray)

    layer_1 = output_dir+"/visualisation/"+ structure+"_cutaway_back_edit.png"
    layer_2 = output_dir+"/visualisation/"+ structure+"_cutaway_front.png"
    layers = [layer_1, layer_2]
    overlay_output_name = output_dir+"/visualisation/"+ structure+"_cutaway.png"

    cutaway = image_functions.recursive_alpha_composite(layers, (3000,3000), overlay_output_name)

    ###################################################################################################
    # Clip cut-away P1 focus    (for omicron)  i.e. Figure 5                                          #
    ###################################################################################################

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

    pymol.cmd.scene(structure+"_pMHC_P1_old_Lys_focus_surface", "store")
    ray_functions.rayTime(output_dir+"/visualisation/"+ structure+"_P1_old_focus_cutaway_surface.png", do_ray)

    # coverer
    image = Image.new('RGBA', (3000, 3000))
    draw = ImageDraw.Draw(image)
    draw.rectangle((1750, 1400, 2250, 1850), fill = (155,155,155), outline =None)
    image.save(output_dir+"/visualisation/"+ structure+'_old_zoom_coverer.png')

    layer_1 = output_dir+"/visualisation/"+ structure+"_P1_old_focus_cutaway_surface.png"
    layer_2 = output_dir+"/visualisation/"+ structure+"_old_zoom_coverer.png"
    layers = [layer_1, layer_2]
    overlay_output_name = output_dir+"/visualisation/"+ structure+"_P1_old_focus_cutaway_surface_edit.png"

    cutaway_zoom_surface_edit = image_functions.recursive_alpha_composite(layers, (3000,3000), overlay_output_name)

    pymol.cmd.set_view("\
        -0.782565236,    0.267936289,   -0.561967373,\
         0.592841208,    0.596321344,   -0.541243196,\
         0.190093905,   -0.756709278,   -0.625500619,\
        -0.000273535,    0.000260269,  -88.982124329,\
        71.976158142,   20.260526657,   29.843364716,\
       -97.736183167,  155.500396729,   20.000000000" )

    pymol.cmd.hide("surface" , "all")
    pymol.cmd.select("custom_near_P1", structure+" and chain "+chains_dict["MHCA"]+" and resi 7,24,25,26,31,32,43,48 or "+structure+" and chain "+chains_dict["MHCB"]+" and resi 85,89,90,153")

    pymol.cmd.show("sticks", "custom_near_P1")
    pymol.cmd.util.cnc("custom_near_P1")

    pymol.cmd.set("stick_transparency" , 0.5 , structure+" and chain "+chains_dict["MHCA"]+" and resi 31,32,43")

    pymol.cmd.scene(structure+"_pMHC_P1_old_Lys_focus_pocket", "store")
    ray_functions.rayTime(output_dir+"/visualisation/"+ structure+"_P1_old_focus_cutaway_pocket.png", do_ray)

    pymol.cmd.hide("sticks" , "all")

    pymol.cmd.show("sticks" , structure+"_p_obj")

    pymol.cmd.scene(structure+"_pMHC_P1_old_Lys_focus_peptide", "store")
    ray_functions.rayTime(output_dir+"/visualisation/"+ structure+"_P1_old_focus_cutaway_peptide.png", do_ray)

    layer_1 = output_dir+"/visualisation/"+ structure+"_P1_old_focus_cutaway_surface_edit.png"
    layer_2 = output_dir+"/visualisation/"+ structure+"_P1_old_focus_cutaway_pocket.png"
    layer_3 = output_dir+"/visualisation/"+ structure+"_P1_old_focus_cutaway_peptide.png"
    layers = [layer_1, layer_2, layer_3]
    overlay_output_name = output_dir+"/visualisation/"+ structure+"_P1_old_focus_cutaway.png"

    cutaway_zoom = image_functions.recursive_alpha_composite(layers, (3000,3000), overlay_output_name)

    ###################################################################################################
    # Clip cut-away pockets start up                                                                  #
    ###################################################################################################
    
    # These are x,y,z co-ordinates of the Ca atom that generally sits in the groove (based off of 1fyt.pdb template)
    pocket_locs = 		{	"P1"	:	("57.942",  "8.268",   "8.166"),
    						"P2"	:	("54.838",  "10.048",  "6.909"),
    						"P3"	:	("52.624",  "12.228",  "9.059"),
    						"P4"	:	("48.914",  "12.718",  "8.650"),
    						"P5"	:	("47.810",  "16.335",  "8.235"),
    						"P6"	:	("45.527",  "17.809",  "10.826"),
    						"P7"	:	("42.320",  "19.326",  "9.473"),
    						"P8"	:	("41.135",  "22.824",  "10.309"),
    						"P9"	:	("37.875",  "22.873",  "12.235")
    					}

    print("Identifying cores by looking in P1 to P9 pockets..\n")

    pockets = {}

    pymol.cmd.remove("not alt ''+A")
    pymol.cmd.alter("all", "alt=''")

    for pocket in pocket_locs:
        print("Pocket", pocket)
        pockets[pocket] = pymol_functions.find_nearest_residue(structure, pocket_locs[pocket])
        print(pockets[pocket])

    #################################
    # Clip cut-away P1 focus        #
    #################################
    
    current_pocket = "P1"

    pocket_output_dir = output_dir+"/visualisation/P1_pocket/"
    if os.path.exists(pocket_output_dir) == False:
        os.mkdir(pocket_output_dir)

    pymol.cmd.hide("sticks" , "all")
    pymol.cmd.hide("cartoon" , "all")

    pymol.cmd.set("ray_interior_color", "grey80")
    pymol.cmd.set("stick_transparency" , 0.0, structure+" and chain "+chains_dict["MHCA"])
    pymol.cmd.set("stick_transparency" , 0.0, structure+" and chain "+chains_dict["MHCB"])

    pymol.cmd.show("surface" , structure+"_MHCgroove_obj")

    pymol.cmd.set_view("\
        -0.782565236,    0.267936289,   -0.561967373,\
         0.592841208,    0.596321344,   -0.541243196,\
         0.190093905,   -0.756709278,   -0.625500619,\
        -0.000251163,    0.000273665,  -88.982124329,\
        74.280212402,   22.479646683,   32.407936096,\
        55.663780212,  155.100418091,   20.000000000")

    pymol.cmd.scene(structure+"_pMHC_P1_Lys_focus_surface", "store")
    ray_functions.rayTime(pocket_output_dir+ structure+"_P1_focus_cutaway_surface.png", do_ray)

    # coverer
    image = Image.new('RGBA', (3000, 3000))
    draw = ImageDraw.Draw(image)
    draw.rectangle((1750, 1400, 2250, 1850), fill = (155,155,155), outline =None)
    image.save(pocket_output_dir+ structure+'zoom_coverer.png')

    layer_1 = pocket_output_dir+ structure+"_P1_focus_cutaway_surface.png"
    layer_2 = pocket_output_dir+ structure+"zoom_coverer.png"
    layers = [layer_1, layer_2]
    overlay_output_name = pocket_output_dir+ structure+"_P1_focus_cutaway_surface_edit.png"

    cutaway_zoom_surface_edit = image_functions.recursive_alpha_composite(layers, (3000,3000), overlay_output_name)

    pymol.cmd.set_view("\
        -0.782565236,    0.267936289,   -0.561967373,\
         0.592841208,    0.596321344,   -0.541243196,\
         0.190093905,   -0.756709278,   -0.625500619,\
        -0.000273535,    0.000260269,  -88.982124329,\
        71.976158142,   20.260526657,   29.843364716,\
       -97.736183167,  155.500396729,   20.000000000" )

    pymol.cmd.hide("surface" , "all")

    pymol.cmd.select("back_near_P1", structure+" and chain "+chains_dict["MHCA"]+" and resi 24,31,32,43,48 or "+structure+" and chain "+chains_dict["MHCB"]+" and resi 153")
    
    pymol.cmd.show("sticks", "back_near_P1")
    pymol.cmd.util.cnc("back_near_P1")

    pymol.cmd.scene(structure+"_pMHC_P1_Lys_focus_back_pocket", "store")
    ray_functions.rayTime(pocket_output_dir+structure+"_P1_focus_cutaway_back_pocket.png", do_ray)

    pymol.cmd.hide("sticks" , "all")

    pymol.cmd.set("stick_transparency" , 0.0, structure+" and chain "+chains_dict["MHCA"])
    pymol.cmd.set("stick_transparency" , 0.0, structure+" and chain "+chains_dict["MHCB"])

    pymol.cmd.select("front_near_P1", structure+" and chain "+chains_dict["MHCA"]+" and resi 7,25,26 or "+structure+" and chain "+chains_dict["MHCB"]+" and resi 82,85,89,90")

    pymol.cmd.show("sticks", "front_near_P1")
    pymol.cmd.util.cnc("front_near_P1")

    pymol.cmd.set("stick_transparency" , 0.6 , structure+" and chain "+chains_dict["MHCB"]+" and resi 82")

    pymol.cmd.scene(structure+"_pMHC_P1_Lys_focus_front_pocket", "store")
    ray_functions.rayTime(pocket_output_dir+structure+"_P1_focus_cutaway_front_pocket.png", do_ray)

    pymol.cmd.hide("sticks" , "all")

    pymol.cmd.show("cartoon" , structure+" and chain "+chains_dict["PEPTIDE"])
    
    pymol.cmd.set("cartoon_transparency", 0.8, structure+" and chain "+chains_dict["PEPTIDE"])
    
    pymol.cmd.show("sticks" , structure+" and chain "+chains_dict["PEPTIDE"]+" and resi "+str(pockets[current_pocket]["resnum"]))
    
    pymol.cmd.set("stick_transparency", 0.0, structure+" and chain "+chains_dict["PEPTIDE"]+" and resi "+str(pockets[current_pocket]["resnum"]))

    pymol.cmd.scene(structure+"_pMHC_P1_Lys_focus_peptide", "store")
    ray_functions.rayTime(pocket_output_dir+ structure+"_P1_focus_cutaway_peptide.png", do_ray)

    layer_1 = pocket_output_dir+ structure+"_P1_focus_cutaway_surface_edit.png"
    layer_2 = pocket_output_dir+ structure+"_P1_focus_cutaway_back_pocket.png"
    layer_3 = pocket_output_dir+ structure+"_P1_focus_cutaway_peptide.png"
    layer_4 = pocket_output_dir+ structure+"_P1_focus_cutaway_front_pocket.png"

    layers = [layer_1, layer_2, layer_3, layer_4]
    overlay_output_name = pocket_output_dir+ structure+"_P1_focus_cutaway.png"

    cutaway_zoom = image_functions.recursive_alpha_composite(layers, (3000,3000), overlay_output_name)

    #################################
    # Clip cut-away P4 focus        #
    #################################
    current_pocket = "P4"

    pocket_output_dir = output_dir+"/visualisation/P4_pocket/"
    if os.path.exists(pocket_output_dir) == False:
        os.mkdir(pocket_output_dir)

    pymol.cmd.set("ray_interior_color", "grey80")

    pymol.cmd.hide("sticks" , "all")
    pymol.cmd.hide("cartoon" , "all")

    pymol.cmd.set("stick_transparency" , 0.0, structure+" and chain "+chains_dict["MHCA"])
    pymol.cmd.set("stick_transparency" , 0.0, structure+" and chain "+chains_dict["MHCB"])

    pymol.cmd.show("surface" , structure+"_MHCgroove_obj")

    pymol.cmd.set_view("\
        -0.782565236,    0.267936289,   -0.561967373,\
         0.592841208,    0.596321344,   -0.541243196,\
         0.190093905,   -0.756709278,   -0.625500619,\
        -0.000213172,    0.000276900,  -89.006164551,\
        68.204772949,   28.622497559,   32.589332581,\
        55.663780212,  155.100418091,   20.000000000")

    pymol.cmd.scene(structure+"_pMHC_P4_focus_surface", "store")
    ray_functions.rayTime(pocket_output_dir+ structure+"_P4_focus_cutaway_surface.png", do_ray)

    image = Image.new('RGBA', (3000, 3000))
    draw = ImageDraw.Draw(image)
    draw.rectangle(((980, 1625), (1400, 1965)), fill = (155,155,155), outline =None)
    image.save(pocket_output_dir+ structure+'_P4_zoom_coverer.png')

    layer_1 = pocket_output_dir+ structure+"_P4_focus_cutaway_surface.png"
    layer_2 = pocket_output_dir+ structure+"_P4_zoom_coverer.png"
    layers = [layer_1, layer_2]
    overlay_output_name = pocket_output_dir+ structure+"_P4_focus_cutaway_surface_edit.png"

    cutaway_zoom_surface_edit = image_functions.recursive_alpha_composite(layers, (3000,3000), overlay_output_name)

    pymol.cmd.set_view("\
        -0.782565236,    0.267936289,   -0.561967373,\
         0.592841208,    0.596321344,   -0.541243196,\
         0.190093905,   -0.756709278,   -0.625500619,\
        -0.000213172,    0.000276900,  -89.006164551,\
        68.204772949,   28.622497559,   32.589332581,\
        10.423193932,  200.340988159,   20.000000000")

    pymol.cmd.hide("surface" , "all")
    pymol.cmd.select("back_near_P4", structure+" and chain "+chains_dict["MHCA"]+" and resi 9 or "+structure+" and chain "+chains_dict["MHCB"]+" and resi 13,26,78")

    pymol.cmd.show("sticks", "back_near_P4")
    pymol.cmd.util.cnc("back_near_P4")

    pymol.cmd.scene(structure+"_pMHC_P4_Lys_focus_back_pocket", "store")
    ray_functions.rayTime(pocket_output_dir+ structure+"_P4_focus_cutaway_back_pocket.png", do_ray)

    pymol.cmd.hide("sticks" , "all")

    pymol.cmd.select("front_near_P4", structure+" and chain " +chains_dict["MHCB"]+" and resi 74")

    pymol.cmd.show("sticks", "front_near_P4")
    pymol.cmd.util.cnc("front_near_P4")

    pymol.cmd.scene(structure+"_pMHC_P4_Lys_focus_front_pocket", "store")
    ray_functions.rayTime(pocket_output_dir+ structure+"_P4_focus_cutaway_front_pocket.png", do_ray)

    pymol.cmd.hide("sticks" , "all")

    pymol.cmd.show("cartoon" , structure+" and chain "+chains_dict["PEPTIDE"])
    pymol.cmd.set("cartoon_transparency", 0.8, structure+" and chain "+chains_dict["PEPTIDE"])

    pymol.cmd.show("sticks" , structure+" and chain "+chains_dict["PEPTIDE"]+" and resi "+str(pockets[current_pocket]["resnum"]))
    pymol.cmd.set("stick_transparency", 0.0, structure+" and chain "+chains_dict["PEPTIDE"]+" and resi "+str(pockets[current_pocket]["resnum"]))

    pymol.cmd.scene(structure+"_pMHC_P4_Lys_focus_peptide", "store")
    ray_functions.rayTime(pocket_output_dir+ structure+"_P4_focus_cutaway_peptide.png", do_ray)

    layer_1 = pocket_output_dir+ structure+"_P4_focus_cutaway_surface_edit.png"
    layer_2 = pocket_output_dir+ structure+"_P4_focus_cutaway_back_pocket.png"
    layer_3 = pocket_output_dir+ structure+"_P4_focus_cutaway_peptide.png"
    layer_4 = pocket_output_dir+ structure+"_P4_focus_cutaway_front_pocket.png"

    layers = [layer_1, layer_2, layer_3, layer_4]
    overlay_output_name = pocket_output_dir+ structure+"_P4_focus_cutaway.png"

    cutaway_zoom = image_functions.recursive_alpha_composite(layers, (3000,3000), overlay_output_name)

    #################################
    # Clip cut-away P6 focus        #
    #################################
    current_pocket = "P6"

    pocket_output_dir = output_dir+"/visualisation/P6_pocket/"
    if os.path.exists(pocket_output_dir) == False:
        os.mkdir(pocket_output_dir)

    pymol.cmd.set("ray_interior_color", "grey80")

    pymol.cmd.hide("sticks" , "all")
    pymol.cmd.hide("cartoon" , "all")
    pymol.cmd.set("stick_transparency" , 0.0, structure+" and chain "+chains_dict["MHCA"])
    pymol.cmd.set("stick_transparency" , 0.0, structure+" and chain "+chains_dict["MHCB"])

    pymol.cmd.show("surface" , structure+"_MHCgroove_obj")

    pymol.cmd.set_view("\
        -0.782565236,    0.267936289,   -0.561967373,\
         0.592841208,    0.596321344,   -0.541243196,\
         0.190093905,   -0.756709278,   -0.625500619,\
        -0.000175398,    0.000281050,  -89.025260925,\
        62.106658936,   33.530090332,   33.852058411,\
        55.663780212,  155.100418091,   20.000000000")
    ### cut above here and paste into script ###)

    pymol.cmd.scene(structure+"_pMHC_P6_focus_surface", "store")
    ray_functions.rayTime(pocket_output_dir+ structure+"_P6_focus_cutaway_surface.png", do_ray)

    image = Image.new('RGBA', (3000, 3000))
    draw = ImageDraw.Draw(image)
    draw.rectangle(((220, 1650), (630, 2030)), fill = (155,155,155), outline =None)
    image.save(pocket_output_dir+ structure+'_P6_zoom_coverer.png')

    layer_1 = pocket_output_dir+ structure+"_P6_focus_cutaway_surface.png"
    layer_2 = pocket_output_dir+ structure+"_P6_zoom_coverer.png"
    layers = [layer_1, layer_2]
    overlay_output_name = pocket_output_dir+ structure+"_P6_focus_cutaway_surface_edit.png"

    cutaway_zoom_surface_edit = image_functions.recursive_alpha_composite(layers, (3000,3000), overlay_output_name)

    pymol.cmd.set_view("\
        -0.782565236,    0.267936289,   -0.561967373,\
         0.592841208,    0.596321344,   -0.541243196,\
         0.190093905,   -0.756709278,   -0.625500619,\
        -0.000175398,    0.000281050,  -89.025260925,\
        62.106658936,   33.530090332,   33.852058411,\
        41.172477722,  169.591690063,   20.000000000")

    pymol.cmd.hide("surface" , "all")
    pymol.cmd.select("back_near_P6", structure+" and chain "+chains_dict["MHCA"]+" and resi 11,62,65,66")

    pymol.cmd.show("sticks", "back_near_P6")
    pymol.cmd.util.cnc("back_near_P6")

    pymol.cmd.scene(structure+"_pMHC_P6_Lys_focus_pocket", "store")
    ray_functions.rayTime(pocket_output_dir+ structure+"_P6_focus_cutaway_back_pocket.png", do_ray)

    pymol.cmd.hide("sticks" , "all")

    pymol.cmd.select("front_near_P6", structure+" and chain "+chains_dict["MHCB"]+" and resi 11")

    pymol.cmd.show("sticks", "front_near_P6")
    pymol.cmd.util.cnc("front_near_P6")

    pymol.cmd.scene(structure+"_pMHC_P6_Lys_front_pocket", "store")
    ray_functions.rayTime(pocket_output_dir+ structure+"_P6_focus_cutaway_front_pocket.png", do_ray)

    pymol.cmd.hide("sticks" , "all")

    pymol.cmd.show("cartoon" , structure+" and chain "+chains_dict["PEPTIDE"])
    pymol.cmd.set("cartoon_transparency", 0.8, structure+" and chain "+chains_dict["PEPTIDE"])
    pymol.cmd.show("sticks" , structure+" and chain "+chains_dict["PEPTIDE"]+" and resi "+str(pockets[current_pocket]["resnum"]))
    pymol.cmd.set("stick_transparency", 0.0, structure+" and chain "+chains_dict["PEPTIDE"]+" and resi "+str(pockets[current_pocket]["resnum"]))

    pymol.cmd.scene(structure+"_pMHC_P6_Lys_focus_peptide", "store")
    ray_functions.rayTime(pocket_output_dir+ structure+"_P6_focus_cutaway_peptide.png", do_ray)

    pymol.cmd.hide("sticks" , "all")

    layer_1 = pocket_output_dir+ structure+"_P6_focus_cutaway_surface_edit.png"
    layer_2 = pocket_output_dir+ structure+"_P6_focus_cutaway_back_pocket.png"
    layer_3 = pocket_output_dir+ structure+"_P6_focus_cutaway_peptide.png"
    layer_4 = pocket_output_dir+ structure+"_P6_focus_cutaway_front_pocket.png"

    layers = [layer_1, layer_2, layer_3, layer_4]
    overlay_output_name = pocket_output_dir+ structure+"_P6_focus_cutaway.png"

    cutaway_zoom = image_functions.recursive_alpha_composite(layers, (3000,3000), overlay_output_name)

    #################################
    # Clip cut-away P9 focus        #
    #################################
    current_pocket = "P9"

    pocket_output_dir = output_dir+"/visualisation/P9_pocket/"
    if os.path.exists(pocket_output_dir) == False:
        os.mkdir(pocket_output_dir)

    pymol.cmd.set("ray_interior_color", "grey80")

    pymol.cmd.hide("sticks" , "all")
    pymol.cmd.hide("cartoon" , "all")

    pymol.cmd.set("stick_transparency" , 0.0, structure+" and chain "+chains_dict["MHCA"])
    pymol.cmd.set("stick_transparency" , 0.0, structure+" and chain "+chains_dict["MHCB"])

    pymol.cmd.show("surface" , structure+"_MHCgroove_obj")

    pymol.cmd.set_view("\
        -0.782565236,    0.267936289,   -0.561967373,\
         0.592841208,    0.596321344,   -0.541243196,\
         0.190093905,   -0.756709278,   -0.625500619,\
        -0.000139549,    0.000291333,  -89.021316528,\
        55.100147247,   36.606567383,   37.478519440,\
        55.663780212,  155.100418091,   20.000000000")

    pymol.cmd.scene(structure+"_pMHC_P9_focus_surface", "store")

    ray_functions.rayTime(pocket_output_dir+ structure+"_P9_focus_cutaway_surface.png", do_ray)

    pymol.cmd.set_view("\
        -0.782565236,    0.267936289,   -0.561967373,\
         0.592841208,    0.596321344,   -0.541243196,\
         0.190093905,   -0.756709278,   -0.625500619,\
        -0.000139549,    0.000291333,  -89.021316528,\
        55.100147247,   36.606567383,   37.478519440,\
        17.303056717,  193.461135864,   20.000000000" )
    ### cut above here and paste into script ###)

    pymol.cmd.hide("surface" , "all")
    pymol.cmd.select("back_near_P9", structure+" and chain "+chains_dict["MHCA"]+" and resi 69,72,73 or "+structure+" and chain "+chains_dict["MHCB"]+" and resi 9")

    pymol.cmd.show("sticks", "back_near_P9")
    pymol.cmd.util.cnc("back_near_P9")

    pymol.cmd.scene(structure+"_pMHC_P9_Lys_focus_back_pocket", "store")
    ray_functions.rayTime(pocket_output_dir+ structure+"_P9_focus_cutaway_back_pocket.png", do_ray)

    pymol.cmd.hide("sticks" , "all")

    pymol.cmd.select("back_near_P9", structure+" and chain "+chains_dict["MHCA"]+" and resi 76 or "+structure+" and chain "+chains_dict["MHCB"]+" and resi 30,57,60,61")
    pymol.cmd.show("sticks", "back_near_P9")
    pymol.cmd.util.cnc("back_near_P9")

    pymol.cmd.set("stick_transparency" , 0.6 , structure+" and chain "+chains_dict["MHCB"]+" and resi 60,61")

    pymol.cmd.scene(structure+"_pMHC_P9_Lys_focus_front_pocket", "store")
    ray_functions.rayTime(pocket_output_dir+ structure+"_P9_focus_cutaway_front_pocket.png", do_ray)

    pymol.cmd.hide("sticks" , "all")

    pymol.cmd.show("cartoon" , structure+" and chain "+chains_dict["PEPTIDE"])
    pymol.cmd.set("cartoon_transparency", 0.8, structure+" and chain "+chains_dict["PEPTIDE"])
    pymol.cmd.show("sticks" , structure+" and chain "+chains_dict["PEPTIDE"]+" and resi "+str(pockets[current_pocket]["resnum"]))
    pymol.cmd.set("stick_transparency", 0.0, structure+" and chain "+chains_dict["PEPTIDE"]+" and resi "+str(pockets[current_pocket]["resnum"]))

    pymol.cmd.scene(structure+"_pMHC_P9_Lys_focus_peptide", "store")
    ray_functions.rayTime(pocket_output_dir+ structure+"_P9_focus_cutaway_peptide.png", do_ray)

    layer_1 = pocket_output_dir+ structure+"_P9_focus_cutaway_surface.png"
    layer_2 = pocket_output_dir+ structure+"_P9_focus_cutaway_back_pocket.png"
    layer_3 = pocket_output_dir+ structure+"_P9_focus_cutaway_peptide.png"
    layer_4 = pocket_output_dir+ structure+"_P9_focus_cutaway_front_pocket.png"

    layers = [layer_1, layer_2, layer_3, layer_4]
    overlay_output_name = pocket_output_dir+ structure+"_P9_focus_cutaway.png"

    cutaway_zoom = image_functions.recursive_alpha_composite(layers, (3000,3000), overlay_output_name)

    ########### Graphical abstract peptide-HLA (named this PILtoon) ###########

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
    pymol.cmd.show("cartoon" ,  structure+"_p_obj")
    pymol.cmd.show("sticks" ,  structure+"_p_obj and sc."+" or "+ structure+"_p_obj and name ca")
    pymol.cmd.hide("everything", "hetatm")

    pymol.cmd.set("cartoon_loop_radius",0.4)
    pymol.cmd.set("stick_radius" , 0.4)

    pymol.cmd.set_view("\
        -0.782565236,    0.267936289,   -0.561967373,\
         0.592841208,    0.596321344,   -0.541243196,\
         0.190093905,   -0.756709278,   -0.625500619,\
        -0.000302098,    0.000437962, -193.702819824,\
        64.444366455,   28.357858658,   36.271736145,\
      -5579.197753906, 5999.861328125,   20.000000000" )

    colour_functions.set_new_colour("cell_3", list(cell_press_colours["cell_3"]["rgb"]))
    pymol.cmd.color("cell_3", structure+"_p_obj")
    ray_functions.rayTime(output_dir+"/visualisation/"+ structure+"_peptide_spheres.png", do_ray)

    groove_front = image_functions.image_single_colour_filter(output_dir+"/visualisation/"+ structure+"_cutaway_back_green_edit.png" , greenscreen_color_rgb, cell_press_colours["cell_5"]["rgb"]).save(output_dir+"/visualisation/"+ structure+"_PILtoon1.png")
    groove_all = image_functions.image_not_transparent_filter(output_dir+"/visualisation/"+ structure+"_cutaway_back_green_edit.png" , cell_press_colours["cell_2"]["rgb"]).save(output_dir+"/visualisation/"+ structure+"_PILtoon2.png")
    peptide_blob = image_functions.image_not_transparent_filter(output_dir+"/visualisation/"+ structure+"_peptide_spheres.png", cell_press_colours["cell_3"]["rgb"]).save(output_dir+"/visualisation/"+ structure+"_PILtoon3.png")

    layer_1 = output_dir+"/visualisation/"+ structure+"_PILtoon2.png"
    layer_2 = output_dir+"/visualisation/"+ structure+"_PILtoon1.png"
    layer_3 = output_dir+"/visualisation/"+ structure+"_PILtoon3.png"

    overlay_output_name = output_dir+"/visualisation/"+ structure+"_PILtoon.png"

    layers = [layer_1, layer_2, layer_3]

    piltoon = image_functions.recursive_alpha_composite(layers, (3000,3000), overlay_output_name)

    if peptide_name in ["X2_A2_omicron", "SC2_6_omicron"]:
        pymol.cmd.hide("everything", "all")
        if peptide_name == "X2_A2_omicron":

            pymol.cmd.show("sticks" ,  structure+"_p_obj and resi 8,11,13 and sc."+" or "+ structure+"_p_obj and resi 8,11,13 and name ca")
            pymol.cmd.color(peptide_name, structure+"_p_obj")

        if peptide_name == "SC2_6_omicron":

            pymol.cmd.show("sticks" ,  structure+"_p_obj and sc."+" or "+ structure+"_p_obj and name ca")
            pymol.cmd.show("cartoon" , structure+"_p_obj")
            pymol.cmd.color(peptide_name, structure+"_p_obj")

        ray_functions.rayTime(output_dir+"/visualisation/"+ structure+"_peptide_spheres_omicron.png", do_ray)
        peptide_residue_blob = image_functions.image_not_transparent_filter(output_dir+"/visualisation/"+ structure+"_peptide_spheres_omicron.png", (198,99,143)).save(output_dir+"/visualisation/"+ structure+"_PILtoon4.png")

        layer_4 = output_dir+"/visualisation/"+ structure+"_PILtoon4.png"
        layers = [layer_1, layer_2, layer_3, layer_4]
        overlay_output_name = output_dir+"/visualisation/"+ structure+"_PILtoon_omicron.png"

        piltoon_omicron = image_functions.recursive_alpha_composite(layers, (3000,3000), overlay_output_name)

    ###########
    # This is the end. We finally want to save the PyMOL session file and quit pymol. 
    # If we have asked to view the results in PyMOL, we will finish by opening up the session file.
    ###########

    pymol.cmd.save(output_dir+"/"+pdb_name+"_pocket_analysis.pse")

    pymol.cmd.quit()
    if view == True:
        subprocess.call(["pymol", output_dir+"/"+pdb_name+"_pocket_analysis.pse"])


if __name__ == '__main__':
    main()
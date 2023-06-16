#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written using Sublime Text & Spyder

@author:    Bruce J MacLachlan, Division of Infection & Immunity, Cardiff University
            ORCID ID: 0000-0002-2685-2733
            Github: https://github.com/brucemaclachlan
            email: maclachlanb@cardiff.ac.uk
            
            Analysis scripts for HLA-DR1 SARS-CoV-2 structural analyses
"""

description = "\
This script uses PyMOL to make figures relating to a peptide-HLA structure.\
These figures include general overall views on the peptide-HLA complex and \
visualisations of electron density maps focusing on the peptide atoms.\
An omit map (peptide atoms removed) analysis is also performed and visualised.\
"

import os
import sys
import time
import subprocess
import pymol
import itertools
import argparse
import colorsys
import glob
import shutil

import bin.data.viewSet as viewSet
import bin.data.colourSet as colourSet
import bin.data.peptide_colours as peptide_colours

import bin.colour_functions as colour_functions
import bin.ray_functions as ray_functions
import bin.image_functions as image_functions
import bin.pymol_functions as pymol_functions

def parse_args():
    """parse arguments for command line"""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--ray', dest = 'do_ray',   action='store_true', required = False, help ='Do you want to render and save images using ray? Add flag to generate images (time consuming).')
    parser.add_argument('--no_ray', dest = 'do_ray',   action='store_false', required = False, help ='Do you want to render and save images? Add flag to generate no images (quicker).')
    parser.add_argument('--view', dest = 'do_view',   action='store_true', required = False, help ='Do you want to view the pymol session on end? Add flag to open the pymol session at end.')
    parser.add_argument('--mhc_class', dest = 'mhc_class',   type=str, required = True, help ='What is the MHC class of the structure (options: I or II')
    parser.add_argument('--peptide', dest = 'peptide',   type=str, required = True, help ='What is the name of the peptide. This is required to look up the peptide colour. Please ensure the peptide name is listed in peptide_colours.py')
    parser.add_argument('--chains', dest = 'chains',   type=str, required = False, help ='What are the chains of the MHC complex, i.e. A = DR1a, B = DR1b, C = peptide. ABC as here is default.')

    parser.add_argument('--omit_mode', type=str, required = False, help ='How do you want to run the omit map refinement (options: refmac or phenix) fast or slow respectively. Phenix more robust (simmulated annealing)')
    parser.add_argument('--pdb', dest = 'pdb',   type=str, required = True, help ='Input .pdb file')
    parser.add_argument('--mtz', dest = 'mtz',   type=str, required = True, help ='Input .mtz file')

    parser.set_defaults(do_ray=False)
    parser.set_defaults(do_view=False)
    parser.set_defaults(omit_mode="refmac")
    parser.set_defaults(chains="ABC")
    
    args = parser.parse_args()
    return args

def figure1(structure, do_ray, peptide_name, chains_dict, output_dir):
    
    # Set colours
    pymol.cmd.color("gray90",structure+" and chain "+chains_dict["MHCA"])
    pymol.cmd.color("gray90",structure+"_MHCa1_obj")
    pymol.cmd.color("gray90",structure+"_MHCa2_obj")
    pymol.cmd.color("gray90",structure+"_MHCgroove_obj")
    pymol.cmd.color("gray60",structure+" and chain "+chains_dict["MHCB"])
    
    pymol.cmd.color(peptide_name,structure+" and chain "+chains_dict["PEPTIDE"])
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
    ray_functions.rayTime(output_dir+"/visualisation/"+ structure+"_pMHC.png", do_ray)
    
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
    ray_functions.rayTime(output_dir+"/visualisation/"+ structure+"_pMHCa1.png", do_ray)
    
    #### IMAGE 3 ####
    
    pymol.cmd.set("cartoon_transparency", 0.80, structure+"_MHCgroove_obj")
    pymol.cmd.set("cartoon_transparency", 0.95, structure+"_MHCa2_obj")
    
    pymol.cmd.hide("everything", "all")
    
    pymol.cmd.show("cartoon",structure+"_MHCgroove_obj")
    pymol.cmd.show("cartoon",structure+"_MHCa2_obj")
    pymol.cmd.show("sticks",structure+" and chain "+chains_dict["PEPTIDE"]) 
    
    pymol.cmd.set_view(viewSet.pMHC_2)
    pymol.cmd.scene(structure+"_pMHCgroove", "store")

    ray_functions.rayTime(output_dir+"/visualisation/"+ structure+"_pMHCgroove.png", do_ray)

    return None

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
    omit_mode = args.omit_mode

    pdb = args.pdb
    pdb_name = pdb.split(".")[0]
    pdb = os.path.realpath(pdb)

    mtz = args.mtz
    mtz_name = mtz.split(".")[0]
    mtz = os.path.realpath(mtz)

    bin_dir = os.getcwd()+"/bin"
    home_dir = os.getcwd()

    output_dir = os.getcwd()+"/"+pdb_name+"_"+chains

    # Create ouput directories if not present
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

    pymol.cmd.remove(structure+ " and NOT chain "+chains_dict["MHCA"]+","+chains_dict["MHCB"]+","+chains_dict["PEPTIDE"])
    pymol.cmd.remove("hetatm and NOT resn DHL")

    pymol_functions.pMHC_create_objects(structure, mhc_class, chains_dict)

    peptide_colour = colour_functions.hex_to_fraction(peptide_colours.peptide_colours[peptide_name])
    colour_functions.set_new_colour(peptide_name, peptide_colour)

    figure1(structure, do_ray, peptide_name, chains_dict, output_dir)

    pymol.cmd.save(output_dir+"/"+pdb_name+".pse")


    ################
    # Maps         #
    ################

    pymol.cmd.reinitialize()

    import bin.omit_maps as omit
    if omit_mode == "refmac":
        edmap, diffmap = omit.pMHC_omit_maps_refmac(pdb=pdb, mtz=mtz, mhc_class=mhc_class, chains_dict=chains_dict, file_name=pdb_name, output_dir=output_dir)

    elif omit_mode == "phenix":
        edmap, diffmap = omit.pMHC_omit_maps_phenix(pdb=pdb, mtz=mtz, mhc_class=mhc_class, chains_dict=chains_dict, file_name=pdb_name, output_dir=output_dir)

    pymol.cmd.reinitialize()
    omit.pMHC_omit_maps_visualise(pdb=pdb, edmap=edmap, diffmap=diffmap, mhc_class=mhc_class, ray=do_ray, chains_dict=chains_dict, output_dir=output_dir, file_name=pdb_name, peptide_name=peptide_name,)
    pymol.cmd.save(output_dir+"/"+pdb_name+"_omit.pse")
    
    ##########
    # This is the end. We finally want to save the PyMOL session file and quit pymol. 
    # If we have asked to view the results in PyMOL, we will finish by opening up the session file.
    ###########

    pymol.cmd.quit()
    if view == True:
        subprocess.call(["pymol", output_dir+"/"+pdb_name+".pse"])

if __name__ == '__main__':
    main()
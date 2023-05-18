#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 11:09:10 2019
Python 2.7 Anaconda Distribution recommended
Written using Spyder

@author:    Bruce J MacLachlan, Monash Biomedicine Discovery Institute, Clayton, VIC, 3800
            ORCID ID: 0000-0002-2685-2733
            Github: https://github.com/brucemaclachlan
            email: bruce.maclachlan@monash.edu
"""
description = "\
"

import pymol
import os
import subprocess
import sys
import itertools
import time
import argparse
import shutil
import urllib
import bin.colourSet as colourSet
import bin.viewSet as viewSet


def parse_args():
    """parse arguments for command line"""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--pdb', dest = "pdb", type=str, required=True)
    parser.add_argument('--mtz', dest = "mtz", type=str, required=True)
    parser.add_argument('--MHCclass', dest = "MHCclass", type=str, required=True)
    parser.add_argument('--chains', dest = "chains", type=str, required=False, default="ABCDE")
    parser.add_argument('--d', dest = "dir", type=str, required=True)
    parser.add_argument('--ray', dest = 'do_ray',   action='store_true', required = False, help ='Do you want to render and save images? Add flag to generate images (time consuming).')
    parser.add_argument('--no_ray', dest = 'do_ray',   action='store_false', required = False, help ='Do you want to render and save images? Add flag to generate no images (quicker).')
    parser.add_argument('--view', dest = 'do_view',   action='store_true', required = False, help ='Do you want to view the pymol session on end? Add flag to leave open at end.')
    parser.add_argument('--peptide_colour', dest = 'peptide_colour',type=str, required=False, default="", help = "Add a colour for the peptide, if blank, will default to magenta. Supplied colour must be a pymol keyword colour i.e. palegreen")

    parser.set_defaults(do_ray=False)
    parser.set_defaults(do_view=False)
    
    args = parser.parse_args()
    return args



def initialisePymol():
    '''
    This function asks python to start a new pymol session and apply a set of parameters related pymol renders the molecules.
    i.e. I don't like shadows, so they are turned off.
    This helps to keep all figures consistent.
    '''
    print "\nInitialising pymol...\n"
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

def wait4ray(query):  
    counter = 0
    spinner = itertools.cycle(['-', '/', '|', '\\'])
    while not os.path.exists(query):
        toWrite=spinner.next() + " Time elapsed: "+str(counter)+" seconds"
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
        print "Outputting image.. This may take a few seconds.."
        if os.path.exists(saveas):
            print "Removing "+saveas+" as it already exists!"
            os.remove(saveas)
        time.sleep(5)
        pymol.cmd.png(saveas,ray=do_ray,width=3000,height=3000, dpi=200)
        wait4ray(saveas) 
        print "Done! "+str(saveas)+ " was outputted" 
        return None
    
###### BODY ########
        
args = parse_args()
pdb = os.path.realpath(args.pdb)
mtz = os.path.realpath(args.mtz)
chains = args.chains
file_name = args.dir
MHCclass = args.MHCclass
ray = args.do_ray
view = args.do_view
peptide_colour=args.peptide_colour

# need to add this line as pymol.png wants 0 or 1 not bool
if ray == True:
    do_ray = 1
if ray == False:
    do_ray = 0


pdb_name = pdb.rsplit('.', 1)[0].lower()

# Sort chains
MHCachain, MHCbchain, peptidechain, TCRachain, TCRbchain = chains[0], chains[1], chains[2], chains[3], chains[4]

# Make output folder #

if not os.path.exists(file_name):
    print "Creating Directory " + file_name
    os.makedirs(file_name)

if not os.path.exists(file_name + "/visualisation"):
    print "Creating Directory " + file_name + "/visualisation"
    os.makedirs(file_name + "/visualisation")

if not os.path.exists(file_name + "/maps"):
    print "Creating Directory " + file_name + "/maps"
    os.makedirs(file_name + "/maps")

if not os.path.exists(file_name + "/pdbs"):
    print "Creating Directory " + file_name + "/pdbs"
    os.makedirs(file_name + "/pdbs")

if mtz != "ebi":
    print "A map.mtz file was provided!", mtz, "will be moved to", file_name + "/electrostatics/" + file_name + ".mtz"

    shutil.copy(mtz, file_name + "/electrostatics/" + file_name + ".mtz")

if mtz == "ebi":
    if not os.path.exists(file_name + "/electrostatics/" + file_name + ".mtz"):
        print "Downloading map.mtz for entry", pdb_name, "from the PDBe (EBI)"

        try:
            urllib.urlretrieve("http://www.ebi.ac.uk/pdbe/coordinates/files/" + pdb_name + "_map.mtz",
                           file_name + "/electrostatics/" + file_name + ".mtz")
        except:
            print "Could not retrieve url. Please try again, making sure you either supply a file," \
                  " or your file shares its name with one on PDB"
            # quit early, get rid of pymol
            pymol.cmd.quit()
            sys.exit()


tempTxt = "rmchain " + peptidechain + "\n" + "END"
temp = open(file_name + "/pdbs/" + file_name + "_pdbcurPARAM.tmp", "w")
temp.write(tempTxt)
temp.close()

# delete chain
os.system("pdbcur XYZIN " + pdb + " XYZOUT " + file_name + "/pdbs/" + file_name +
          "_nopeptide.pdb" + " < " + file_name + "/pdbs/" + file_name + "_pdbcurPARAM.tmp")
#os.remove(file_name + "/pdbs/" + file_name + "_pdbcurPARAM.tmp")

mtz = file_name + "/electrostatics/" + file_name + ".mtz"


#Run one cycle of refmac without peptide in the groove
os.system("refmac5 XYZIN " + file_name + "/pdbs/" + file_name + "_nopeptide.pdb" + " XYZOUT " +
          file_name + "/pdbs/" + file_name + "_refmac5_omitmap.pdb" + " HKLIN " + file_name + "/electrostatics/"
          + file_name + ".mtz" + " HKLOUT " + file_name + "/electrostatics/" + file_name + "_refmac5_omitmap.mtz" + " LIBOUT "
          + file_name + "/electrostatics/" + file_name + "_refmac_omitmap.cif" + " < " + "bin/refmacOMITparams.tmp")

os.system("mtzdump " + " HKLIN " + file_name + "/electrostatics/" + file_name + "_refmac5_omitmap.mtz END")
mtz = file_name+"/electrostatics/" + file_name + "_refmac5_omitmap.mtz"
# Use CCP4 to generate map

os.system("fft HKLIN " + mtz + " MAPOUT " + file_name + "/electrostatics/" +
          file_name + ".map1.tmp" + " < " + "bin/EDMparam1.tmp")
os.system("mapmask MAPIN " + file_name + "/electrostatics/" + file_name + ".map1.tmp" +
          " MAPOUT " + file_name + "/electrostatics/" + file_name + ".map.ccp4" + " XYZIN " + pdb + " < " + "bin/EDMparam2.tmp")
os.system("fft HKLIN " + mtz + " MAPOUT " + file_name + "/electrostatics/" +
          file_name + ".map3.tmp" + " < " + " bin/EDMparam3.tmp")
os.system("mapmask MAPIN " + file_name + "/electrostatics/" + file_name + ".map3.tmp" +
          " MAPOUT " + file_name + "/electrostatics/" + file_name + ".difference_map.ccp4" +
          " XYZIN " + pdb + " < " + "bin/EDMparam4.tmp")

os.remove(file_name + "/electrostatics/" + file_name + ".map1.tmp")
os.remove(file_name + "/electrostatics/" + file_name + ".map3.tmp")

edmap = file_name + "/electrostatics/" + file_name + ".map.ccp4"

diffmap = file_name + "/electrostatics/" + file_name + ".difference_map.ccp4"

# Sort chains
MHCachain, MHCbchain, peptidechain, TCRachain, TCRbchain = chains[0], chains[1], chains[2], chains[3], chains[4]

# Find the MHC helices
if MHCclass == "I":
    a1locs = range(50, 86)
    MHCa1 = ["MHCa"] + a1locs
    a2locs = range(140, 176)
    MHCa2 = ["MHCa"] + a2locs

if MHCclass == "II":
    a1locs = range(46, 78)
    MHCa1 = ["MHCa"] + a1locs
    a2locs = range(54, 91)
    MHCa2 = ["MHCb"] + a2locs

# Let's get started

initialisePymol()

pymol.cmd.load(pdb, "complex")
pymol.cmd.load(file_name + "/pdbs/" + file_name + "_refmac5_omitmap.pdb", "omitxyz")
pymol.cmd.load(edmap, file_name + "_map")
pymol.cmd.load(diffmap, file_name + "_dmap")

# align to template
print "\nAligning file to template...\n"
pymol.cmd.load("bin/" + MHCclass + "_cdr_template.pdb")
pymol.cmd.align("complex", MHCclass + "_cdr_template")
pymol.cmd.align("omitxyz", MHCclass + "_cdr_template")
pymol.cmd.matrix_copy("omitxyz", "complex")
pymol.cmd.matrix_copy("omitxyz", file_name + "_map")
pymol.cmd.matrix_copy("omitxyz", file_name + "_dmap")

pymol.cmd.delete(MHCclass + "_cdr_template")
print "\nAlignment to " + MHCclass + "_cdr_template.pdb  complete!\n"

# Make chains objects
pymol.cmd.select("MHCas", selection="chain " + MHCachain)
pymol.cmd.select("MHCbs", selection="chain " + MHCbchain)
pymol.cmd.select("ps", selection="chain " + peptidechain)
pymol.cmd.select("TCRas", selection="chain " + TCRachain)
pymol.cmd.select("TCRbs", selection="chain " + TCRbchain)

pymol.cmd.hide("all")
pymol.cmd.show("sticks", "ps")

if peptide_colour == "":
    pymol.cmd.color(colourSet.generalColourSet["p"], "ps")
else:
    pymol.cmd.color(peptide_colour, "ps")
pymol.cmd.util.cnc("ps")

# Select the MHCa helices
#MHCa1, MHCa2 = None, None

locs = '+'.join(str(x) for x in MHCa1[1:])
pymol.cmd.select("MHCa1", selection=MHCa1[0] + "s and resi " + locs)
locs = '+'.join(str(x) for x in MHCa2[1:])
pymol.cmd.select("MHCa2", selection=MHCa2[0] + "s and resi " + locs)

# Make electron density map
pymol.cmd.map_double(file_name + "_map", -1)
pymol.cmd.isomesh("p_map_1sigma", file_name + "_map", 1.0, "ps", carve=1.6)
pymol.cmd.isomesh("p_map_05sigma", file_name + "_map", 0.5, "ps", carve=1.6)
pymol.cmd.set("mesh_width", 0.5)

pymol.cmd.hide("mesh", "all")

pymol.cmd.set_view(viewSet.peptideView)

pymol.cmd.show("mesh", "p_map_1sigma")
pymol.cmd.color("grey50", "p_map_1sigma")

# Photo op here
pymol.cmd.scene(key="PeptideEdm1sig", action="store")
PeptideEdm1sig = file_name + "/visualisation/" + "PeptideEdm1sig.png"
rayTime(PeptideEdm1sig, do_ray)

pymol.cmd.hide("mesh", "all")

pymol.cmd.show("mesh", "p_map_05sigma")
pymol.cmd.color("grey50", "p_map_05sigma")

# Photo op here
pymol.cmd.scene(key="PeptideEdm05sig", action="store")
PeptideEdm05sig = file_name + "/visualisation/" + "PeptideEdm05sig.png"
rayTime(PeptideEdm05sig, do_ray)

# Make difference map

pymol.cmd.hide("mesh", "all")
pymol.cmd.isomesh("posdiffmesh", file_name + "_dmap", 3.0, "ps", carve=1.6)
pymol.cmd.color("green", "posdiffmesh")
pymol.cmd.show("mesh", "posdiffmesh")
pymol.cmd.isomesh("negdiffmesh", file_name + "_dmap", -3.0, "ps", carve=1.6)
pymol.cmd.color("red", "negdiffmesh")
pymol.cmd.show("mesh", "negdiffmesh")

# Photo op here
pymol.cmd.scene(key="differencemap", action="store")
differencemap = file_name + "/visualisation/" + "differencemap.png"
rayTime(differencemap, do_ray)

# pMHC helices
pymol.cmd.hide("mesh", "all")
pymol.cmd.show("cartoon", "MHCa1")
pymol.cmd.color(colourSet.generalColourSet["MHCa"], "MHCa1")
pymol.cmd.set("cartoon_transparency", 0.5)

# Photo op here
pymol.cmd.scene(key="MHChelixPeptide1", action="store")
MHChelixPeptide1 = file_name + "/visualisation/" + "MHChelixPeptid1e.png"
rayTime(differencemap, do_ray)

# Copy above scene but with a pink peptide
pymol.cmd.color("magenta", "ps")
pymol.cmd.util.cnc("ps")


# Save the session
pymol.cmd.save(file_name + "/visualisation/" + file_name + "_peptideMHCvis.pse")

# Quit pymol
pymol.cmd.quit()

if view == True:
    subprocess.call(["pymol", file_name + "/visualisation/" + file_name + "_peptideMHCvis.pse"])


print('     ~  End peptideMHCvisualisation.py v0.1 BETA  ~')

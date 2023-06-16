#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written using Sublime Text & Spyder

@author:    Bruce J MacLachlan, Division of Infection & Immunity, Cardiff University
            ORCID ID: 0000-0002-2685-2733
            Github: https://github.com/brucemaclachlan
            email: maclachlanb@cardiff.ac.uk
            
            omit_maps.py
            Omit map analysis via Refmac or phenix and the subsequent visualisation.
"""

def pMHC_omit_maps_refmac(pdb, mtz, mhc_class, chains_dict, output_dir, file_name):
    import urllib
    import os
    import shutil
    import subprocess
    import sys

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

    if not os.path.exists(output_dir + "/pdbs"):
        print("Creating Directory " + output_dir + "/pdbs")
        os.makedirs(output_dir + "/pdbs")


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

    return edmap, diffmap

def pMHC_omit_maps_phenix(pdb, mtz, mhc_class, chains_dict, output_dir, file_name):
    import urllib
    import os
    import shutil
    import subprocess
    import sys

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

    if not os.path.exists(output_dir + "/pdbs"):
        print("Creating Directory " + output_dir + "/pdbs")
        os.makedirs(output_dir + "/pdbs")


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

    # delete chain using phenix.pdbtools
    command = " ".join(["phenix.pdbtools", pdb, 'remove=\"chain '+peptidechain+'\"', 'filename='+output_dir + "/pdbs/" + file_name + "_nopeptide.pdb"])
    print("\n")
    print(command)
    subprocess.call([command], shell=True)

    
    # Run default phenix.refine with sim anneal without peptide in the groove
    model = output_dir + "/pdbs/" + file_name + "_nopeptide.pdb"
    
    command = " ".join(["phenix.refine", mtz, model, 'simulated_annealing=True', "main.number_of_macro_cycles=2", "tls=False", "--overwrite"])
    print("\n")
    print(command)
    subprocess.call([command], shell=True)
    
    #  I don't really have a handle on how phenix cmd line spits outputs so for now I'll just move them to tidy up
    
    shutil.move(file_name + "_nopeptide"+  "_refine_002.def", output_dir + "/pdbs/" + file_name + "_nopeptide"+  "_refine_002.def")
    shutil.move(file_name + "_nopeptide"+  "_refine_001.mtz", output_dir + "/maps/" + file_name + "_nopeptide"+  "_refine_001.mtz")
    shutil.move(file_name + "_nopeptide"+  "_refine_001.log", output_dir + "/pdbs/" + file_name + "_nopeptide"+  "_refine_001.log")
    shutil.move(file_name + "_nopeptide"+  "_refine_001.pdb", output_dir + "/pdbs/" + file_name + "_nopeptide"+  "_refine_001.pdb")
    shutil.move(file_name + "_nopeptide"+  "_refine_001.cif", output_dir + "/pdbs/" + file_name + "_nopeptide"+  "_refine_001.cif")
    shutil.move(file_name + "_nopeptide"+  "_refine_001.geo", output_dir + "/pdbs/" + file_name + "_nopeptide"+  "_refine_001.geo")
    shutil.move(file_name + "_nopeptide"+  "_refine_001.eff", output_dir + "/pdbs/" + file_name + "_nopeptide"+  "_refine_001.eff")

    
    omit_pdb = output_dir + "/pdbs/" + file_name + "_nopeptide_refine_001.pdb"
    omit_mtz = output_dir + "/maps/" + file_name + "_nopeptide_refine_001.mtz"
    
    # Run phenix.mtz2map which spits out ccp4 format of each map coefficient dataset in mtz fike
    
    command = " ".join(["phenix.mtz2map", omit_mtz, pdb])
    print("\n")
    print(command)
    subprocess.call([command], shell=True)
    
    edmap = file_name + "_nopeptide_refine_001_2mFo-DFc.ccp4"
    diffmap = file_name + "_nopeptide_refine_001_mFo-DFc.ccp4"
    nofillmap = file_name + "_nopeptide_refine_001_2mFo-DFc_no_fill_no_fill.ccp4"
    
    #  I don't really have a handle on how phenix cmd line spits outputs so for now I'll just move them to tidy up
    
    shutil.move(edmap, output_dir + "/maps/" + edmap)
    shutil.move(diffmap, output_dir + "/maps/" + diffmap)
    shutil.move(nofillmap, output_dir + "/maps/" + nofillmap)
    
    # we don't need these so remove them to prevent confusion
    os.remove(output_dir + "/maps/" + edmap)
    os.remove(output_dir + "/maps/" + nofillmap)
    
    # Use phenix.mtz2map to make ccp4 format maps of refined map (i.e. DO NOT USE THE OMIT MAP HERE)
    command = " ".join(["phenix.mtz2map", mtz, pdb])
    print("\n")
    print(command)
    subprocess.call([command], shell=True)
    
    edmap = file_name + "_2mFo-DFc.ccp4"
    diffmap = file_name + "_mFo-DFc.ccp4"
    nofillmap = file_name + "_2mFo-DFc_no_fill_no_fill.ccp4"
    
    #  I don't really have a handle on how phenix cmd line spits outputs so for now I'll just move them to tidy up
    
    shutil.move(edmap, output_dir + "/maps/" + edmap)
    shutil.move(diffmap, output_dir + "/maps/" + diffmap)
    shutil.move(nofillmap, output_dir + "/maps/" + nofillmap)
    
    
    # we don't need these so remove them to prevent confusion
    os.remove(output_dir + "/maps/" + nofillmap)
    
    # COMMENT OUT UP TO HERE TO SPEED UP #
    # This is the important bit.    The edmap must be from the OG .mtz
    #                               The diffmap is from the omit map. This is poor nomenclature and needs to be cleaned up.
    edmap = output_dir + "/maps/" + file_name + "_2mFo-DFc.ccp4"
    diffmap = output_dir + "/maps/" + file_name + "_nopeptide_refine_001_mFo-DFc.ccp4"

    return edmap, diffmap

def pMHC_omit_maps_visualise(pdb, edmap, diffmap, mhc_class, ray, chains_dict, output_dir, file_name, peptide_name):
    import pymol

    import bin.data.viewSet as viewSet
    import bin.data.colourSet as colourSet
    import bin.data.peptide_colours as peptide_colours

    import bin.colour_functions as colour_functions
    import bin.ray_functions as ray_functions
    import bin.image_functions as image_functions
    # Sort chains
    MHCachain, MHCbchain, peptidechain = chains_dict["MHCA"] , chains_dict["MHCB"], chains_dict["PEPTIDE"]
    pdb_name = file_name.lower()

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

    pymol.cmd.load(pdb, "structure")
    pymol.cmd.load(edmap, file_name + "_map")
    pymol.cmd.load(diffmap, file_name + "_dmap")

    pymol.cmd.remove("structure and NOT chain "+chains_dict["MHCA"]+","+chains_dict["MHCB"]+","+chains_dict["PEPTIDE"])
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


    peptide_colour = colour_functions.hex_to_fraction(peptide_colours.peptide_colours[peptide_name])
    colour_functions.set_new_colour(peptide_name, peptide_colour)

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
    pymol.cmd.isomesh("p_map_1sigma", file_name + "_map", 1.0, "ps", carve=2.0)
    pymol.cmd.isomesh("p_map_05sigma", file_name + "_map", 0.5, "ps", carve=2.0)
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
        ray_functions.rayTime(PeptideEdm1sig, do_ray=ray)

    pymol.cmd.hide("mesh", "all")

    pymol.cmd.show("mesh", "p_map_05sigma")
    pymol.cmd.color("grey50", "p_map_05sigma")


    pymol.cmd.set_view(viewSet.pMHC_2)
    # Photo op here
    pymol.cmd.scene(key="PeptideEdm05sig", action="store")
    if ray:
        PeptideEdm05sig = output_dir + "/visualisation/" + "PeptideEdm05sig.png"
        ray_functions.rayTime(PeptideEdm05sig, do_ray=ray)

    # Make difference map 
    pymol.cmd.hide("mesh", "all")
    pymol.cmd.isomesh("posdiffmesh", file_name + "_dmap", 3.0, "ps", carve=2.0)
    pymol.cmd.color("green", "posdiffmesh")
    pymol.cmd.show("mesh", "posdiffmesh")
    pymol.cmd.isomesh("negdiffmesh", file_name + "_dmap", -3.0, "ps", carve=2.0)
    pymol.cmd.color("red", "negdiffmesh")
    pymol.cmd.show("mesh", "negdiffmesh")

    # Photo op here
    pymol.cmd.set_view(viewSet.pMHC_2)
    pymol.cmd.scene(key="differencemap", action="store")
    if ray:
        differencemap = output_dir + "/visualisation/" + "differencemap.png"
        ray_functions.rayTime(differencemap, do_ray=ray)

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
        ray_functions.rayTime(MHChelixPeptide1, do_ray=ray)

    # Copy above scene but with a pink peptide
    pymol.cmd.color("magenta", "ps")
    pymol.cmd.util.cnc("ps")

    # Photo op here
    pymol.cmd.set_view(viewSet.pMHC_2)

    pymol.cmd.scene(key="MHChelixPeptide2", action="store")
    if ray:
        MHChelixPeptide2 = output_dir + "/visualisation/" + "MHChelixPeptide2.png"
        ray_functions.rayTime(MHChelixPeptide2, do_ray=ray)


    # Copy above scene but with a pink peptide
    pymol.cmd.color("magenta", "ps")
    pymol.cmd.util.cnc("ps")

    # Photo op here
    pymol.cmd.set_view(viewSet.pMHC_2)

    pymol.cmd.scene(key="MHChelixPeptide3", action="store")
    if ray:
        MHChelixPeptide3 = output_dir + "/visualisation/" + "MHChelixPeptide3.png"
        ray_functions.rayTime(MHChelixPeptide3, do_ray=ray)

    # Save the session
    pymol.cmd.save(output_dir + "/visualisation/" + file_name + "_peptideMHCvis.pse")

    print('     ~  End omit_maps.pMHC_omit_maps_visualise  ~')

    return None

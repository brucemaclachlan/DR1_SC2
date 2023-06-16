#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written using Sublime Text & Spyder

@author:    Bruce J MacLachlan, Division of Infection & Immunity, Cardiff University
            ORCID ID: 0000-0002-2685-2733
            Github: https://github.com/brucemaclachlan
            email: maclachlanb@cardiff.ac.uk
            
            pymol_functions.py

            This contains pymol functions which I commonly use that are associated with peptide-HLA structures.
"""

def initialisePymol():
    '''
    This function asks python to start a new pymol session and apply a set of parameters related pymol renders the molecules.
    i.e. I don't like shadows, so they are turned off.
    This helps to keep all figures consistent.
    '''
    import pymol
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

def find_nearest_residue(structure, coords, i=0.1, max=3.0):
    '''
    Input a structure name, coords (3 value tuple) and output name and it returns the residue closest to those atoms.
    The coords should be near the ca.
    
    '''
    import pymol
    from Bio.PDB.PDBParser import PDBParser
    
    output = {}
    
    distance = 0.1
    atoms_found = 0
    pymol.cmd.pseudoatom("coord_pseudo", pos=coords)
    
    out_id = "out_id"
    
    while atoms_found == 0:
        atoms_found = pymol.cmd.select(out_id, structure+" within "+str(distance)+" of coord_pseudo")
        if atoms_found > 1:
            
            pymol.stored.list=[]
            pymol.cmd.iterate("("+out_id+")","stored.list.append((resi,resn))")
            found = pymol.stored.list
            
            resis = []
            for atom in found:
                resis.append(atom[0])
            
            if all(x == resis[0] for x in resis) == True:
                None
            else:
                raise Exception("Error in find_nearest_residue:", structure, "distance was:", distance, atoms_found, coords, found)

        distance += i
      
    pymol.stored.list=[]
    pymol.cmd.iterate("("+out_id+")","stored.list.append((resi,resn))")
    resi = pymol.stored.list[0][0]
    resn = pymol.stored.list[0][1]
    
    if distance < max:
        if resn != "HOH":
            output["chain"] = pymol.cmd.get_chains(out_id)[0]
            output["resnum"] = int(resi)
            output["resid"] = resn
        else:
            output["chain"] = None
            output["resnum"] = None
            output["resid"] = None
            
    else:
        output["chain"] = None
        output["resnum"] = None
        output["resid"] = None
    
    pymol.cmd.delete("coord_pseudo")
    pymol.cmd.delete(out_id)
    
                    
    return output

def pMHC_load(structure, file_path):
    import pymol
    pymol.cmd.load(file_path, structure)
    # Let's remove hydrogens from all structures (if they have them) otherwise they can mess with alignments even if we don't see them
    pymol.cmd.remove("hydro")
    return None

def pMHC_align(mobile, target, mhc_class, chains_dict):
    
    import pymol
    # This aligns the a1 and a2 helices of our moving (mobile) and reference (target) MHC molecules

    #class I
    if mhc_class == "I":
        pymol.cmd.align(mobile+" and chain "+chains_dict["MHCA"]+" and resi 137-181,49-85", target+" and chain A and resi 137-181,49-85")

    #class II
    if mhc_class == "II":
        pymol.cmd.align(mobile+" and chain "+chains_dict["MHCA"]+" and resi 46-78 or "+mobile+" and chain "+chains_dict["MHCB"]+" and resi 54-91", target+" and chain A and resi 46-78 or "+target+" and chain B and resi 54-91")

    return None
        
def pMHC_create_objects(structure, mhc_class, chains_dict):

    import pymol
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
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
def wait4ray(query):  
    counter = 0
    print("Did I make it in here?!?")
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
        print("before wait 5")
        time.sleep(5)
        print("after wait 5")
        pymol.cmd.png(saveas,ray=do_ray,width=3000,height=3000, dpi=200)
        print("after cmd.png")
        wait4ray(saveas) 
        print("Done! "+str(saveas)+ " was outputted" )
        return None

def initialisePymol():
    '''
    This function asks python to start a new pymol session and apply a set of parameters related pymol renders the molecules.
    i.e. I don't like shadows, so they are turned off.
    This helps to keep all figures consistent.
    '''
    print("\nInitialising pymol...\n")
    import pymol
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

def magickray(layers, saveas, view, pdb, do_ray):
    
    if do_ray == 1:
        
        pymol_dir = pdb.split(".")[0]+'/pymol'
        
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
        subprocess.call(["convert", pymol_dir+"/"+saveas+"/"+saveas+"*", "-background", "none", "-flatten", pymol_dir+"/"+saveas+".png"])
        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Done! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        
    pymol.cmd.hide("everything", "all")
    
    for layer in layers:
        for ob in layer:
            for representation in layer[ob]:
                pymol.cmd.show(representation, ob)
    
    pymol.cmd.set_view(view)
    pymol.cmd.scene(key=saveas, action="store")
    return None   

def pMHC_load(structure, file_path):
    import pymol
    pymol.cmd.load(file_path, structure)
    # Let's remove hydrogens from all structures (if they have them) otherwise they can mess with alignments even if we don't see them
    pymol.cmd.remove("hydro")
    return None

def pMHC_align(mobile, target, mhc_class):
    import pymol
    # This aligns the a1 and a2 helices of our moving (mobile) and reference (target) MHC molecules

    #class I
    if mhc_class == "I":
        pymol.cmd.align(mobile+" and chain A and resi 137-181,49-85", target+" and chain A and resi 137-181,49-85")

    #class II
    if mhc_class == "II":
        pymol.cmd.align(mobile+" and chain A and resi 46-78 or "+mobile+" and chain B and resi 54-91", target+" and chain A and resi 46-78 or "+target+" and chain B and resi 54-91")

    return None
        
def pMHC_create_objects(structure, mhc_class):
    import pymol

    if mhc_class == 'I':

        pymol.cmd.select(structure+"_MHCa1a2", selection=structure+" and chain A and resi 137-181,49-85")
        pymol.cmd.create(structure+"_MHCa1a2_obj", selection=structure+"_MHCa1a2")

        pymol.cmd.select(structure+"_MHCa1", selection=structure+" and chain A and resi 49-85")
        pymol.cmd.create(structure+"_MHCa1_obj", selection=structure+"_MHCa1")

        pymol.cmd.select(structure+"_MHCa2", selection=structure+" and chain A and resi 137-181")
        pymol.cmd.create(structure+"_MHCa2_obj", selection=structure+"_MHCa2")
        
        pymol.cmd.select(structure+"_MHCgroove", selection=structure+" and chain A and resi 1-136")
        pymol.cmd.create(structure+"_MHCgroove_obj", selection=structure+"_MHCgroove")
      

    if mhc_class == 'II':

        pymol.cmd.select(structure+"_MHCa1a2", selection=structure+" and chain A and resi 46-78 or "+structure+" and chain B and resi 54-91")
        pymol.cmd.create(structure+"_MHCa1a2_obj", selection=structure+"_MHCa1a2")

        pymol.cmd.select(structure+"_MHCa1", selection=structure+" and chain A and resi 46-78")
        pymol.cmd.create(structure+"_MHCa1_obj", selection=structure+"_MHCa1")

        pymol.cmd.select(structure+"_MHCa2", selection=structure+" and chain B and resi 54-91")
        pymol.cmd.create(structure+"_MHCa2_obj", selection=structure+"_MHCa2")
        
        pymol.cmd.select(structure+"_MHCgroove", selection=structure+" and chain A and resi 1-78 or "+structure+" and chain B and resi 6-91")
        pymol.cmd.create(structure+"_MHCgroove_obj", selection=structure+"_MHCgroove")


    pymol.cmd.select(structure+"_p", selection=structure+" and chain C")
    pymol.cmd.create(structure+"_p_obj", selection=structure+"_p")

    pymol.cmd.delete("hetatms")

    return None

def magick_outline(layer, view, image_name, subimage_name, dir, do_ray):
    '''
    Hacky way to generate an outline of an object then overlay it on the image stack.
    subimage_name should be image_name+i (can only do layers on top).
    This code isn't very flexible right now.
    '''
    import subprocess
    if do_ray == True:
        
        pymol_dir = pdb.split(".")[0]+'/pymol'
        
        for lay in layer:
            selection = lay
            representation = layer[lay][0]
        
        pymol.cmd.hide("everything", "all")
        pymol.cmd.create(selection+"o",selection)
        
        pymol.cmd.show(representation,selection+"o")
        pymol.cmd.set_view(view)
        
        rayTime(pymol_dir+"/"+image_name+"/"+subimage_name+".png", do_ray)
        time.sleep(2)
        subprocess.call(["convert", pymol_dir+"/"+image_name+"/"+subimage_name+".png", "-alpha", "extract", "-threshold", "0", "-negate", "-transparent", "white", pymol_dir+"/"+image_name+"/"+subimage_name+"_1.png"], shell=True)
        subprocess.call(["convert", pymol_dir+"/"+image_name+"/"+subimage_name+"_1.png", "-negate", "-threshold", "1", "-edge", "10", pymol_dir+"/"+image_name+"/"+subimage_name+"_2.png"], shell=True)
        subprocess.call(["convert", pymol_dir+"/"+image_name+"/"+subimage_name+"_2.png", "-negate", pymol_dir+"/"+image_name+"/"+subimage_name+"_3.png"], shell=True)
        subprocess.call(["convert", pymol_dir+"/"+image_name+"/"+subimage_name+"_3.png", "-fuzz", "20%", "-transparent", "white", pymol_dir+"/"+image_name+"/"+subimage_name+"_4.png"], shell=True)
        
        os.remove(pymol_dir+"/"+image_name+"/"+subimage_name+".png")
        os.remove(pymol_dir+"/"+image_name+"/"+subimage_name+"_1.png")
        os.remove(pymol_dir+"/"+image_name+"/"+subimage_name+"_2.png")
        os.remove(pymol_dir+"/"+image_name+"/"+subimage_name+"_3.png")
        
        subprocess.call(["convert", pymol_dir+"/"+image_name+"/*.png", "-background", "none", "-flatten", pymol_dir+"/"+image_name+".png"], shell=True)
        
    return None

def apbs_electrostatics(pdb_name, pdb, do_apbs, do_ray, template, mhc_class, neg_colour, pos_colour):
    import subprocess
    import os
    import pymol

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
    pymol.cmd.delete('align_target')
    pymol.cmd.remove("hetatm")

    pymol.cmd.select("pMHC", "structure and chain A+B+C")
    pymol.cmd.create("pMHC"+"o", "pMHC") 
    pMHC_create_objects("structure", mhc_class)
    
   

    obj_name = "pMHCo"
    pqr_name = apbs_dir+"/"+pdb_name+"_"+obj_name+".pqr"
    apbs_output_name = pdb_name+"_"+obj_name+"_map"
    apbs_output_name_full = apbs_dir+"/"+pdb_name+"_"+obj_name+"_map"
    apbs_input_name = pdb_name+"_"+obj_name+"_apbs.in"
    prepped_pdb = apbs_dir+"/"+pdb_name+"_"+obj_name+".pdb"

    pymol.cmd.save(prepped_pdb, obj_name)
    
    if do_apbs == True:
        
        print("Converting pdb file of object to .pqr using subprocess.call to p2b2pqr")
        command = " ".join(["pdb2pqr30", prepped_pdb, pqr_name,"--ff=AMBER"])
        print(command)
        subprocess.call([command], shell=True, stdout=logger, stderr=logger)

        print("######## END pdb2pqr30 ##########")

        with open("bin/data/apbs_template.txt", 'r') as file :
            filedata = file.read()
        
        filedata = filedata.replace('PQR_IN', os.path.realpath(pqr_name))
        filedata = filedata.replace('APBS_OUT', os.path.realpath(apbs_output_name))
        
        with open(apbs_input_name, 'w') as file:
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
    pymol.cmd.set("surface_color", "map_esramp", "pMHCo")    
    pymol.cmd.set("surface_ramp_above_mode", 1)

    pymol.cmd.show("surface" , "pMHCo")
    pymol.cmd.colour("grey90" , "pMHCo")
    pymol.cmd.colour("grey90" , "structure_peptide_obj")

    import bin.data.viewSet as viewSet
    layer1 = {"pMHCo": ["surface"]
             }
    layer2 = {"structure_MHCa1a2": ["cartoon"],
              "structure_peptide_obj" : "sticks"
             }

    layers = [layer1, layer2]

    # MHC birdseye view of electrostatics
    pymol.cmd.disable("*_esramp")
    pymol.cmd.set("cartoon_transparency", 0.7, "mMHCa1a2")
    magickray(layers, "electrostatic_pMHC", viewSet.birdseye, pdb, ray)

    spike1 = {"p_obj": ["surface"]}
    electrostatics.magick_outline(spike1, viewSet.birdseye, "electrostatic_pMHC", "melectrostatic_pMHC3", pdb_name+"/apbs/", ray)
    spike2 = {"structure_MHCgroove and chain A": ["surface"]}
    electrostatics.magick_outline(spike2, viewSet.birdseye, "electrostatic_pMHC", "melectrostatic_pMHC4", pdb_name+"/apbs/", ray)
    spike3 = {"structure_MHCgroove and chain A": ["surface"]}
    electrostatics.magick_outline(spike3, viewSet.birdseye, "electrostatic_pMHC", "melectrostatic_pMHC5", pdb_name+"/apbs/", ray)

    pymol.cmd.save(pdb_name+"/apbs/"+pdb_name+"_electrostatics.pse")

    return None
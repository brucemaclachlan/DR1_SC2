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
import os
import pymol
import time
import itertools
from PIL import Image

description=\
"foobar"

def parse_args():
    import argparse
    """parse arguments for command line"""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--input', dest = 'fasta', type = str, required = True, help ='The fasta file containing sequence of peptide to be drawn')

    args = parser.parse_args()
    return args





def generate_peptide(sequence):

    from PeptideBuilder import Geometry
    from PeptideBuilder import PeptideBuilder

    geo = Geometry.geometry(sequence[0])
    structure = PeptideBuilder.initialize_res(geo)

    for residue in sequence[1:]:
        geo = Geometry.geometry(residue)
        PeptideBuilder.add_residue(structure, geo)

    PeptideBuilder.add_terminal_OXT(structure)

    return structure

def save_modelled_peptide(structure, name):

    import Bio.PDB
    out = Bio.PDB.PDBIO()
    out.set_structure(structure)
    out.save(name)

    return None

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

def main():

    # Handle command line arguments
    args = parse_args()
    print("Running analysis with the following inputs.. ")
    print(args)
    fasta = args.fasta
    initialisePymol()
    output_dir = "drawn_peptides/"


    if not os.path.exists(output_dir ):
        print("Creating Directory "+output_dir )
        os.makedirs(output_dir )

    from Bio import SeqIO

    for record in SeqIO.parse(fasta, "fasta"):
        print(record.seq)
        print(record.id)

        pdb_file_name = output_dir+record.id+".pdb"
        pdb_name = record.id

        structure = generate_peptide(record.seq)
        save_modelled_peptide(structure, name=pdb_file_name)



        pymol.cmd.load(pdb_file_name, "peptide")


        pymol.cmd.hide("everything", "all")
        pymol.cmd.show("cartoon" ,  "peptide")
        pymol.cmd.show("sticks" ,  "peptide and sc. or peptide and name ca")
        pymol.cmd.color("yellow" , "peptide")
        pymol.util.cnc("peptide")

        pymol.cmd.set("cartoon_loop_radius",0.4)
        pymol.cmd.set("stick_radius" , 0.4)


        pymol.cmd.set_view ("\
     0.800060689,    0.550701916,    0.237970576,\
    -0.585625112,    0.802998960,    0.110611588,\
    -0.130176201,   -0.227857754,    0.964953244,\
     0.000000000,    0.000000000, -256.088378906,\
    19.274999619,  -14.034000397,   -0.952000022,\
  -2223.210205078, 2735.386718750,   20.000000000" )


        image_file_name = output_dir+pdb_name+".png"
        rayTime(image_file_name, 1)

        
        pymol.cmd.save("test.pse")

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


        flat_image_file_name = output_dir+pdb_name+"PILtoon.png"

        peptide_blob = image_not_transparent_filter(image_file_name, cell_press_colours["cell_3"]["rgb"]).save(flat_image_file_name)


        if record.id in ["X2_A2_omicron", "SC2_6_omicron"]:
            pymol.cmd.hide("everything" , "all")

            if record.id == "X2_A2_omicron":

                pymol.cmd.show("sticks" ,  "peptide and resi 8,11,13 and sc. or peptide and resi 8,11,13 and name ca")

            if record.id == "SC2_6_omicron":

                 pymol.cmd.show("sticks" ,  "peptide and resi 4 and sc. or peptide and resi 4 and name ca")

            omicron_spike_in = output_dir+pdb_name+"_omicron_spike_in.png"
            rayTime(omicron_spike_in, 1)

            omicron_spike_in_flat = output_dir+pdb_name+"_omicron_spike_in_flat.png"
            peptide_blob_omicron = image_not_transparent_filter(omicron_spike_in, (198,99,143)).save(omicron_spike_in_flat)


            omicron_flat_image_file_name = output_dir+pdb_name+"PILtoon_omicron.png"

            layer_1 = flat_image_file_name
            layer_2 = omicron_spike_in_flat

            layers = [layer_1, layer_2]

            recursive_alpha_composite(layers, (3000,3000), omicron_flat_image_file_name)


        pymol.cmd.delete("all")








if __name__ == '__main__':
    main()
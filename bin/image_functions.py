#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written using Sublime Text & Spyder

@author:    Bruce J MacLachlan, Division of Infection & Immunity, Cardiff University
            ORCID ID: 0000-0002-2685-2733
            Github: https://github.com/brucemaclachlan
            email: maclachlanb@cardiff.ac.uk
            
            Image functions
"""

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
    from PIL import Image
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
    from PIL import Image
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
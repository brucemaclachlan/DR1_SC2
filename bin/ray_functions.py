#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written using Sublime Text & Spyder

@author:    Bruce J MacLachlan, Division of Infection & Immunity, Cardiff University
            ORCID ID: 0000-0002-2685-2733
            Github: https://github.com/brucemaclachlan
            email: maclachlanb@cardiff.ac.uk
            
            ray_functions

            This contains a few functions a use frequently to render/ray images in pymol
"""
def wait4ray(query): 
    import sys
    import time 
    import itertools
    import os
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
    import pymol
    import time
    import pymol
    import os

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

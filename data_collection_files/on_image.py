#!/usr/bin/env python
import os,sys,pdb,scipy,glob
import pandas as pd
from pylab import *


from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u


def get_coordinates(fitsfile,rr=None,dd=None):
    
    with fits.open(fitsfile) as hdul:
        # Get the WCS object
        w = WCS(hdul[1].header)
        coords = SkyCoord(rr, dd, unit='deg')
        # Convert to pixel coordinates
        x, y = w.world_to_pixel(coords)
        if ((x <= 0) |
            (y <= 0) |
            (x > w.pixel_shape[0]) |
            (x > w.pixel_shape[0])
            ): 
            print(f"Filename: %s, Pixel coordinates: x=%.2f, y=%.2f (off axis)" %(fitsfile, x, y))
        else:
            print(f"Filename: %s, Pixel coordinates: x=%.2f, y=%.2f" %(fitsfile, x, y))
    return(x, y)

if __name__=='__main__':
    
    filename = sys.argv[1]
    rr = sys.argv[2]
    dd = sys.argv[3]

    
    
    if filename.startswith('@'):
        filelist = filename.strip('@')
        f = open(filelist,'r')
        lines = f.readlines()
        file_list = [line.strip() for line in lines]
        f.close()
    else:
        file_list = [filename]

    for filename in file_list:
        filename=filename.rstrip()
        filename=filename.strip('*')
        get_coordinates(filename, rr=rr, dd=dd)
        
        

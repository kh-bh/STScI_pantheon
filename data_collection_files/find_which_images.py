#!/usr/bin/env python
'''
this script will query data on disk to determine
(a) if SN coordinates are on image, and
(b) if there are any counts in a small box around that position
the goal would be to return list of images to quick-look to determine
if we can get clean photometry from them.
'''

import os,sys,pdb,scipy,glob
import pandas as pd
from pylab import *

from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u


def get_coordinates(fitsfile,rr=None,dd=None, verbose=False):
    img = os.path.basename(fitsfile)
    with fits.open(fitsfile) as hdul:
        # Get the WCS object
        w = WCS(hdul[1].header)
        coords = SkyCoord(rr, dd, unit='deg')
        # Convert to pixel coordinates
        try:
            x, y = w.world_to_pixel(coords)
            print(f"Successfully converted: RA={rr}, Dec={dd} -> x={x}, y={y}")
        except Exception as e:
            print(f"WCS conversion failed for {img}: {e}")
            print(f"RA={rr}, Dec={dd}")
            status = False
            return(0, 0, status)
        if ((x <= 0) |
            (y <= 0) |
            (x > w.pixel_shape[1]) |  # x should be compared to width (shape[1])
            (y > w.pixel_shape[0])    # y should be compared to height (shape[0])
            ): 
            if verbose: print(f"Filename: %s, Pixel coordinates: x=%.2f, y=%.2f (off axis)" %(img, x, y))
            status = False
        else:
            if verbose: print(f"Filename: %s, Pixel coordinates: x=%.2f, y=%.2f" %(img, x, y))
            status = True
    return(x, y, status)

def get_counts(fitsfile,xx=None,yy=None, box=100):
    image_data = fits.getdata(fitsfile, ext=1)
    region = image_data[int(yy-0.5*box):int(yy+0.5*box),  # Note: y is first dimension
                        int(xx-0.5*box):int(xx+0.5*box)]  # Note: x is second dimension
    return(region.sum())



if __name__=='__main__':

    #datatable = 'pantheon_fits_cleand_hst_data.csv'
    datatable = 'data_collection_files/data_files/fits_cleaned_JWST_data.csv'
    df = pd.read_csv(datatable)
    verbose =False

    ## coordtable = 'probable_galaxy_catalog.csv'
    coordtable = 'data_collection_files/data_files/unmatched_pantheon.csv'
    tf = pd.read_csv(coordtable)
    
    for index, rw in df.iterrows():
        status = False

        snid = rw.SNID
        if verbose: print(snid, index)
        try:
            ## rr = tf[tf.SNID==snid].RA_Galaxy.iloc[0]
            ## dd = tf[tf.SNID==snid].Dec_Galaxy.iloc[0]
            rr = tf[tf.SNID==snid].RA.iloc[0]
            dd = tf[tf.SNID==snid].Dec.iloc[0]
        except:
            print('Help!')
            pdb.set_trace()

        filename = rw.filename
        x,y,status=get_coordinates(filename, rr=rr, dd=dd, verbose=False)
        df.loc[index, 'On_Image']=False
        if (status):
            counts = get_counts(filename, xx=x, yy=y)
            if ((counts > 0) & (not np.isnan(counts))):
                df.loc[index, 'On_Image']=True
                

    of = df[df.On_Image==True]                
    of=of.reset_index()                        
    #outab = 'out_table.csv'
    # outab = 'data_collection_files/JWST_table_sn_in_images.csv'
    # if os.path.isfile(outab): os.remove(outab)
    # of.to_csv(outab)

    #merge of table with 'data_collection_files/data_files/reduced_table_sn_in_images.csv' table

    
    #load the HST_table_sn_im_checklist.csv and merge with the out_table.csv
    sn_in_images = pd.read_csv('data_collection_files/data_files/table_sn_in_images.csv')
    all_sn = pd.read_csv('data_collection_files/data_files/HST_table_sn_im_checklist.csv')
    #merge the two tables on the SNID column
    #merged_table = pd.merge(sn_in_images, all_sn, on='filename', how='left')
    
    #save the merged table
    #merged_table.to_csv('data_collection_files/data_files/reduced_table_sn_in_images.csv', index=False)

    #load the reduced_table_sn_in_images.csv
    merged_table = pd.read_csv('data_collection_files/data_files/reduced_table_sn_in_images.csv')
    # Reset index to ensure sequential indices (0, 1, 2, 3...)
    merged_table = merged_table.reset_index(drop=True)
    #chnage the name of the colums that end with _x to the original name without _x
    merged_table = merged_table.rename(columns={col: col.replace('_y', '') for col in merged_table.columns if col.endswith('_y')})
    #drop the columns that end with _x or _y
    
    merged_table = merged_table.drop(columns=[col for col in merged_table.columns if col.endswith('_x') or col.endswith('_y')])

    
    # Concatenate the cleaned table with the 'of' table (add of at the bottom)
    jwst_hst_table = pd.concat([merged_table, of], ignore_index=True)
    # Append to the existing file instead of overwriting
    jwst_hst_table.to_csv('data_collection_files/data_files/reduced_table_sn_in_images.csv', index=False, mode='a')



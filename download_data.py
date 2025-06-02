import sys,os,glob
from astropy.io import fits
from astropy.table import Table
from astropy.nddata import extract_array
from astropy.coordinates import SkyCoord
from astropy import wcs
from astropy.wcs.utils import skycoord_to_pixel
from astropy import units as u
import numpy as np
import matplotlib.pyplot as plt
from astroquery.mast import Observations
from astropy.visualization import (simple_norm,LinearStretch)
import pandas as pd

def main():
    df = pd.read_csv('unmatched_pantheon.csv')
    vp = pd.read_csv('valid_points_index_list.csv')
    print(vp.dtypes)

    def save_images_supernova(resolved_coord_index):
        """
        If there are images for the supernovae in JWST or HST at the resolved query, create a folder and keep the data there
        Parameters:
            resolved_coord_index: integer, the index number of the query
        Returns:
            Folder of Data at requested file path (pantheon_data_folder/{resolved coord})
        """
        print(resolved_coord_index, df['resolved_coord'][resolved_coord_index])

        #try the the query, if a error is thrown then print "No Data Points" and "skip" this data point, else make folder and store data
        try:
            obs_table = Observations.query_criteria(coordinates=df['resolved_coord'][resolved_coord_index],
                                            radius="0.006 deg",
                                            intentType = 'science',
                                            filters = ['F1*'],
                                            obs_collection=['HST'])
            obs_table = obs_table[obs_table['calib_level']==3]
            data_products = Observations.get_product_list(obs_table)
            data_products = data_products[data_products['calib_level'] == 3]
            data_products = data_products[data_products['productType'] == 'SCIENCE']
        except: 
            print("No Data Points :(")
            exit
        else:
            print("Has Data Points!")
            folder_path = "pantheon_data_folder/{}".format(df['SNID'][resolved_coord_index])
            os.makedirs(folder_path, exist_ok=True)
            manifest = Observations.download_products(data_products, download_dir=folder_path, extension=['fits'])
            print(manifest)
    
    for i in range(len(vp['0'])):
        save_images_supernova(vp['0'][i])

if __name__=="__main__":
    main()
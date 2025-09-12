#!/usr/bin/env python
import subprocess
import os
import pandas as pd
import shutil
import sys

def source_extractor(source_extractor_path = ""):
    """
    Params:
        - source_extractor_path: the path to where the source extractors files are located
    """
    if not os.path.exists(source_extractor_path + "default.conv"):
        print("HEY MAKE A DEFAULT CONV FILE")

    if not os.path.exists(source_extractor_path + "default.param"):
        subprocess.call("sex -dp -> default.param", shell=True)
        print("HEY CHANGE THE DEFAULT PARAMETER FILE")
    
    if not os.path.exists(source_extractor_path + "default.sex"):
        subprocess.call("sex -d -> default.sex", shell=True)
        print("HEY MAKE THE DEFAULT SEX FILE")

def catalog_creator(SNID, filename, base, filter, home_path, sex_file_name = 'default', source_path = "", force=False):
    #Check if the filter is of type detection (full stack image). Don't go through these
    if filter == 'detection':
        return;
    #Check if source extractor was already done. Don't do these (unless force=True)
    if not force:
        done = check_file_exist(home_path, source_path, SNID, base)
        if done == True:
            return;

    dest_path = filename
    # if the instrument is JWST, then we need to specify extension 1
    if base.startswith('jw'):
        exec_string = 'sex ' + dest_path + '[1] -c ' + sex_file_name + '.sex'  # [1] specifies extension 1
    else:
        exec_string = 'sex ' + dest_path + ' -c ' + sex_file_name + '.sex'

    
    cat_file_name = "test.cat"  # This matches the CATALOG_NAME in default.sex
    new_cat_file_name = f"{base}.cat"

    fits_file_name = "check.fits"
    new_fits_file_name = f"{base}.fits"

    
    subprocess.call(exec_string, shell = True);

    os.rename(cat_file_name, new_cat_file_name)
    #os.rename(fits_file_name, new_fits_file_name)
    move_files_group(SNID, base, new_cat_file_name, new_fits_file_name, home_path, source_path)

def check_file_exist(home_path, source_path, SNID, base):
    """
    Makes the process able to be stopped and started at your discretion
    Returns:
        Boolean value of if the file exists
    """
    return os.path.exists(f"{home_path}/{source_path}/{SNID}/{base}")
    

        

def move_files_group(SNid, base, cat_file_name, fits_file_name, home_path, source_path = ''):
    """
    Moves files from one location to another given the filename, home_path, and source path
    Parameters:
        SNid: The SNID of the suoernova
        file_name: the name of the file
        home_path: the path to the file (generally the original home_path)
        source_path: the path of the sub folder that should store files (defaults to no folder but causes problems)
    """
    try:
        destination_path = f"{home_path}/{source_path}/{SNid}/{base}"
        os.makedirs(destination_path, exist_ok=True)
        
        # Remove existing file if it exists (to allow overwriting)
        dest_file = os.path.join(destination_path, cat_file_name)
        if os.path.exists(dest_file):
            os.remove(dest_file)
            
        shutil.move(cat_file_name, destination_path)
        #shutil.move(fits_file_name, destination_path)
    except Exception as e:
        print(f"An unexpected error occurred: {e}")




def main():
    full_dataframe = pd.read_csv('/home/epadill/photometry_ds9/STScI_pantheon/data_collection_files/data_files/fits_summary.csv')
    #get the rows where telescope is JWST
    full_dataframe = full_dataframe[full_dataframe['Telescope'] == 'JWST']
    print(full_dataframe['filename'].str.split('/').str[-1])    #for testing

    #activates the conda enviroment
    # Note: conda activation in subprocess doesn't work well
    # Make sure to run this script in the correct environment
    pass

    df = full_dataframe[start:end]
    #Makes sure that each of the source extractor files are avaliable and found
    source_extractor(source_extractor_path= "/home/epadill/photometry_ds9/STScI_pantheon/")
    for row in df.itertuples():
        catalog_creator(SNID = row.SNID, 
                        filename = row.filename, 
                        base = row.File_key,
                        filter = row.Filter, 
                        home_path = "/astro/armin/bhoomika",
                        sex_file_name = sex_name, 
                        source_path= "source_extractor",
                        force=True)

if __name__=="__main__":
    start = int(sys.argv[1])
    end = int(sys.argv[2])
    sex_name = sys.argv[3]
    main()
                         


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

def catalog_creator(SNID, filename, base, filter, home_path, sex_file_name = 'default', source_path = ""):
    #Check if the filter is of type detection (full stack image). Don't go through these
    if filter == 'detection':
        return;
    #Check if source extractor was already done. Don't do these
    done = check_file_exist(home_path, source_path, SNID, base)
    if done == True:
        return;

    dest_path = filename
    exec_string = 'sex ' + dest_path + ' -c ' + sex_file_name + '.sex'

    
    cat_file_name = f"{sex_file_name}.cat"
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
        shutil.move(cat_file_name, destination_path)
        #shutil.move(fits_file_name, destination_path)
    except Exception as e:
        print(f"An unexpected error occurred: {e}")




def main():
    full_dataframe = pd.read_csv('/home/bkhatri/STScI_pantheon/data_collection_files/data_files/fits_summary.csv')

    #activates the conda enviroment
    conda_act = "conda activate .conda"
    subprocess.call("conda init", shell = True)
    subprocess.call(conda_act, shell=True)

    df = full_dataframe[start:end]
    #Makes sure that each of the source extractor files are avaliable and found
    source_extractor(source_extractor_path= "/home/bkhatri/STScI_pantheon/")
    for row in df.itertuples():
        catalog_creator(SNID = row.SNID, 
                        filename = row.filename, 
                        base = row.File_key,
                        filter = row.Filter, 
                        home_path = "/astro/armin/bhoomika",
                        sex_file_name = sex_name, 
                        source_path= "source_extractor")

if __name__=="__main__":
    start = int(sys.argv[1])
    end = int(sys.argv[2])
    sex_name = sys.argv[3]
    main()
                         


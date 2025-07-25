from astropy.io import ascii
import pandas as pd
import numpy as np
import glob
import os


def create_brightest_index(root_dir, relative_key, output_csv):
    """
    Create a datatable with the BRIGHTEST and LARGEST objects of each source extractor file
    Parameters:
        root_dir (str): The root directory to search.

    Returns:
        dict: A dictionary where keys are relative paths to .cat files,
              and values are dictionaries of the first data row.
    """
    # Find all .cat files recursively
    
    cat_files = glob.glob(os.path.join(root_dir, '**', '*.cat'), recursive=True)

    data_dict = []

    for file_path in cat_files:
        try:
            starcounts = ascii.read(file_path, data_start=0)

            strPath = os.path.realpath(file_path)
            nmFolders = strPath.split(os.path.sep)

            SNID_key = nmFolders[-3]

            b = os.path.basename(file_path)
            filename_without_extension = os.path.splitext(b)[0]
            relative_key += 1

            info = {
                    "SNID": SNID_key,
                    "File_key": filename_without_extension, 
                    "Galaxy_Key": relative_key,
                    "RA_Galaxy": starcounts['ALPHA_J2000'][0],
                    "Dec_Galaxy": starcounts['DELTA_J2000'][0],
                    "Position_Angle": starcounts['THETA_IMAGE'][0],
                    "Major_Axis": starcounts['A_IMAGE'][0],
                    "Minor_Axis": starcounts['B_IMAGE'][0],
                    "Mag": starcounts['MAG_AUTO'][0],
                    "Mag_Error": starcounts['MAGERR_AUTO'][0],
                    "Flux": starcounts['FLUX_AUTO'][0],
                    "Flux_Error": starcounts['FLUXERR_AUTO'][0],
                    "CXX": starcounts['CXX_IMAGE'][0],
                    "CXY": starcounts['CXY_IMAGE'][0],
                    "CYY":starcounts['CYY_IMAGE'][0]
                    }
            data_dict.append(info)

        except Exception as e:
            print(f"Error reading {file_path}: {e}")
        
        # Convert to DataFrame and export
    df = pd.DataFrame(data_dict)
    df.to_csv(output_csv, index=False)
    print(f"CSV exported: {output_csv}")
    return data_dict





if __name__== "__main__":
    source_extractor_path = "/astro/armin/bhoomika/source_extractor"
    #root_dir = os.path.join(home_dir, 'source_extractor')
    #print(root_dir)
    relative_key = 0
    create_brightest_index(source_extractor_path, relative_key, "brightest_galaxy.csv")
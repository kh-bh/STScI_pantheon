'''
This script queries the optical observations of the supernovae from the SDSS database and merges
the data with the reduced table of the supernovae that are in the images. We also create a file that combines the sdss ad 
the hst IR data into a single file. 

We use the file we create to plot the SDSS data and the hst IR data on the same plot. 

'''

import pandas as pd
import numpy as np
from astroquery.sdss import SDSS
from astropy.table import Table
from tqdm import tqdm
import time
from astropy.coordinates import SkyCoord
from astropy import units as u
import sys
import os

# Add the data_collection_files directory to the path
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'data_collection_files'))

# Import everything from magnitude_corrector_calculator to get all required imports
from magnitude_corrector_calculator import *

def calculate_hst_zeropoint_ab(fits_file_path, verbose=False):
    """
    Calculate AB zeropoint for HST FITS files using multiple methods.
    
    Parameters:
    -----------
    fits_file_path : str
        Path to the FITS file
    verbose : bool
        Whether to print debug information
        
    Returns:
    --------
    float
        AB zeropoint value, or 0 if calculation fails
    """
    try:
        with fits.open(fits_file_path) as hdul:
            header = hdul[0].header
            
            # Extract key information
            filter_name = header.get("FILTER")
            instrument = header.get("INSTRUME")
            telescope = header.get("TELESCOP")
            
            if verbose:
                print(f"\n=== Calculating zeropoint for {fits_file_path} ===")
                print(f"  FILTER: {filter_name}")
                print(f"  INSTRUMENT: {instrument}")
                print(f"  TELESCOPE: {telescope}")
            
            # Method 1: Look for direct zeropoint keywords (expanded list)
            phot_zpt = (header.get("PHOTZPT") or header.get("ZPT") or header.get("ZEROPOINT") or 
                       header.get("MAGZPT") or header.get("MAGZERO") or header.get("PHOTZERO") or
                       header.get("ABMAG") or header.get("AB_ZERO") or header.get("VEGAMAG"))
            
            if phot_zpt is not None and phot_zpt > 0:
                if verbose:
                    print(f"  Using direct zeropoint: {phot_zpt}")
                return phot_zpt
            
            # Method 2: Calculate from PHOTFLAM and PHOTPLAM
            phot_flam = header.get("PHOTFLAM") or header.get("FLAM") or header.get("FLUXZERO")
            phot_plam = header.get("PHOTPLAM") or header.get("PLAM") or header.get("PIVOT")
            
            if verbose:
                print(f"  PHOTFLAM/FLAM: {phot_flam}")
                print(f"  PHOTPLAM/PLAM: {phot_plam}")
            
            if phot_flam is not None and phot_plam is not None and phot_flam > 0 and phot_plam > 0:
                zp_ab = -2.5 * np.log10(phot_flam) - 5 * np.log10(phot_plam) - 2.408
                if verbose:
                    print(f"  Calculated from PHOTFLAM/PHOTPLAM: {zp_ab}")
                return zp_ab

            
            
            # Method 3: Use standard HST zeropoints based on filter
            hst_zeropoints = {
                # WFC3/IR filters
                'F160W': 25.96,  # WFC3/IR F160W
                'F105W': 26.27,  # WFC3/IR F105W  
                'F110W': 26.82,  # WFC3/IR F110W
                'F125W': 26.23,  # WFC3/IR F125W
                'F140W': 26.45,  # WFC3/IR F140W
                
                # WFC3/UVIS filters
                'F275W': 24.13,  # WFC3/UVIS F275W
                'F336W': 24.67,  # WFC3/UVIS F336W
                'F435W': 25.68,  # WFC3/UVIS F435W
                'F475W': 26.16,  # WFC3/UVIS F475W
                'F555W': 25.72,  # WFC3/UVIS F555W
                'F606W': 26.40,  # WFC3/UVIS F606W
                'F625W': 25.69,  # WFC3/UVIS F625W
                'F775W': 25.25,  # WFC3/UVIS F775W
                'F814W': 25.96,  # WFC3/UVIS F814W
                'F850LP': 24.87,  # WFC3/UVIS F850LP
                
                # ACS/WFC filters (if present)
                'F435W': 25.68,  # ACS/WFC F435W
                'F475W': 26.16,  # ACS/WFC F475W
                'F555W': 25.72,  # ACS/WFC F555W
                'F606W': 26.40,  # ACS/WFC F606W
                'F625W': 25.69,  # ACS/WFC F625W
                'F775W': 25.25,  # ACS/WFC F775W
                'F814W': 25.96,  # ACS/WFC F814W
                'F850LP': 24.87,  # ACS/WFC F850LP
            }
            
            if filter_name in hst_zeropoints:
                zp_ab = hst_zeropoints[filter_name]
                if verbose:
                    print(f"  Using standard HST zeropoint for {filter_name}: {zp_ab}")
                return zp_ab
            else:
                if verbose:
                    print(f"  Warning: No zeropoint found for filter {filter_name}")
                return 0.0
                
    except Exception as e:
        if verbose:
            print(f"  Error reading {fits_file_path}: {e}")
        return 0.0

def get_sdss_good_coordinates_all(name, RA, Dec, search_radius=0.1):
    """
    Query SDSS DR16 for galaxies with good quality data near specified coordinates.
    
    Parameters:
    -----------
    name : str
        Name/identifier for the target object
    RA : float
        Right Ascension in degrees
    Dec : float
        Declination in degrees
    search_radius : float, default=0.1
        Search radius in degrees (default 0.1° ≈ 6 arcminutes)
    output_file : str, default='sdss_good_galaxies_test.csv'
        File to save the results
        
    Returns:
    --------
    astropy.table.Table or None
        Table with the retrieved data, or None if no data found
    """
    from astropy.coordinates import SkyCoord
    from astropy import units as u
    
    # Convert coordinates to SkyCoord for distance calculations
    target_coord = SkyCoord(ra=RA*u.degree, dec=Dec*u.degree, frame='icrs')
    
    # Calculate search boundaries
    ra_min = RA - search_radius
    ra_max = RA + search_radius
    dec_min = Dec - search_radius
    dec_max = Dec + search_radius
    
    # Handle RA wraparound near 0/360 degrees
    if ra_min < 0:
        ra_min += 360
    if ra_max > 360:
        ra_max -= 360
    
    print(f"Querying SDSS DR16 for galaxies near {name} at RA={RA:.6f}°, Dec={Dec:.6f}°")
    print(f"Search radius: {search_radius}° ({search_radius*60:.1f} arcminutes)")
    
    # Use SDSS cone search instead of SQL query
    

    try:
        # Create coordinate object
        coord = SkyCoord(ra=RA*u.degree, dec=Dec*u.degree, frame='icrs')
        
        # Use SDSS cone search
        result = SDSS.query_region(coord, radius=search_radius*u.degree, 
                                  spectro=True, data_release=16, photoobj_fields=[ 'objID', 'type', 'flags', 'clean', 'score', 'specObjID', 'u', 'err_u', 'g', 'err_g', 'r', 'err_r', 'i', 'err_i', 'z', 'err_z', 'ra', 'dec'])
        # result = SDSS.query_region(coord, radius=search_radius*u.degree, 
        #                          spectro=True, data_release=16)
        
        
        if result is not None and len(result) > 0:
            print(f"Found {len(result)} objects in search area around {name}")
            
            # Filter for galaxies (type = 3) with good quality
            if 'type' in result.colnames:
                galaxy_mask = result['type'] == 3
                result = result[galaxy_mask]
            
            # Additional quality filters if columns exist
            if 'clean' in result.colnames:
                clean_mask = result['clean'] == 1
                result = result[clean_mask]
            
            # if 'score' in result.colnames:
            #     score_mask = result['score'] == 1
            #     result = result[score_mask]
            
            if len(result) > 0:
                # Calculate angular distances
                result_coords = SkyCoord(ra=result['ra']*u.degree, dec=result['dec']*u.degree, frame='icrs')
                angular_distances = target_coord.separation(result_coords)
                
                # Add angular distance column in arcseconds
                result['angular_distance_arcsec'] = angular_distances.to(u.arcsec).value
                result['SNID'] = name
                #add in the ra and dec of the sn 
                result['SN_RA'] = RA
                result['SN_Dec'] = Dec
                result['Galaxy_RA_sdss'] = result['ra']
                result['Galaxy_Dec_sdss'] = result['dec']
                # Sort by angular distance
                result.sort('angular_distance_arcsec')
                
                # Save results
                #result.write(output_file, format='csv', overwrite=True)
                
                # Print summary
                print(f"\nSummary for {name}:")
                print(f"Target coordinates: RA={RA:.6f}°, Dec={Dec:.6f}°")
                print(f"Search radius: {search_radius}° ({search_radius*60:.1f} arcminutes)")
                print(f"Galaxies found: {len(result)}")
                if len(result) > 0:
                    print(f"Closest galaxy: {result['angular_distance_arcsec'][0]:.2f} arcseconds away")
                    print(f"Farthest galaxy: {result['angular_distance_arcsec'][-1]:.2f} arcseconds away")
                
                return result
            else:
                print(f"No galaxies found within {search_radius}° of {name} after filtering")
                return None
        else:
            print(f"No objects found in the search area around {name}")
            return None
            
    except Exception as e:
        print(f"Error during cone search: {e}")
        print("Trying alternative approach...")
        
        # Fallback to simple SQL query without complex joins
        try:
            sql_query = f"""
            SELECT TOP 1000 type, flags, objID, specObjID, u, err_u, g, err_g, r, err_r, i, err_i, z, err_z, ra, dec
            FROM PhotoObj
            WHERE type = 3 
              AND clean = 1 
              AND score = 1 
              AND ra BETWEEN {ra_min} AND {ra_max}
              AND dec BETWEEN {dec_min} AND {dec_max}
            """
            
            result = SDSS.query_sql(sql_query)
            
            if result is not None and len(result) > 0:
                print(f"Found {len(result)} galaxies in search area around {name}")
                
                # Calculate angular distances
                result_coords = SkyCoord(ra=result['ra']*u.degree, dec=result['dec']*u.degree, frame='icrs')
                angular_distances = target_coord.separation(result_coords)
                
                # Filter by search radius
                within_radius = angular_distances < search_radius * u.degree
                result = result[within_radius]
                
                if len(result) > 0:
                    # Add angular distance column in arcseconds
                    result['angular_distance_arcsec'] = angular_distances[within_radius].to(u.arcsec).value
                    
                    # Sort by angular distance
                    result.sort('angular_distance_arcsec')
                    
                    # Save results only if output_file is provided
                    # if output_file is not None:
                    #     result.write(output_file, format='csv', overwrite=True)
                    print(f"Found {len(result)} galaxies within {search_radius}° of {name}")
                    # if output_file is not None:
                    #     print(f"Saved {len(result)} galaxies to {output_file}")
                    
                    return result
                else:
                    print(f"No galaxies found within {search_radius}° of {name}")
                    return None
            else:
                print(f"No galaxies found in the search area around {name}")
                return None
                
        except Exception as e2:
            print(f"Error with fallback query: {e2}")
            return None
    
    #save the result to a csv file
    

def read_source_catalogs(table_sdss, table_initial, sn_name, catalog_type='B'):
    print(f"\n=== Processing {sn_name} ===")
    
    # Handle both Astropy Table and pandas DataFrame
    try:
        if hasattr(table_sdss, 'shape'):
            print(f"Input table_sdss shape: {table_sdss.shape}")
            print(f"Input table_sdss columns: {list(table_sdss.columns)}")
        else:
            print(f"Input table_sdss type: {type(table_sdss)}")
            print(f"Input table_sdss length: {len(table_sdss)}")
            if hasattr(table_sdss, 'colnames'):
                print(f"Input table_sdss columns: {list(table_sdss.colnames)}")
            else:
                print(f"Input table_sdss columns: {list(table_sdss.columns)}")
    except Exception as debug_e:
        print(f"Debug error: {debug_e}")
        print(f"Input table_sdss type: {type(table_sdss)}")
    
    if catalog_type == 'B':
        table_extractor = pd.read_csv('data_collection_files/data_files/brightest_galaxy.csv')    
    elif catalog_type == 'P':
        table_extractor = pd.read_csv('data_collection_files/data_files/probable_galaxy_catalog.csv')

    # pick a subset of the table where the SNID is the same as the sn_name
    table_subset = table_initial[table_initial['SNID_x'] == sn_name]
    print(f"Found {len(table_subset)} rows for {sn_name}")


    # Get unique filters for this SN from the table_subset (from the table_initial with all the images)
    filters = table_subset['Filter_x'].unique()
    #print(f"Processing {len(filters)} filters for {sn_name}: {filters}")
    
    # Process each filter once
    for f in filters:
        print(f"\nProcessing filter: {f}")
        # Get all rows for this specific filter
        table_subset_filter = table_subset[table_subset['Filter_x'] == f]
        filenames_in_f = table_subset_filter['filename']
        #print(f"Found {len(filenames_in_f)} files for filter {f}")
        
        # Extract just the filename from the full path
        filenames_only = filenames_in_f.str.split('/').str[-1].str.split('.').str[0]
        print(f"Filenames: {list(filenames_only)}")
            
        # Calculate zp_ab for each file in this filter
        zp_ab_values = []
        for filename in filenames_in_f:
            # Get the base path from this specific filename
            base_path = '/astro/armin/bhoomika/pantheon_data_folder'
            #print('base_path for file', filename, ':', base_path)
            
            # Extract FITS data for this specific file - returns zp_ab directly
            try:
                #open the fits file here directly and extract the zp_ab
                with fits.open(filename) as hdul:
                    header = hdul[1].header
                    
                    # Debug: Print all header keywords to see what's available
                    print(f"\n=== Header keywords for {filename} ===")
                    for key in sorted(header.keys()):
                        if any(x in key.upper() for x in ['PHOT', 'FILTER', 'INSTRUME', 'TELESCOP', 'ZPT', 'ZERO', 'MAG', 'FLUX']):
                            print(f"  {key}: {header[key]}")
                    
                    # For HST data, look for different keywords
                    filter_name = header.get("FILTER")
                    instrument = header.get("INSTRUME")
                    telescope = header.get("TELESCOP")
                    
                    # Look for zeropoint-related keywords
                    phot_zpt = header.get("PHOTZPT") or header.get("ZPT") or header.get("ZEROPOINT") or header.get("MAGZPT") 
                    
                    # Look for flux-related keywords
                    phot_flam = header.get("PHOTFLAM") or header.get("FLAM") or header.get("FLUXZERO")
                    phot_plam = header.get("PHOTPLAM") or header.get("PLAM") or header.get("PIVOT")
                    
                    print(f"\nExtracted values:")
                    print(f"  FILTER: {filter_name}")
                    print(f"  INSTRUMENT: {instrument}")
                    print(f"  TELESCOPE: {telescope}")
                    print(f"  ZEROPOINT (from header): {phot_zpt}")
                    print(f"  PHOTFLAM/FLAM: {phot_flam}")
                    print(f"  PHOTPLAM/PLAM: {phot_plam}")
                    
                    # For HST, use the direct zeropoint if available, or calculate from other keywords
                    if phot_zpt is not None and phot_zpt > 0:
                        zp_ab = phot_zpt
                        print(f"  Using direct zeropoint: {zp_ab}")
                    elif phot_flam is not None and phot_plam is not None and phot_flam > 0 and phot_plam > 0:
                        zp_ab = -2.5 * np.log10(phot_flam) - 5 * np.log10(phot_plam) - 2.408
                        print(f"  Calculated zp_ab from PHOTFLAM/PHOTPLAM: {zp_ab}")
                    else:
                        # For HST WFC3, use standard zeropoints based on filter
                        hst_zeropoints = {
                            'F160W': 25.96,  # WFC3/IR F160W
                            'F105W': 26.27,  # WFC3/IR F105W  
                            'F110W': 26.82,  # WFC3/IR F110W
                            'F125W': 26.23,  # WFC3/IR F125W
                            'F140W': 26.45,  # WFC3/IR F140W
                            'F435W': 25.68,  # WFC3/UVIS F435W
                            'F475W': 26.16,  # WFC3/UVIS F475W
                            'F555W': 25.72,  # WFC3/UVIS F555W
                            'F606W': 26.40,  # WFC3/UVIS F606W
                            'F625W': 25.69,  # WFC3/UVIS F625W
                            'F775W': 25.25,  # WFC3/UVIS F775W
                            'F814W': 25.96,  # WFC3/UVIS F814W
                            'F850LP': 24.87  # WFC3/UVIS F850LP
                        }
                        
                        if filter_name in hst_zeropoints:
                            zp_ab = hst_zeropoints[filter_name]
                            print(f"  Using standard HST zeropoint for {filter_name}: {zp_ab}")
                            print(f"  WARNING: This is a standard value - may not be correct for your data!")
                        else:
                            print(f"  Warning: No zeropoint found for filter {filter_name}")
                            zp_ab = 0
                    
                    zp_ab_values.append(zp_ab) #zp values for each file in this filter
                #print('zp_ab for file', filename, ':', zp_ab)
            except Exception as e:
                print(f"Error extracting FITS data for {filename}: {e}")
                zp_ab_values.append(0)
        
        # Calculate average zp_ab for this filter
        avg_zp_ab = np.mean(zp_ab_values) if zp_ab_values else 0
        print('Average zp_ab for filter', f, ':', avg_zp_ab)
        
        # Add columns to table_sdss for this filter (INSIDE the filter loop)
        matching_rows = table_extractor[table_extractor['File_key'].isin(filenames_only)].copy()
        if len(matching_rows) > 0:
            # Convert string columns to numeric
            matching_rows['Mag'] = pd.to_numeric(matching_rows['Mag'], errors='coerce')
            matching_rows['Mag_Error'] = pd.to_numeric(matching_rows['Mag_Error'], errors='coerce')
            matching_rows['Flux'] = pd.to_numeric(matching_rows['Flux'], errors='coerce')
            matching_rows['Flux_Error'] = pd.to_numeric(matching_rows['Flux_Error'], errors='coerce')
            
            if len(table_subset_filter) > 1:
                # For multiple files, average the values
                instrumental_mag = matching_rows['Mag'].mean()
                table_sdss[f + '_Mag'] = instrumental_mag + avg_zp_ab
                print(f"instrumental magnitude before correction: {instrumental_mag}")
                table_sdss[f + '_Mag_Error'] = matching_rows['Mag_Error'].mean() 
                table_sdss[f + '_Flux'] = matching_rows['Flux'].mean()
                table_sdss[f + '_Flux_Error'] = matching_rows['Flux_Error'].mean()
                ra_sex, dec_sex = matching_rows['RA_Galaxy'].mean(), matching_rows['Dec_Galaxy'].mean()
                ra_sdss, dec_sdss = matching_rows['Galaxy_RA_sdss'].mean(), matching_rows['Galaxy_Dec_sdss'].mean()
                #get the distance between the sex and sdss
                distance = angular_distance(ra_sex, dec_sex, ra_sdss, dec_sdss)
                table_sdss['Galaxy_Distance_sdds_sex'] = distance
            else:
                # For single file, get the first matching value
                instrumental_mag = matching_rows['Mag'].iloc[0]
                table_sdss[f + '_Mag'] = instrumental_mag + avg_zp_ab
                print(f"instrumental magnitude before correction: {instrumental_mag}")
                table_sdss[f + '_Mag_Error'] = matching_rows['Mag_Error'].iloc[0]
                table_sdss[f + '_Flux'] = matching_rows['Flux'].iloc[0]
                table_sdss[f + '_Flux_Error'] = matching_rows['Flux_Error'].iloc[0]
                ra_sex, dec_sex = matching_rows['RA_Galaxy'].iloc[0], matching_rows['Dec_Galaxy'].iloc[0]
                ra_sdss, dec_sdss = matching_rows['Galaxy_RA_sdss'].iloc[0], matching_rows['Galaxy_Dec_sdss'].iloc[0]
                #get the distance between the sex and sdss
                distance = angular_distance(ra_sex, dec_sex, ra_sdss, dec_sdss)
                table_sdss['Galaxy_Distance_sdds_sex'] = distance
        else:
            print(f"No matching rows found for filter {f} with filenames: {filenames_only}")
        
        print(f"Added columns for filter {f}: {f}_Mag, {f}_Mag_Error, {f}_Flux, {f}_Flux_Error")

    # Reorganize the table: put SNID first and rename columns to clarify these are galaxies
    try:
        # Convert to pandas DataFrame for easier manipulation
        if hasattr(table_sdss, 'to_pandas'):
            df = table_sdss.to_pandas()
        else:
            df = table_sdss.copy()
        
        # Rename columns to clarify these are galaxies
        column_mapping = {
            'objID': 'Galaxy_objID',
            'ra': 'Galaxy_RA',
            'dec': 'Galaxy_Dec',
            'u': 'Galaxy_u_mag',
            'g': 'Galaxy_g_mag', 
            'r': 'Galaxy_r_mag',
            'i': 'Galaxy_i_mag',
            'z': 'Galaxy_z_mag',
            'err_u': 'Galaxy_u_err',
            'err_g': 'Galaxy_g_err',
            'err_r': 'Galaxy_r_err', 
            'err_i': 'Galaxy_i_err',
            'err_z': 'Galaxy_z_err',
            'angular_distance_arcsec': 'Galaxy_distance_arcsec',
            'type': 'Galaxy_type',
            'flags': 'Galaxy_flags',
            'clean': 'Galaxy_clean',
            'score': 'Galaxy_score',
            'specObjID': 'Galaxy_specObjID'
        }
        
        # Apply the column renaming
        df = df.rename(columns=column_mapping)
        
        # Reorder columns to put SNID first, then galaxy info, then HST data
        base_columns = ['SNID']
        galaxy_columns = [col for col in df.columns if col.startswith('Galaxy_')]
        hst_columns = [col for col in df.columns if any(filter_name in col for filter_name in ['F160W', 'F125W', 'F110W', 'F105W', 'F140W', 'F435W', 'F475W', 'F555W', 'F606W', 'F625W', 'F775W', 'F814W', 'F850LP', 'F275W', 'F336W'])]
        other_columns = [col for col in df.columns if col not in base_columns + galaxy_columns + hst_columns]
        
        # Create the new column order
        new_column_order = base_columns + galaxy_columns + hst_columns + other_columns
        df = df[new_column_order]
        
        # Convert back to Astropy Table
        from astropy.table import Table
        table_sdss = Table.from_pandas(df)
        
        print(f"Final table_sdss type: {type(table_sdss)}")
        print(f"Final table_sdss length: {len(table_sdss)}")
        print(f"Final table_sdss columns: {list(table_sdss.colnames)}")
        print(f"Final table_sdss:\n{table_sdss}")
        
    except Exception as debug_e:
        print(f"Final debug error: {debug_e}")
        print(f"Final table_sdss type: {type(table_sdss)}")
        
    return table_sdss

    



def read_ir_data(table):
    """
    Query SDSS for galaxies near each supernova in the table
    """
    all_results = []
    output_file = 'sdss_crossmatched_galaxies_data.csv'
    processed_sns = set()  # Track which supernovae we've already processed
    
    # Print column names to debug
    #print("Available columns:", table.columns.tolist())
    
    for i in range(len(table)):
        # Try different possible column names for SNID
        if 'SNID_x' in table.columns:
            name = table['SNID_x'][i]
        elif 'SNID' in table.columns:
            name = table['SNID'][i]
        else:
            print(f"Warning: No SNID column found. Available columns: {table.columns.tolist()}")
            continue
            
        # Skip if we've already processed this supernova
        if name in processed_sns:
            print(f"Skipping {name} - already processed")
            continue
            
        ra = table['RA_SN'][i]
        dec = table['Dec_SN'][i]
        search_radius = 0.05  # 0.05 degrees = 3 arcminutes (SDSS limit)
        
        print(f"\n--- Processing {i+1}/{len(table)}: {name} ---")
        
        # Mark this supernova as processed    
        try:
            # Don't pass output_file to avoid overwriting
            result = get_sdss_good_coordinates_all(name, ra, dec, search_radius)
            if result is not None:
                print(f"Successfully found {len(result)} galaxies near {name}")
                # Add SNID column to identify which supernova this belongs to
                result['SNID'] = name
                
                # Check which type of source extraction we need to use
                try:
                    if 'B_P_y' in table.columns:
                        if table['B_P_y'][i] == 'B':
                            type = 'B'
                        elif table['B_P_y'][i] == 'P':
                            type = 'P'
                        elif table['B_P_y'][i] == 'BP' or table['B_P_y'][i] == 'PB':
                            type = 'B'  # let use the brightest galaxy for now
                        else:
                            type = 'B'  # default
                    else:
                        type = 'B'  # default if column doesn't exist
                    
                    result = read_source_catalogs(result, table, name, catalog_type=type)
                    processed_sns.add(name)
                    all_results.append(result)
                except Exception as e:
                    print(f"Error in read_source_catalogs for {name}: {e}")
                    # Still add the result even if source catalog processing fails
                    all_results.append(result)
            else:
                print(f"No galaxies found near {name}")
        except Exception as e:
            print(f"Error processing {name}: {e}")
            continue
        
        # Add a small delay to avoid overwhelming the SDSS server
        time.sleep(1)
    
    # Combine all results and save once at the end
    if all_results:
        from astropy.table import vstack
        combined_results = vstack(all_results)
        #create the sdss file 
        combined_results.write('SDSS_data/'+output_file, format='csv', overwrite=True)
        print(f"\nSaved {len(combined_results)} total galaxies to {output_file}")
        return combined_results
    else:
        print("No galaxies found for any supernovae")
        return None


def plot_sdss_and_hst(file_name):
    """
    Plot SDSS and HST data for galaxies, only showing filters that have data
    """
    import matplotlib.pyplot as plt
    import numpy as np
    
    # Read the file
    sdss_data = pd.read_csv(file_name)
    print(f"Loaded data with {len(sdss_data['SNID'].unique())} galaxies")
    print(f"Loaded data with {sdss_data['SNID'].unique()}")
   
    #print(f"Available columns: {list(sdss_data.columns)}")
    
    # Define filter wavelengths (in Angstroms)
    filter_wavelengths = {
        # SDSS filters
        'Galaxy_u_mag': 3551,   # SDSS u-band
        'Galaxy_g_mag': 4686,   # SDSS g-band  
        'Galaxy_r_mag': 6166,   # SDSS r-band
        'Galaxy_i_mag': 7480,   # SDSS i-band
        'Galaxy_z_mag': 8932,   # SDSS z-band
        
        # HST WFC3/IR filters
        'F105W_Mag': 10552,         # WFC3/IR F105W
        'F110W_Mag': 11534,         # WFC3/IR F110W
        'F125W_Mag': 12486,         # WFC3/IR F125W
        'F140W_Mag': 13923,         # WFC3/IR F140W
        'F160W_Mag': 15369,         # WFC3/IR F160W
        
        # HST WFC3/UVIS filters
        'F275W_Mag': 2704,          # WFC3/UVIS F275W
        'F336W_Mag': 3354,          # WFC3/UVIS F336W
        'F435W_Mag': 4328,          # WFC3/UVIS F435W
        'F475W_Mag': 4747,          # WFC3/UVIS F475W
        'F555W_Mag': 5308,          # WFC3/UVIS F555W
        'F606W_Mag': 5921,          # WFC3/UVIS F606W
        'F625W_Mag': 6311,          # WFC3/UVIS F625W
        'F775W_Mag': 7693,          # WFC3/UVIS F775W
        'F814W_Mag': 8024,          # WFC3/UVIS F814W
        'F850LP_Mag': 9055,         # WFC3/UVIS F850LP
    }
    
    # Find which filters have data (non-null, non-zero values)
    available_filters = {}
    for filter_name, wavelength in filter_wavelengths.items():
        if filter_name in sdss_data.columns:
            # Check if there's any valid data (not NaN, not 0)
            valid_data = sdss_data[filter_name].dropna()
            valid_data = valid_data[valid_data != 0]  # Remove zeros
            
            if len(valid_data) > 0:
                available_filters[filter_name] = wavelength
                #print(f"Found data for {filter_name}: {len(valid_data)} valid measurements")
            else:
                print(f"No valid data for {filter_name}")
        else:
            print(f"Column {filter_name} not found in data")
    
    if not available_filters:
        print("No valid filter data found!")
        return
    
    # Create the plot
    plt.figure(figsize=(12, 8))
    
    # Plot each galaxy as a separate line
    for idx, galaxy in sdss_data.iterrows():
        x_values = []
        y_values = []
        filters = []
        total_galaxies = []
        
        for filter_name, wavelength in available_filters.items():
            if pd.notna(galaxy[filter_name]) and galaxy[filter_name] != 0:
                x_values.append(wavelength)
                y_values.append(galaxy[filter_name])
                filters.append(filter_name)
        
        # Check if there's any IR data (wavelength >= 10000) before plotting
        has_ir_data = any(wavelength >= 10000 for wavelength in x_values)
        
        if len(x_values) > 0 and has_ir_data:  # Only plot if there's data AND IR data
            plt.figure(figsize=(10, 6))
            plt.title(f'Galaxy :{galaxy['objID']}, SN :{galaxy['SNID']}')
            plt.plot(x_values, y_values, 'o-', alpha=0.7, markersize=4, marker='o', ls=None, label = filters)
            plt.gca().invert_yaxis()
            plt.xlabel('Wavelength (Angstroms)')
            plt.ylabel('Magnitude')
            plt.grid(True, alpha=0.3)
            plt.legend()
            plt.tight_layout()
            plt.show(block=False)
            total_galaxies.append(galaxy['SNID'])
            
            # Wait for user input to continue
            input(f"Press Enter to continue to next galaxy (current: {galaxy['objID']})...")
            plt.close()
            
        elif len(x_values) > 0:
            print(f"Skipping galaxy {galaxy['Galaxy_objID']} - no IR data (only wavelengths < 10000)")
    
    print("All plots completed!")
    print(f"Total galaxies: {len(total_galaxies.unique())}")
    print(f"Total galaxies: {total_galaxies.unique()}")
    



if __name__ == "__main__":

    if not os.path.exists('SDSS_data/sdss_crossmatched_galaxies_data.csv'):
        plot_sdss_and_hst('SDSS_data/sdss_crossmatched_galaxies_data.csv')
    else:
        print("File does not exist, creating it")
        #read_ir_data(filtered_table)

        # Run example when script is executed directly
        table1 = pd.read_csv('data_collection_files/data_files/reduced_table_sn_in_images.csv')
        #reneme the SNID_x column to SNID
        table1 = table1.rename(columns={'SNID_x': 'SNID'})
        table2 = pd.read_csv('data_collection_files/data_files/fits_summary_combined.csv')
        
        # Use inner join to only get SNIDs that exist in both tables
        # This prevents the Cartesian product issue
        table = pd.merge(table1, table2, on='filename', how='inner')
        
        # Filter for 'Good_Bad_y' starting with 'T', handling NaN values
        filtered_table = table[table['Good_Bad'].fillna('').str.startswith('T')]

        #check how many galaxies are matched with the brightest galaxy
        print(f"Filtered table length: {len(filtered_table)}")
        print(f"Number of galaxies matched with the brightest galaxy: {len(filtered_table[filtered_table['B_P'] == 'B'])}")
        print(f"The galaxies that are matched with the brightest galaxy: {filtered_table[filtered_table['B_P'] == 'B']['SNID_x'].unique()}")
     
        
        # If you want to keep all entries from table1 but only matching ones from table2:
        # table = pd.merge(table1, table2, on='SNID', how='left')
        filtered_table.reset_index(drop=True, inplace=True)
        read_ir_data(filtered_table)
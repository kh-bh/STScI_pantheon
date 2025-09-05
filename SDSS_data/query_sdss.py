'''
This scrops queries the optical observations of the supernovae from the SDSS database and merges
the data with the reduced table of the supernovae that are in the images

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

# Import the function
from magnitude_corrector_calculator import extract_fits_by_snid

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
    

def read_source_catalogs(table_sdss, table_initial, sn_name, type ='B'):
    if type == 'B':
        table_extractor = pd.read_csv('data_collection_files/data_files/brightest_galaxy.csv')    
    elif type == 'P':
        table_extractor = pd.read_csv('data_collection_files/data_files/probable_galaxy_catalog.csv')

    # pick a subset of the table where the SNID is the same as the sn_name
    table_subset = table_initial[table_initial['SNID_x'] == sn_name]

    for i in range(len(table_subset)):
        filter = table_subset['Filter_x'].unique()
         #since it starts at 1
        for f in filter:
            #if table_subset[filter] is greater than 1, then we need to average the mag and error
            table_subset_filter = table_subset[table_subset['Filter_x'] == f]
            filenames_in_f = table_subset_filter['filename']
            # Extract just the filename from the full path
            filenames_only = filenames_in_f.str.split('/').str[-1].str.split('.').str[0]
            print('filenames_only', filenames_only)
            
            # Calculate zp_ab for each file in this filter
            zp_ab_values = []
            for filename in filenames_in_f:
                # Get the base path from this specific filename
                base_path = os.path.dirname(filename)
                print('base_path for file', filename, ':', base_path)
                
                # Extract FITS data for this specific file - returns zp_ab directly
                try:
                    #open the fits file here directly and extract the zp_ab
                    with fits.open(filename) as hdul:
                        header = hdul[0].header
                        phot_mode = header.get("PHOTMODE")
                        phot_flam = header.get("PHOTFLAM")
                        phot_zpt = header.get("PHOTZPT")
                        phot_plam = header.get("PHOTPLAM")
                    zp_ab = extract_fits_by_snid([sn_name], base_path=base_path, output_csv=None, recursive=True, data_list=[])
                    zp_ab_values.append(zp_ab)
                    print('zp_ab for file', filename, ':', zp_ab)
                except Exception as e:
                    print(f"Error extracting FITS data for {filename}: {e}")
                    zp_ab_values.append(0)
            
            # Calculate average zp_ab for this filter
            avg_zp_ab = np.mean(zp_ab_values) if zp_ab_values else 0
            print('Average zp_ab for filter', f, ':', avg_zp_ab)
            
            if len(table_subset_filter) > 1:
                # For multiple files, average the values
                table_sdss[f + '_Mag'] = (table_extractor[table_extractor['File_key'].isin(filenames_only)]['Mag'] + zp_ab_values).mean()
                table_sdss[f + '_Mag_Error'] = table_extractor[table_extractor['File_key'].isin(filenames_only)]['Mag_Error'].mean() 
                table_sdss[f + '_Flux'] = table_extractor[table_extractor['File_key'].isin(filenames_only)]['Flux'].mean()
                table_sdss[f + '_Flux_Error'] = table_extractor[table_extractor['File_key'].isin(filenames_only)]['Flux_Error'].mean()
            else:
                # For single file, get the first matching value
                matching_rows = table_extractor[table_extractor['File_key'].isin(filenames_only)]
                if len(matching_rows) > 0:
                    table_sdss[f + '_Mag'] = matching_rows['Mag'].iloc[0] + zp_ab_values
                    table_sdss[f + '_Mag_Error'] = matching_rows['Mag_Error'].iloc[0]
                    table_sdss[f + '_Flux'] = matching_rows['Flux'].iloc[0]
                    table_sdss[f + '_Flux_Error'] = matching_rows['Flux_Error'].iloc[0]

        
    return table

    



def read_ir_data(table):
    """
    Query SDSS for galaxies near each supernova in the table
    """
    all_results = []
    output_file = 'sdss_crossmatched_galaxies_data.csv'
    
    # Print column names to debug
    print("Available columns:", table.columns.tolist())
    
    for i in range(len(table)):
        # Try different possible column names for SNID
        if 'SNID_x' in table.columns:
            name = table['SNID_x'][i]
        elif 'SNID' in table.columns:
            name = table['SNID'][i]
        else:
            print(f"Warning: No SNID column found. Available columns: {table.columns.tolist()}")
            continue
            
        ra = table['RA_SN'][i]
        dec = table['Dec_SN'][i]
        search_radius = 0.05  # 0.05 degrees = 3 arcminutes (SDSS limit)
        
        print(f"\n--- Processing {i+1}/{len(table)}: {name} ---")
        
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
                    
                    result = read_source_catalogs(result, table, name, type=type)
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
        combined_results.write('SDSS_DATA/'+output_file, format='csv', overwrite=True)
        print(f"\nSaved {len(combined_results)} total galaxies to {output_file}")
        return combined_results
    else:
        print("No galaxies found for any supernovae")
        return None

    
if __name__ == "__main__":
    # Run example when script is executed directly
    table1 = pd.read_csv('data_collection_files/data_files/reduced_table_sn_in_images.csv')
    #reneme the SNID_x column to SNID
    table1 = table1.rename(columns={'SNID_x': 'SNID'})
    table2 = pd.read_csv('data_collection_files/data_files/fits_summary_combined.csv')
    
    # Use inner join to only get SNIDs that exist in both tables
    # This prevents the Cartesian product issue
    table = pd.merge(table1, table2, on='filename', how='inner')
    
    # Filter for 'Good_Bad_y' starting with 'T', handling NaN values
    filtered_table = table[table['Good_Bad_y'].fillna('').str.startswith('T')]
    print(f"Merged table length: {len(table)}")
    print(f"Filtered table length: {len(filtered_table)}")
    
    # If you want to keep all entries from table1 but only matching ones from table2:
    # table = pd.merge(table1, table2, on='SNID', how='left')
    filtered_table.reset_index(drop=True, inplace=True)
    read_ir_data(filtered_table[0:20])
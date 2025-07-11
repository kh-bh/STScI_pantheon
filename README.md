# Guidelines to Files and Order to be Run in

Information for pantheon files

## Order to Run

Step 1. Run data_collection_Hubble.py and data_collection_JWST

Step 2. Run download_data.ipynb or download_data.py (make sure to manually change from HST to JWST)

Step 3. Run fits_data_table.py

Step 4. Run catalog_extractor.py

Step 5. Run extractor_brightest_info_data_table.py and extractor_radius_info_data_table.py

Step 6. Run magnitude_corrector_calculator.py

## File Info

### Folder: data_collection_files
- Where all the files related to data collection are

#### data_collection_HST.py: 
If there are images for the supernovae in HST at the resolved query, 
create a file that stores the point of the list in the order "SNID, RA, DEC"
    - File name that is returned back is valid_points_index_list_HST.csv

#### data_collection_JWST.py: 
If there are images for the supernovae in HST at the resolved query, 
create a file that stores the point of the list in the order "SNID, RA, DEC"
    - File name that is returned back is valid_points_index_list_JWST.csv

#### download_data.ipynb:
Downloads the data from the obeservations of either the JWST or HST valid points 
  - **IMPORTANT** must change from HST or JWST in corrected

#### download_data.py:
Same as download_data.ipynb but in a python file

#### extractor_brightest_info_data_table.py:
Create a datatable with the BRIGHTEST and LARGEST objects of each source extractor file
    - File returned back is brightest_galaxy.csv

#### extractor_radius_info_data_table.py:
Create a datatable with the RADIUS of objects of around the RA and DEC of each SN1 source extractor file
    - File returned back is radius_galaxy.csv

#### fits_data_table.py:
Create a datatable with most important information from the fits file 
    - File returned back is fits_summary.csv

#### magnitude_corrector_calculator.py:
Create a datatable that corrects for each magnitude of the brightest galaxy or the radius galaxy
    - File returned is corrected_mag.csv

### Folder: data_files  
- Where all of the data files end up going

#### brightest_galaxy.csv:
- SNID,File_key,Galaxy_Key,RA_Galaxy,Dec_Galaxy,Position_Angle,Major_Axis,Minor_Axis,Mag,Mag_Error,Flux,Flux_Error
- A list of the brightest galaxies from each of the MAST queries

#### fits_summary.csv:
- SNID,filename,Date,Telescope,Instrument,Filter,CD1_1,CD1_2,CD2_1,CD2_2
- Records the fits information of each image you have

#### HST.csv:
- RA,Dec,SNID,has_f1
- Este's HST image catalog from May 

#### magnitude_calculated.csv:
- SNID,filename,Date,Photometric Calibration,Inverse Sensitivity,Zeropoint,Wavelength,Zeropoint AB
- A list of the folders and their correcected magnitudes

#### unmatched_pantheon.csv
- SNID,IAUC,host,RA,Dec,RA_host,Dec_host,zhel,zcmb,zhelerr,zHD,zHDerr,PV,vpecerr,RA_group,Dec_group,zhel_group,zcmb_group,zHD_group,PV_group,in_group,has_host,is_SNz,resolved_coord
- The orginal file list AND resolved_coord
  - resolved_coord is a string of 'RA Dec' form


#### unmatched_pantheon(in).csv:
- SNID,IAUC,host,RA,Dec,RA_host,Dec_host,zhel,zcmb,zhelerr,zHD,zHDerr,PV,vpecerr,RA_group,Dec_group,zhel_group,zcmb_group,zHD_group,PV_group,in_group,has_host,is_SNz
- The true orginal file

#### valid_points_index_list_HST.csv:
- SNID,RA,Dec
- Supernovae with the points that have IR or NIR data for HST (as of June 2025)

#### valid_points_index_list_JWST.csv:
- SNID,RA,Dec
- Supernovae with the points that have *data* for JWST (as of June 2025)

### Folder: pantheon_data_folder
- Each supernova has a folder with all of the MAST downloads
- Ideally, it will be split by HST, HLA, and JWST
- To get this folder, run download_data
- **IMPORTANT** this folder takes up a lot of space 

### Folder: source_extractor
- For each supernova image in the pantheon_data_folder, this holds all of the source extractor files

### Main folder:

#### default.cov:
- The configuration file for source extractor 

#### default.param:
- the parameter file for source extractor

#### default.sex:
- the source extractor file

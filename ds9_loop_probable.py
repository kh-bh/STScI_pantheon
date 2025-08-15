#!/usr/bin/env python

import sys, os, glob, re, subprocess, shutil
import argparse
import pandas as pd
import numpy as np
import math 
import sep
from astropy.io import fits
from pathlib import Path 




def define_args():
    parser = argparse.ArgumentParser(description="Generate DS9 region files", conflict_handler='resolve')
    
    if 'MASS_STEP_OUT_ROOTDIR' in os.environ:
        outrootdir = os.path.abspath(f'{os.environ["MASS_STEP_OUT_ROOTDIR"]}')
        print(outrootdir)
    else:
        outrootdir = '.'

    if 'MASS_STEP_DATA_ROOTDIR' in os.environ:
        datadir = os.path.abspath(f'{os.environ["MASS_STEP_DATA_ROOTDIR"]}')
        print(datadir)
    else:
        datadir = '.'
    if 'MASS_STEP_TABLES_ROOTDIR' in os.environ:
        tabledir = os.path.abspath(f'{os.environ["MASS_STEP_TABLES_ROOTDIR"]}')
        print(tabledir)
    else:
        tabledir = '.'

    parser.add_argument(
        "supernova_names", 
        type=str, 
        nargs="+", 
        help="List of supernova names to process."
    )

    print('GGGGGG',tabledir)
    parser.add_argument(
        "--fits_summary_path", 
        type=str, 
        default=f'{tabledir}/data_collection_files/data_files/fits_cleaned_HST_data.csv', 
        help="Path to the FITS summary CSV file (default=%(default)s)"
    )

    parser.add_argument(
        "--unmatched_pantheon_path", 
        type=str, 
        default=f'{tabledir}/data_collection_files/data_files/unmatched_pantheon.csv', 
        help="Path to the list of supernovas CSV file (default=%(default)s)"
    )

    parser.add_argument(
        "--probable_galaxy_path", 
        type=str, 
        default=f'{tabledir}/data_collection_files/data_files/probable_galaxy_catalog.csv', 
        help="Path to the probable galaxy CSV file (default=%(default)s)"
    )

    parser.add_argument(
        "--all_output_filename", 
        type=str, 
        default="ds9_cmds_all.reg", 
        help="Filename for the combined region file output (default=%(default)s)"
    )

    parser.add_argument(
        "--out_dir", 
        type=str, 
        default=f'{outrootdir}', 
        help="Directory where individual region files will be saved (default=%(default)s)"
    )
    parser.add_argument(
        "--datadir", 
        type=str, 
        default=f'{datadir}', 
        help="Directory where the data is (default=%(default)s)"
    )

    parser.add_argument(
        "--download",
        action="store_true",
        help="Dowloads and fpacks images into local machine (default=%(default)s)"
    )

    return parser.parse_args()



# ds9:
def make_ellipse_command_for_position(ra,dec,position_angle,major_axis,minor_axis,name=None,color='red',width=4):
    print(f'\nMaking ellipse command: RA_Galaxy {ra}, Dec_Galaxy {dec}, name {name}, Position_Angle {position_angle}, Major_Axis {major_axis}, Minor_Axis {minor_axis}')
    s = f'ellipse({ra},{dec},{minor_axis}",{major_axis}",{position_angle}) # color={color} width={width}'
    if name is not None:
        name='{'+f'{name}'+'}'
        s+=f' text={name}'

    #s = f'ellipse({ra},{dec},{major_axis}",{minor_axis}",{position_angle}) # color={color} width={width} text={name}'
    print(s)
    return(s)

def make_circle_command_for_position(ra,dec,radius_arcsec=1,name=None,color='blue',width=5):
    print(f'\nMaking circle command: RA {ra}, Dec {dec}, name {name},')
    name='{'+f'{name}'+'}'
    s = f'circle({ra},{dec},{radius_arcsec}") # color={color} width={width} text={name}'
    print(s)
    return(s)

# def get_region_info(cd_matrix, cov_matrix):
#     # 3. World-space covariance matrix (deg²)
#     cov_world = cd_matrix @ cov_matrix @ cd_matrix.T

#     # 4. Eigen-decomposition
#     eigvals, eigvecs = np.linalg.eigh(cov_world)  # always returns sorted values
    
#     a_arcsec = 3600 * np.sqrt(eigvals[1])  # major axis
#     b_arcsec = 3600 * np.sqrt(eigvals[0])  # minor axis

#     # 5. Position angle (in degrees, from North to East)
#     vec = eigvecs[:, 1]  # eigenvector of major axis
#     theta_rad = np.arctan2(vec[0], vec[1])  # careful: note the order (RA, Dec)
#     theta_deg = np.degrees(theta_rad)

#     theta_deg = (theta_deg + 360) % 180

#     return theta_deg, a_arcsec, b_arcsec


def get_region_info2(table):
    for ix, row in table.iterrows():
        cxx = row['CXX']
        cyy = row['CYY']
        cxy = row['CXY']

        (a, b, theta_rad) = sep.ellipse_axes(cxx, cyy, cxy)
        table.loc[ix,'a_pixels'] = a
        table.loc[ix,'b_pixels'] = b

      

        #theta_rad = 0.5 * math.atan2(2 * cxy, cxx - cyy)
        theta_deg = math.degrees(theta_rad)
        alpha_deg = np.degrees(np.arctan2(-row['CD1_2'], row['CD2_2']))
        table.loc[ix,'theta_deg'] = theta_deg
        table.loc[ix,'alpha_deg'] = alpha_deg
        table.loc[ix,'PA_deg'] = theta_deg-90.0-alpha_deg
        table.loc[ix,'theta_rad'] = theta_rad

        #print (f'theta_deg {theta_deg} {a} {b} {theta_rad}')
        #print (f'alpha_deg {alpha_deg}')

        # Pixel scale (arcsec/pixel)
        pixscale = 3600.0 * np.sqrt(math.fabs((row['CD1_1']*row['CD2_2']) - row['CD1_2']*row['CD2_1']))
        #pixscale = 3600.0 * np.sqrt(CD1_1*CD2_2 - CD1_2*CD2_1)
        table.loc[ix,'pixscale'] = pixscale

        #sys.exit(0)

        #a_arcsec = a*pixscale
        #b_arcsec = b*pixscale
        a_arcsec = a * pixscale * row["Kron_Radius"]
        b_arcsec = b * pixscale * row["Kron_Radius"]
 


        # a_arcsec = math.fabs(row['CD1_1'])*3600.0*a
        # b_arcsec = math.fabs(row['CD1_1'])*3600.0*b
        table.loc[ix,'a_arcsec'] = a_arcsec
        table.loc[ix,'b_arcsec'] = b_arcsec


        #continue

        # trace = cxx + cyy
        # delta = (cxx - cyy) / 2
        # root_term = math.sqrt(delta**2 + cxy**2)
        # lambda1 = (trace / 2) + root_term
        # lambda2 = (trace / 2) - root_term

     
        # a = math.sqrt(abs(max(lambda1, lambda2)))
        # b = math.sqrt(abs(min(lambda1, lambda2)))
        # table.loc[ix,'a_pixels'] = a
        # table.loc[ix,'b_pixels'] = b

        # theta_rad = 0.5 * math.atan2(2 * cxy, cxx - cyy)
        # theta_deg = math.degrees(theta_rad)
        # table.loc[ix,'theta_deg'] = theta_deg
        # table.loc[ix,'theta_rad'] = theta_rad

        # a_arcsec = math.fabs(row['CD1_1'])*a
        # b_arcsec = math.fabs(row['CD1_1'])*b
        # table.loc[ix,'a_arcsec'] = a
        # table.loc[ix,'b_arcsec'] = b

        #print(f'PA(deg) and a, b of ellipse:{PA_deg}, {a}, {b}')
    return(table)


def save_region_file(filepath, output):
    print(f'Writing region file to {filepath}')
    with open(filepath, "w") as file:
        for line in output:
            #print(line, 'printing line heeeeerrrrreeeeeee')
            file.write(str(line) + "\n")
    file.close()

#combining fits_summary and probable_galaxy
def combine_fs_and_bg(supernova_names, fits_summary, probable_galaxy):
    for sn_name in supernova_names:
        print("SN name: ",sn_name)

        sn_ix_bg = probable_galaxy[probable_galaxy["SNID"] == sn_name].index.values
        print("probable galaxy indices: ",sn_ix_bg)
        if len(sn_ix_bg) == 0:
            ######
            print(f'WARNING! could not find indices for galaxy {sn_name}, skipping')
            continue
#            raise RuntimeError(f'could not find indices for galaxy {sn_name}')

        sn_ix = np.where(fits_summary["SNID"] == sn_name)[0]
        print("Fits summary indices: ",sn_ix)

        for index in sn_ix_bg:
            filename = probable_galaxy.loc[index,'File_key']
            filename_ix_fs = sn_ix[np.where(fits_summary.loc[sn_ix, "filenameshort"] == filename+".fits")[0]]
            print(f'{filename} {filename_ix_fs}')
            if len(filename_ix_fs)==0:
                continue
            if len(filename_ix_fs)>1:
                raise RuntimeError(f"filename_ix_fs is returning multiple matches for filenameshort in fits_summary: {filename_ix_fs}")

            for column_name in probable_galaxy.columns:
                if column_name not in fits_summary.columns:
                    fits_summary[column_name] = np.nan
                fits_summary.at[filename_ix_fs[0], column_name] = probable_galaxy.at[index, column_name]
    return fits_summary



def combine_fs_and_up(supernova_names, fits_summary, unmatched_pantheon):
    for sn_name in supernova_names:
        print("SN name: ",sn_name)

        sn_ix_up = unmatched_pantheon[unmatched_pantheon["SNID"] == sn_name].index.values
        print(f'TTTTTT',sn_ix_up)
        if len(sn_ix_up)==0:
            raise RuntimeError(f'Could not find SN {sn_name}')
        elif len(sn_ix_up)>1:
            raise RuntimeError(f'more than one entry for SN {sn_name}')
        print('TEST1',unmatched_pantheon.loc[sn_ix_up[0],'RA'],unmatched_pantheon.loc[sn_ix_up[0],'Dec'])

        sn_ix_fs = fits_summary[fits_summary["SNID"] == sn_name].index.values
        if len(sn_ix_fs)==0:
            raise RuntimeError(f'Could not find entries in fits_summary for SN {sn_name}')

        for ix in sn_ix_fs:
            fits_summary.loc[ix,'RA_SN']=unmatched_pantheon.loc[sn_ix_up[0],'RA']
            fits_summary.loc[ix,'Dec_SN']=unmatched_pantheon.loc[sn_ix_up[0],'Dec']
            

    return fits_summary





def get_region_commands(index, fits_summary, ellipse_color='red', ellipse_width=1, sn_name=None):

    file_key = fits_summary.at[index, "File_key"]
    if isinstance(file_key, float) and np.isnan(file_key):
        return None, None, None

    ra_galaxy = fits_summary.at[index, "RA_Galaxy"]
    dec_galaxy = fits_summary.at[index, "Dec_Galaxy"]
   

    PA_deg, a_arcsec, b_arcsec = fits_summary.loc[index,['PA_deg','a_arcsec','b_arcsec']]

    cmds = []
    cmds.append(make_ellipse_command_for_position(ra_galaxy, dec_galaxy, PA_deg, a_arcsec, b_arcsec, name=sn_name, color=ellipse_color, width=ellipse_width))
    cmds.append(make_circle_command_for_position(fits_summary.at[index, "RA_SN"], fits_summary.at[index, "Dec_SN"], name=sn_name))

    return cmds, ra_galaxy, dec_galaxy



def fpack_download(fitsfilename, regionfilename, clobber=False, scale=100.0):
    # ensure absolute paths
    #fitsfilename = os.path.abspath(fitsfilename)
    #regionfilename = os.path.abspath(regionfilename)

    # build output names (keep your style)
    outfile = fitsfilename[:-5] + '_sci_int32.fits'   # e.g. *_drz_sci_int32.fits
    out_fz  = outfile + '.fz'                         # after fpack

    # --- read SCI image (robust selection) ---
    with fits.open(fitsfilename, memmap=False) as hdul:
        try:
            sci = hdul['SCI', 1]            # preferred
        except Exception:
            if len(hdul) > 1:
                sci = hdul[1]               # common DRZ/DRC case
            else:
                sci = hdul[0]               # single-HDU fallback
        data = np.array(sci.data, dtype=np.float32, copy=True)
        hdr  = sci.header.copy()

    # --- scale to int32 with BSCALE so DS9 still sees floats ---
    #int_data = np.rint(data * scale).astype(np.int32)
    phdu = fits.PrimaryHDU()
    hdu_img = fits.ImageHDU(data, header=hdr)
    hdu_img.header['BITPIX'] = 32
    hdu_img.header['BZERO']  = 0.0
    hdu_img.header['BSCALE'] = 1.0/scale
    fits.HDUList([phdu, hdu_img]).writeto(outfile, overwrite=True)

    if os.path.exists(outfile):
        print('file created successfully:', outfile)
   
       # 2) If .fz already exists, decide whether to reuse or clobber
    if os.path.exists(out_fz):
        if clobber:
            try:
                os.unlink(out_fz)
            except OSError:
                pass
        else:
            # Reuse existing compressed file
            # Clean up the temporary uncompressed file and return
            try:
                os.unlink(outfile)
            except OSError:
                pass
            return out_fz, regionfilename

    # 3) Run fpack via the dedicated env
    try:
        subprocess.run(
            ["conda", "run", "-n", "fpack", "fpack", "-Y", outfile],
            check=True
        )
    except subprocess.CalledProcessError as e:
        # If fpack complained but the .fz exists anyway, treat as success
        if os.path.exists(out_fz):
            print(f"[fpack] non-zero exit ({e.returncode}) but {out_fz} exists; using it.")
        else:
            # real failure -> re-raise so caller can skip this file
            raise

    # 4) Remove the temporary uncompressed file
    try:
        os.unlink(outfile)
    except OSError:
        pass

    # return both files so the caller can add to the manifest
    return out_fz, regionfilename






def fpack_download2(fitsfilename, relative_regionfilename):
    # 1) Read SCI HDU (usually ext=1), keep its header (contains WCS)
    print(fitsfilename, 'fitsfilename')
    outfile = fitsfilename[:-5]+'sci.fits'
    out_fz = fitsfilename[:-5]+'sci.fz'

    with fits.open(fitsfilename, memmap=False) as hdul:
        sci = hdul[1]    # or hdul['SCI',1]
        print('sci data:', sci)
        data = sci.data.astype(np.float32).copy()
        hdr  = sci.header.copy()

    # 2) Choose a scale (here 100 -> 0.01 precision)
    scale = 100.0
    int_data = np.round(data * scale).astype(np.int32)

    # 3) Write as integer with BSCALE/BZERO so viewers recover floats = int*BSCALE+BZERO
    hdu_img = fits.ImageHDU(int_data, header=hdr)
    hdu_img.header['BITPIX'] = 32
    hdu_img.header['BZERO']  = 0.0
    hdu_img.header['BSCALE'] = 1.0/scale
    # (optional) clarify units
    if 'BUNIT' in hdu_img.header:
        hdu_img.header['BUNIT'] = hdu_img.header['BUNIT']  # unchanged physically, scaling is via BSCALE

    # Add a minimal PRIMARY HDU (no data)
    phdu = fits.PrimaryHDU()
    fits.HDUList([phdu, hdu_img]).writeto(outfile, overwrite=True)
    if outfile.exists():
        print('file created succesfully' , outfile)
    if shutil.which("fpack") is None:
        raise RuntimeError("fpack not found in PATH on the server")
    subprocess.rn(["fpack", "-Y", str(outfile)], check =True)

    #to delete the science image
    if out_fz.exists():
        try: outfile.unlin()
        except Exception: pass

    return str(out_fz)


if __name__ == "__main__":

    args = define_args()
    images_to_download=[]


    # pandas reads the table at the path args.fits_summary_path and then returns it to the variable fits_summary
    print (f'Loading fits summary: {args.fits_summary_path}')
    fits_summary = pd.read_table(args.fits_summary_path,sep=',') #original


    if args.supernova_names[0].lower() == 'all':
        supernova_names = sorted(fits_summary['SNID'].unique())
    else:
        supernova_names = args.supernova_names

    print(f'supernova_names: {supernova_names}')

    # reading table at args.probable_galaxy_path and assigning it to probable_galaxy variable
    print (f'Loading probable galaxy: {args.probable_galaxy_path}')
    probable_galaxy = pd.read_table(args.probable_galaxy_path,sep=',') #original

    # reading table at args.probable_galaxy_path and assigning it to probable_galaxy variable
    print (f'Loading unmatched_pantheon: {args.unmatched_pantheon_path}')
    unmatched_pantheon = pd.read_table(args.unmatched_pantheon_path,sep=',') #original
    fits_summary = combine_fs_and_up(supernova_names, fits_summary, unmatched_pantheon)

    # receiving end: the "filenameshort" column in the fits_summary table
    fits_summary['filenameshort'] = fits_summary['filename'].str.replace(r'.*/', '', regex=True)

    # combines fits_summary and probable_galaxy tables
    print("Combining fits_summary and probable_galaxy tables...")
    fits_summary = combine_fs_and_bg(supernova_names, fits_summary, probable_galaxy)

    #print(fits_summary.head().to_string())
    fits_summary_combined_filename = re.sub('\.csv$','_combined.csv',args.fits_summary_path)
    if fits_summary_combined_filename == args.fits_summary_path:
        raise RuntimeError(f'EERRRRRRROORRRR: {fits_summary_combined_filename}')

    get_region_info2(fits_summary)
    print(f'Saving combined file into {fits_summary_combined_filename}')
    fits_summary.to_csv(fits_summary_combined_filename)

    pd.set_option('display.max_colwidth', None)
    



    if not os.path.exists(args.out_dir):
        os.mkdir(args.out_dir)

    ds9cmds_list = []

    for sn_name in supernova_names:
        print("\nSN name: ",sn_name)

        sn_ix = np.where(fits_summary["SNID"] == sn_name)[0]

        #ra_galaxy = None # fits_summary['RA_GALAXY][sn_ix]
        #dec_galaxy = None ## fits_summary['Dec_GALAXY][sn_ix]

        all_output = ['global color=blue dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1']
        all_output.append('fk5')
        
        #sn_rootdir = f'{args.datadir}/source_extractor/{sn_name}'
        sn_rootdir = f'{args.datadir}/pantheon_data_folder/{sn_name}/mastDownload/HST'

        ds9cmd = f'\ncd {sn_rootdir} \n ds9 -zscale '

        # if sn_name=='2010ai':

        #     print(fits_summary.loc[sn_ix,['File_key','Position_Angle','Major_Axis','Minor_Axis','a_pixels','b_pixels','theta_deg','PA_deg','a_arcsec','b_arcsec']])
        #     sys.exit(0)


        #CHANGE TO GREEN 

        green_ellipse_list=[]
        for index in sn_ix:

            if fits_summary.loc[index,"Telescope"]=='JWST':
                print(f'WARNING! skipping since telescope={fits_summary.loc[index,"Telescope"]}')
                continue

            if fits_summary.loc[index,"Filter"]=='detection':
                print(f'WARNING! skipping since filter={fits_summary.loc[index,"Filter"]}')
                continue
                

            #print(f'OH DEAR!',fits_summary.columns)
            #print(f'Telescope: {fits_summary.loc[index,"Telescope"]}, Filter:{fits_summary.loc[index,"Filter"]}')
            #sys.exit(0)

            cmds, ra_galaxy, dec_galaxy = get_region_commands(index, fits_summary, ellipse_color='green', ellipse_width=1, sn_name=None)
            if (cmds is None) or (cmds is []):
                continue

            green_ellipse_list.append(cmds[0])
            # print (f' GIRL I HOPE THIS WORKS {green_ellipse_list}')
            # print (f'LOOOOOOOOOOK{cmds}')
            # print(f'HAHAHHAHAHHAHAHAH anyway... {cmds[0]}')
            

        print(f'PRINTING THE ELLIPSE LIST THAT WILL BE IN GREEN... {green_ellipse_list}')
        
        

        for index in sn_ix:

            if fits_summary.loc[index,"Telescope"]=='JWST':
                print(f'WARNING! skipping since telescope={fits_summary.loc[index,"Telescope"]}')
                continue

            if fits_summary.loc[index,"Filter"]=='detection':
                print(f'WARNING! skipping since filter={fits_summary.loc[index,"Filter"]}')
                continue
                

            #print(f'OH DEAR!',fits_summary.columns)
            #print(f'Telescope: {fits_summary.loc[index,"Telescope"]}, Filter:{fits_summary.loc[index,"Filter"]}')
            #sys.exit(0)

            cmds, ra_galaxy, dec_galaxy = get_region_commands(index, fits_summary, ellipse_width=4, sn_name=sn_name)
            if (cmds is None) or (cmds is []):
                continue

            # make region file
            output = ['global color=blue dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1']
            output.append('fk5')

            for cmd in cmds:
                output.append(cmd)
                all_output.append(cmd)

            output.extend(green_ellipse_list)
            # sys.exit(0)

            #ADD HERE

            #fitsfilename = f'{sn_rootdir}/{fits_summary.at[index, "File_key"]}/{fits_summary.at[index, "File_key"]}.fits'

            (path1,fname1) = os.path.split(fits_summary.at[index, "filename"])
            (path2,subdir1) = os.path.split(path1)
            fitsfilename = f'{sn_rootdir}/{subdir1}/{fname1}'

            #print('bbbb',fitsfilename_old,'\n',fitsfilename)
            #sys.exit(0)

            print(f'finding directory for {fitsfilename}')
            if os.path.isdir(fitsfilename):
                print(f'fitsfilename {fitsfilename} directory exists!!!')
            else:
                if os.path.isfile(fitsfilename):
                    print(f'fitsfilename {fitsfilename} file exists!!!')
                else:
                    print(f'ERROR: fitsfilename {fitsfilename} does not exist!!!')


            regionfilename = re.sub('\.fits$','.reg',fitsfilename)
            if regionfilename==fitsfilename:
                raise RuntimeError(f'Could not replace .fits with .reg in {fitsfilename}')
            
            #filepath_regionfile = f"{args.out_dir}/{fits_summary.at[index, 'File_key']}_region_file.reg"

            print('fitsfilename',fitsfilename)
            print('regionfilename',regionfilename)

            if not os.path.isfile(fitsfilename):
                print(f'WARNING!!!!! fitsfile {fitsfilename} does not exist!!!')
                #sys.exit(0)
                continue

            relative_fitsfilename = os.path.relpath(fitsfilename, start=sn_rootdir)
            relative_regionfilename = os.path.relpath(regionfilename, start=sn_rootdir)

            ds9cmd += f' {relative_fitsfilename} -regionfile {relative_regionfilename}'
            print(ds9cmd)

            save_region_file(regionfilename, output)

            if args.download:
                try:
                    print('print trying to download',fitsfilename, regionfilename)
                    packed_path, reg_path = fpack_download(fitsfilename, regionfilename)
                    images_to_download.extend([packed_path, reg_path])
                    print(f"[download] created {packed_path}")
                    print(f"[download] include region {reg_path}")

                except Exception as e:
                    import traceback; traceback.print_exc()
                    print(f"[download] skipping {fitsfilename}: {e}")
                    print(f"images failed to download for {relative_fitsfilename}:{e}")
        
        if sn_name == '1997bq':
            print(fits_summary.loc[sn_ix,['File_key','theta_deg','alpha_deg','PA_deg','pixscale']])
        print(ds9cmd)


        ds9cmds_list.append(ds9cmd)

    filepath = f"{args.out_dir}/{args.all_output_filename}"
    save_region_file(filepath, all_output)

    ds9cmd_filename = f"{args.out_dir}/ds9cmds.txt"
    ds9cmds = "\n".join(ds9cmds_list)

    print (f'writing ds9 commands to {ds9cmd_filename}')

    with open(ds9cmd_filename, "w") as f:
        f.write(ds9cmds)
    f.close()
     
    print(images_to_download, 'images to donwload') 
    if args.download and images_to_download:
        manifest = os.path.join(args.out_dir, "download_manifest.txt")
        print(manifest, 'manifest')
        with open(manifest, "w") as f:
             for p in images_to_download:
                f.write(os.path.abspath(p) + "\n")
        print(f"\n[download] Wrote manifest: {manifest}")

        # Ready-to-paste commands on YOUR LAPTOP terminal:
        print("\n=== COPY THESE TO YOUR LAPTOP TERMINAL ===")
        print("# Option A: rsync using a file list (fast, keeps names)")
        print('mkdir -p ~/HST_quicklook')
        print(f'rsync -av --info=progress2 -e "ssh -J aperturesshduo.stsci.edu:8222" '
          f'--files-from="{manifest}" plhstproc6.stsci.edu:/ ~/HST_quicklook/')
        print("\n# Option B: loop over files (if rsync with --files-from is blocked):")
        print('while IFS= read -r p; do '
          'scp -o ProxyJump=aperturesshduo.stsci.edu:8222 '
          '"plhstproc6.stsci.edu:${p}" ~/HST_quicklook/ ; '
          f'done < "{manifest}"')
        print("=== END ===\n")

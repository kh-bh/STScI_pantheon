#!/usr/bin/env python
'''
Input:
1. FITS image
2. Coordinates of SN on image (in ra,dec u.deg)

Looking to get
1. Probability of host for any nearby sources
2. d_DLR for same nearby sources
'''
import os,sys,pdb,scipy,glob
import pandas as pd
from pylab import *
from scipy.integrate import quad
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.io import ascii
import pandas as pd



def read_sexcat(filename, cmtchar='#'):
    col_names =[]
    with open(filename, 'r') as f:
        for i, line in enumerate(f):
            if line.startswith(cmtchar):
                col_names.append(line.split()[2])
    df = pd.read_csv(filename, comment=cmtchar, header=None, sep=r'\s+')
    df.columns=col_names
    return(df)

def make_regfile(regfile,xx, yy, df,nclose=6):
    ## select only the sources with d_DLR < 5
    df=df[df['d_DLR']<4]
    
    f=open(regfile, 'w')
    f.write('global color=green dashlist=8 3 width=1 ')
    f.write('font="helvetica 10 normal roman" select=1 ')
    f.write('highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
    f.write('image\n')
    for ii, rr in df.iloc[:nclose].iterrows():
        f.write('ellipse(%d,%d,%d,%d,%d)\n' %(rr['X_IMAGE'],
                                              rr['Y_IMAGE'],
                                              rr['A_IMAGE']*rr['KRON_RADIUS'],
                                              rr['B_IMAGE']*rr['KRON_RADIUS'],
                                              rr['THETA_IMAGE'])
                )
        f.write('# text(%.1f,%.1f) textangle=0, text={%d}\n'
                %(rr['X_IMAGE']+10, rr['Y_IMAGE']+10, rr['NUMBER']))
        f.write('point (%.1f, %.1f) # point=cross color=red\n' %(xx, yy))
    f.close()
    return()

def dv(x,*p):
    k,n=p
    I_Ie=exp(-k*(x**(1./n)-1.0))
    return(I_Ie)

def prob_of_gxy(offset,n):
    p0 = [7.67, n]
    frac_of_light = quad(dv, 0, offset, args=tuple(p0))[0]/quad(dv, 0, 1e3, args=tuple(p0))[0]
    return(frac_of_light)

def create_prob_index(root_dir, output_csv, top_range = 0, bottom_range = 0):
    """
    Create a datatable with the MOST PROBABLE objects of each source extractor file
    Parameters:
        root_dir (str): The root directory to search.

    Returns:
        dict: A dictionary where keys are relative paths to .cat files,
              and values are dictionaries of the first data row.
    """
    valid_filter = pd.read_csv('/home/bkhatri/STScI_pantheon_Github/data_collection_files/data_files/fits_cleaned_HST_data.csv')
    imfile = valid_filter['filename']
    SNID = valid_filter['SNID']
    file_key = valid_filter['File_key']
    data_dict = []
    sn1_data = pd.read_csv('/home/bkhatri/STScI_pantheon_Github/data_collection_files/data_files/unmatched_pantheon.csv')
    i = 0

    if not os.path.exists(output_csv):
        try:
            with open(output_csv, 'x') as f:
                f.write("This is a new file.")
            print(f"File '{output_csv}' created successfully.")
        except FileExistsError:
            # This block is theoretically not reached if os.path.exists is accurate,
            # but provides a robust fallback for race conditions.
            print(f"File '{output_csv}' already exists")
    else:
        print(f"File '{output_csv}' already exists, ignoring creation.")




    if top_range == 0 and bottom_range == 0:
        bottom_range = len(file_key)    
    
    for index in range(top_range, bottom_range):
        rr = sn1_data.loc[sn1_data['SNID'] == SNID[index], 'RA'].values[0]
        dd = sn1_data.loc[sn1_data['SNID'] == SNID[index], 'Dec'].values[0]
        with fits.open(imfile[index]) as hdul:
            # Get the WCS object
            w = WCS(hdul[1].header)
            coords = SkyCoord(rr, dd, unit='deg')
            # Convert to pixel coordinates
            xx, yy = w.world_to_pixel(coords)
        
        path = glob.glob(os.path.join(root_dir, SNID[index], file_key[index], '*.cat'), recursive=True)
        file_path = path[0]
        df = read_sexcat(file_path)


        to_rad = pi/180.0
        Dx = xx-df['X_IMAGE']
        Dy = yy-df['Y_IMAGE']
        df['OFFSET']=sqrt(Dx**2+Dy**2)
        phi = arctan(Dy/Dx)-to_rad*df['THETA_IMAGE'] ##in rad
        xp = df['A_IMAGE']*cos(phi)*cos(to_rad*df['THETA_IMAGE'])-df['B_IMAGE']*sin(phi)*sin(to_rad*df['THETA_IMAGE'])
        yp = df['A_IMAGE']*cos(phi)*sin(to_rad*df['THETA_IMAGE'])-df['B_IMAGE']*sin(phi)*cos(to_rad*df['THETA_IMAGE'])
        df['RE']= df['KRON_RADIUS']*sqrt(xp**2.+yp**2.)
        df['d_DLR']=df['OFFSET']/df['RE']

        df=df.sort_values(by='OFFSET', ignore_index=True)

        for nn in list([1,2,4]):
            f_hosts=[]
            for item in list(df['d_DLR']):
                tmp = 1.0-round(prob_of_gxy(item, n=nn), 5)
                if tmp < 0: pdb.set_trace()
                f_hosts.append(tmp)
            df['F_HOST_%d'%nn]=f_hosts
        ## pdb.set_trace()
        df=df.sort_values(by='F_HOST_4', ignore_index=True, ascending=False)
        probable = 0
        info = {
                "SNID": SNID[index],
                "File_key": file_key[index], 
                "RA_Galaxy": df['ALPHA_J2000'][probable],
                "Dec_Galaxy": df['DELTA_J2000'][probable],
                "Position_Angle": df['THETA_IMAGE'][probable],
                "Major_Axis": df['A_IMAGE'][probable],
                "Minor_Axis": df['B_IMAGE'][probable],
                "Kron_Radius": df['KRON_RADIUS'][probable],
                "Mag": df['MAG_AUTO'][probable],
                "Mag_Error": df['MAGERR_AUTO'][probable],
                "Flux": df['FLUX_AUTO'][probable],
                "Flux_Error": df['FLUXERR_AUTO'][probable],
                "ISOAREA_IMAGE": df['ISOAREA_IMAGE'][probable],
                "X_IMAGE": df['X_IMAGE'][probable],
                "Y_IMAGE": df['Y_IMAGE'][probable],
                "CXX": df['CXX_IMAGE'][probable],
                "CXY": df['CXY_IMAGE'][probable],
                "CYY":df['CYY_IMAGE'][probable],
                "OFFSET": df['OFFSET'][probable],
                "RE": df['RE'][probable],
                "d_DLR": df['d_DLR'][probable],
                "F_HOST_4": df['F_HOST_4'][probable]
                }
        i = i+1

        if i%10 == 0:
            df = pd.DataFrame(data_dict)
            df.to_csv('most_probable_gal.csv', mode = 'a', index=False)
        data_dict.append(info)

    


    
if __name__=='__main__':
    tr = int(sys.argv[1])
    br = int(sys.argv[2])
    source_extractor_path = "/astro/armin/bhoomika/source_extractor"
    #root_dir = os.path.join(home_dir, 'source_extractor')
    #print(root_dir)
    relative_key = 0
    create_prob_index(source_extractor_path, 'most_probable_gal.csv',tr, br)
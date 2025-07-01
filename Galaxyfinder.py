#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys, os, glob
import argparse
import pandas as pd

# ds9:
#ellipse(186.8194200,9.4313379,10.000",30.000",44.999994) # color=red width=4
#ellipse(RA, Dec, minor axis lenght, major axis length, Poisition Angle PA in deg)
def mk_regioncmd_for_position(ra,dec,position_angle,major_axis,minor_axis,name=None,symbol='circle',color='red',width=4,size='10"'):
	print(f'hello RA_Galaxy:{ra} Dec_Galaxy:{dec} name:{name} Position_Angle{position_angle} Major_Axis{major_axis} Minor_Axis{minor_axis}')
	name='{'+f'{name}'+'}'
	s = f'circle({ra},{dec},{size},{position_angle},{major_axis},{minor_axis}) # color={color} width={width} text={name}'
	print(s)
	return(s)

if __name__ == '__main__':
	parser = argparse.ArgumentParser(conflict_handler='resolve')
	parser.add_argument('galaxytablefilename', type=str, help=f'Specify the filename of the SN table.')
	parser.add_argument('regionfile', type=str, help=f'Specify the region output filename.')
	parser.add_argument('-c','--color', type=str, default='green', choices = ['red','green','magenta'],help=f'Specify the color.')
	parser.add_argument('-s','--size_arcsec', type=float, default=5.0,help=f'Specify the size in arcsec')
	args = parser.parse_args()

	print(f'galaxytablefilename: {args.galaxytablefilename}')
	print(f'regionfile: {args.regionfile}')

	mytable = pd.read_table(args.galaxytablefilename,sep=',')
	#mytable = pd.read_table(args.SNtablefilename)
	print(mytable)
	print(mytable.columns)

	print('VVVV')
	ixs = [1,2,3]
	columns = ['SNID', 'File_key', 'Galaxy_Key', 'RA_Galaxy', 'Dec_Galaxy',
	   'Position_Angle', 'Major_Axis', 'Minor_Axis', 'Mag', 'Mag_Error',
	   'Flux', 'Flux_Error', 'CXX_IMAGE','CYY_IMAGE','CXY_IMAGE']
	print(mytable.loc[ixs,columns])
	sys.exit(0)

	output = ['global color=blue dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1']
	output.append('fk5')

	for ix, row in df.iterrows():
		ra = row['RA_Galaxy']
		dec = row['Dec_Galaxy']
		name = row['SNID'] # name = row.get('SNID', f'Obj{ix}')

		# Get covariance matrix elements
		cxx = row['CXX_IMAGE']
		cyy = row['CYY_IMAGE']
		cxy = row['CXY_IMAGE']

		# Eigenvalues
		trace = cxx + cyy
		delta = (cxx - cyy) / 2
		root_term = math.sqrt(delta**2 + cxy**2)
		lambda1 = (trace / 2) + root_term
		lambda2 = (trace / 2) - root_term
		
		# If fixed size requested
		if args.size_arcsec:
			a = b = args.size_arcsec
		else:
			# Axes lengths (in pixels; assuming pixel scale is 1 arcsec/pix unless specified)
			a = math.sqrt(max(lambda1, lambda2))
			b = math.sqrt(min(lambda1, lambda2))

		# Position angle (in degrees)
		theta_rad = 0.5 * math.atan2(2 * cxy, cxx - cyy)
		theta_deg = math.degrees(theta_rad)

		region = mk_regioncmd_for_position(ra, dec, theta_deg, a, b, name=name, color=args.color)
		output.append(region)

	print(f'Writing region file to {args.regionfile}')
	with open(args.regionfile, "w") as file:
		for line in output:
			file.write(str(line) + "\n")
	file.close()	



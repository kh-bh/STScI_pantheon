#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys, os, glob
import argparse
import pandas as pd

def mk_regioncmd_for_position(ra,dec,name=None,symbol='circle',color='red',width=4,size='10"'):
	print(f'hello ra:{ra} dec:{dec} name:{name}')
	name='{'+f'{name}'+'}'
	s = f'circle({ra},{dec},{size}) # color={color} width={width} text={name}'
	print(s)
	return(s)

if __name__ == '__main__':
	parser = argparse.ArgumentParser(conflict_handler='resolve')
	parser.add_argument('SNtablefilename', type=str, help=f'Specify the filename of the SN table.')
	parser.add_argument('regionfile', type=str, help=f'Specify the region output filename.')
	parser.add_argument('-c','--color', type=str, default='magenta', choices = ['red','green','magenta'],help=f'Specify the color.')
	parser.add_argument('-s','--size_arcsec', type=float, default=5.0,help=f'Specify the size in arcsec')
	args = parser.parse_args()

	print(f'FFFFF {args.SNtablefilename}')

	a = 9
	print(f'HELLO FAE FAE {a}')

	mytable = pd.read_table(args.SNtablefilename,sep=',')
	#mytable = pd.read_table(args.SNtablefilename)
	print(mytable)
	print(mytable.columns)

	print('VVVV')
	ixs = [1,2,3]
	columns = ['RA','Dec']
	print(mytable.loc[ixs,columns])

	output = ['global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1']
	output.append('fk5')

	ixs = mytable.index.values
	for ix in ixs:
		print(mytable.loc[ix,'SNID'])

		output.append(mk_regioncmd_for_position(mytable.loc[ix,'RA'],
			mytable.loc[ix,'Dec'],name=mytable.loc[ix,'SNID'],
			color=args.color,
			size=f'{args.size_arcsec}"'))

	print(f'Writing region file to {args.regionfile}')
	with open(args.regionfile, "w") as file:
	    for line in output:
	        file.write(str(line) + "\n")
	file.close()	



#!/usr/bin/python

#===========================================================
#Copyright(c)2015, IBMC, CNRS
#All rights reserved.
#NAME:		format_check.py
#ABSTRACT:	
#DATE:		Mon Oct 17 20:09:56 2016
#Usage:		
#VERSION: 	0.01
#AUTHOR: 	Miao Zhichao
#CONTACT: 	chichaumiau AT gmail DOT com
#NOTICE: This is free software and the source code is freely
#available. You are free to redistribute or modify under the
#conditions that (1) this notice is not removed or modified
#in any way and (2) any modified versions of the program are
#also available for free.
#		** Absolutely no Warranty **
#===========================================================

#Input: A PDB file that need to be checked and a reference PDB of good format
#$python format_check.py 2gdi.pdb 2gdi.out.pdb

import sys,os
from Bio.PDB import *

#fpdb is the file need to be checked
#fref is the file of good format generated by rna_puzzles_format.py
def format_check(fpdb,fref):
	parser= PDBParser(PERMISSIVE=1)
	try:
		struct_ref=parser.get_structure('SI',fref)
	except Exception:
		sys.stderr.write('ERROR: structure file format error <%s>\n'%fref)
		print(Exception)
		return False
	try:
		struct=parser.get_structure('SI',fpdb)
	except Exception:
		sys.stderr.write('ERROR: structure file format error <%s>\n'%fpdb)
		print(Exception)
		return False
	atm_list=[]
	for chn in struct_ref[0]:
		for res in chn:
			for atm in res:
				atm_list.append('%s%s%s'%(chn.id,res.__repr__(),atm.__repr__()))
	atm_map={}
	for i,atm in enumerate(atm_list):
		atm_map[atm]=i
	
	n_err=0
	err_log=''
	if len(struct.child_list)<1:
		open('%s.format_check.txt'%fpdb,'w').write("0 MODEL in the submission!")
		return False
	for i,mdl in enumerate(struct.child_list):
		if i >=5:break
		flag=[False for x in range(len(atm_list))]
		for chn in mdl:
			for res in chn:
				for atm in res:
					name='%s%s%s'%(chn.id,res.__repr__(),atm.__repr__())
					n=atm_map.get(name,-1)
					if n>-1:
						flag[n]=True
		for j,fg in enumerate(flag):
			if not fg:
				n_err+=1
				err_log+='ERROR: ATOM [%s] not found in model %d\n'%(atm_list[j],i)
	
	if n_err >0:
		open('%s.format_check.txt'%fpdb,'w').write(err_log)
		return False
	return True

if __name__ == '__main__':
	if format_check(sys.argv[1],sys.argv[2]):
		print('Congratulations! Format check passed!')
	else:
		print('Format error! Please check again and re-submit!')

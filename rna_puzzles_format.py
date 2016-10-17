#!/usr/bin/python

#===========================================================
#Copyright(c)2013, IBMC, CNRS
#All rights reserved.
#NAME:		rnatemplate.py
#ABSTRACT:	input a RNA sequence (fasta format), output the standard PDB format
#DATE:		Tue Sep 17 15:08:40 2013
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


import sys

Usage="""rnatemplate.py usage:

input a RNA sequence (fasta format), output the standard PDB format

./rnatemplate.py fasta.file number_of_model(optional) >output.pdb
fasta.file example:
>RNA1 A length1
UGCGAUGAGAAGAAGAGUAUUAAGGAUUUACUAUGAUUAGCGACUCUAGGAUAGUGAAAG
     CUAGAGGAUAGUAACCUUAAGAAGGCACUUCGAGCA
>RNA2 B length2
GCGGAAGUAGUUCAGUGGUAGAACACCACCUUGCCAAGGUGGGGGUCGCGGGUUCGAAUC
     CCGUCUUCCGCUCCA
"""

A_temp="""ATOM  %5d  P     A %c%4d       0.000   0.000   0.000  1.00  0.00           P
ATOM  %5d  OP1   A %c%4d       0.000   0.000   0.000  1.00  0.00           O
ATOM  %5d  OP2   A %c%4d       0.000   0.000   0.000  1.00  0.00           O
ATOM  %5d  O5'   A %c%4d       0.000   0.000   0.000  1.00  0.00           O
ATOM  %5d  C5'   A %c%4d       0.000   0.000   0.000  1.00  0.00           C
ATOM  %5d  C4'   A %c%4d       0.000   0.000   0.000  1.00  0.00           C
ATOM  %5d  O4'   A %c%4d       0.000   0.000   0.000  1.00  0.00           O
ATOM  %5d  C3'   A %c%4d       0.000   0.000   0.000  1.00  0.00           C
ATOM  %5d  O3'   A %c%4d       0.000   0.000   0.000  1.00  0.00           O
ATOM  %5d  C2'   A %c%4d       0.000   0.000   0.000  1.00  0.00           C
ATOM  %5d  O2'   A %c%4d       0.000   0.000   0.000  1.00  0.00           O
ATOM  %5d  C1'   A %c%4d       0.000   0.000   0.000  1.00  0.00           C
ATOM  %5d  N9    A %c%4d       0.000   0.000   0.000  1.00  0.00           N
ATOM  %5d  C8    A %c%4d       0.000   0.000   0.000  1.00  0.00           C
ATOM  %5d  N7    A %c%4d       0.000   0.000   0.000  1.00  0.00           N
ATOM  %5d  C5    A %c%4d       0.000   0.000   0.000  1.00  0.00           C
ATOM  %5d  C6    A %c%4d       0.000   0.000   0.000  1.00  0.00           C
ATOM  %5d  N6    A %c%4d       0.000   0.000   0.000  1.00  0.00           N
ATOM  %5d  N1    A %c%4d       0.000   0.000   0.000  1.00  0.00           N
ATOM  %5d  C2    A %c%4d       0.000   0.000   0.000  1.00  0.00           C
ATOM  %5d  N3    A %c%4d       0.000   0.000   0.000  1.00  0.00           N
ATOM  %5d  C4    A %c%4d       0.000   0.000   0.000  1.00  0.00           C
"""

G_temp="""ATOM  %5d  P     G %c%4d       0.000   0.000   0.000  1.00  0.00           P
ATOM  %5d  OP1   G %c%4d       0.000   0.000   0.000  1.00  0.00           O
ATOM  %5d  OP2   G %c%4d       0.000   0.000   0.000  1.00  0.00           O
ATOM  %5d  O5'   G %c%4d       0.000   0.000   0.000  1.00  0.00           O
ATOM  %5d  C5'   G %c%4d       0.000   0.000   0.000  1.00  0.00           C
ATOM  %5d  C4'   G %c%4d       0.000   0.000   0.000  1.00  0.00           C
ATOM  %5d  O4'   G %c%4d       0.000   0.000   0.000  1.00  0.00           O
ATOM  %5d  C3'   G %c%4d       0.000   0.000   0.000  1.00  0.00           C
ATOM  %5d  O3'   G %c%4d       0.000   0.000   0.000  1.00  0.00           O
ATOM  %5d  C2'   G %c%4d       0.000   0.000   0.000  1.00  0.00           C
ATOM  %5d  O2'   G %c%4d       0.000   0.000   0.000  1.00  0.00           O
ATOM  %5d  C1'   G %c%4d       0.000   0.000   0.000  1.00  0.00           C
ATOM  %5d  N9    G %c%4d       0.000   0.000   0.000  1.00  0.00           N
ATOM  %5d  C8    G %c%4d       0.000   0.000   0.000  1.00  0.00           C
ATOM  %5d  N7    G %c%4d       0.000   0.000   0.000  1.00  0.00           N
ATOM  %5d  C5    G %c%4d       0.000   0.000   0.000  1.00  0.00           C
ATOM  %5d  C6    G %c%4d       0.000   0.000   0.000  1.00  0.00           C
ATOM  %5d  O6    G %c%4d       0.000   0.000   0.000  1.00  0.00           O
ATOM  %5d  N1    G %c%4d       0.000   0.000   0.000  1.00  0.00           N
ATOM  %5d  C2    G %c%4d       0.000   0.000   0.000  1.00  0.00           C
ATOM  %5d  N2    G %c%4d       0.000   0.000   0.000  1.00  0.00           N
ATOM  %5d  N3    G %c%4d       0.000   0.000   0.000  1.00  0.00           N
ATOM  %5d  C4    G %c%4d       0.000   0.000   0.000  1.00  0.00           C
"""

U_temp="""ATOM  %5d  P     U %c%4d       0.000   0.000   0.000  1.00  0.00           P
ATOM  %5d  OP1   U %c%4d       0.000   0.000   0.000  1.00  0.00           O
ATOM  %5d  OP2   U %c%4d       0.000   0.000   0.000  1.00  0.00           O
ATOM  %5d  O5'   U %c%4d       0.000   0.000   0.000  1.00  0.00           O
ATOM  %5d  C5'   U %c%4d       0.000   0.000   0.000  1.00  0.00           C
ATOM  %5d  C4'   U %c%4d       0.000   0.000   0.000  1.00  0.00           C
ATOM  %5d  O4'   U %c%4d       0.000   0.000   0.000  1.00  0.00           O
ATOM  %5d  C3'   U %c%4d       0.000   0.000   0.000  1.00  0.00           C
ATOM  %5d  O3'   U %c%4d       0.000   0.000   0.000  1.00  0.00           O
ATOM  %5d  C2'   U %c%4d       0.000   0.000   0.000  1.00  0.00           C
ATOM  %5d  O2'   U %c%4d       0.000   0.000   0.000  1.00  0.00           O
ATOM  %5d  C1'   U %c%4d       0.000   0.000   0.000  1.00  0.00           C
ATOM  %5d  N1    U %c%4d       0.000   0.000   0.000  1.00  0.00           N
ATOM  %5d  C2    U %c%4d       0.000   0.000   0.000  1.00  0.00           C
ATOM  %5d  O2    U %c%4d       0.000   0.000   0.000  1.00  0.00           O
ATOM  %5d  N3    U %c%4d       0.000   0.000   0.000  1.00  0.00           N
ATOM  %5d  C4    U %c%4d       0.000   0.000   0.000  1.00  0.00           C
ATOM  %5d  O4    U %c%4d       0.000   0.000   0.000  1.00  0.00           O
ATOM  %5d  C5    U %c%4d       0.000   0.000   0.000  1.00  0.00           C
ATOM  %5d  C6    U %c%4d       0.000   0.000   0.000  1.00  0.00           C
"""

C_temp="""ATOM  %5d  P     C %c%4d       0.000   0.000   0.000  1.00  0.00           P
ATOM  %5d  OP1   C %c%4d       0.000   0.000   0.000  1.00  0.00           O
ATOM  %5d  OP2   C %c%4d       0.000   0.000   0.000  1.00  0.00           O
ATOM  %5d  O5'   C %c%4d       0.000   0.000   0.000  1.00  0.00           O
ATOM  %5d  C5'   C %c%4d       0.000   0.000   0.000  1.00  0.00           C
ATOM  %5d  C4'   C %c%4d       0.000   0.000   0.000  1.00  0.00           C
ATOM  %5d  O4'   C %c%4d       0.000   0.000   0.000  1.00  0.00           O
ATOM  %5d  C3'   C %c%4d       0.000   0.000   0.000  1.00  0.00           C
ATOM  %5d  O3'   C %c%4d       0.000   0.000   0.000  1.00  0.00           O
ATOM  %5d  C2'   C %c%4d       0.000   0.000   0.000  1.00  0.00           C
ATOM  %5d  O2'   C %c%4d       0.000   0.000   0.000  1.00  0.00           O
ATOM  %5d  C1'   C %c%4d       0.000   0.000   0.000  1.00  0.00           C
ATOM  %5d  N1    C %c%4d       0.000   0.000   0.000  1.00  0.00           N
ATOM  %5d  C2    C %c%4d       0.000   0.000   0.000  1.00  0.00           C
ATOM  %5d  O2    C %c%4d       0.000   0.000   0.000  1.00  0.00           O
ATOM  %5d  N3    C %c%4d       0.000   0.000   0.000  1.00  0.00           N
ATOM  %5d  C4    C %c%4d       0.000   0.000   0.000  1.00  0.00           C
ATOM  %5d  N4    C %c%4d       0.000   0.000   0.000  1.00  0.00           O
ATOM  %5d  C5    C %c%4d       0.000   0.000   0.000  1.00  0.00           C
ATOM  %5d  C6    C %c%4d       0.000   0.000   0.000  1.00  0.00           C
"""

#=======================================================================
#Read fasta file
#Input: fasta file
#format:
#~ >name1 chain1 length1
#~ sequence1
#~ >name2 chain2 length2
#~ sequence2
#~ ........
def readfasta(fp):
	name='xxx'
	chains=[]
	seqs=[]
	seq=''
	f=file(fp)
	for line in f:
		if len(line)<2:continue
		if line[0] == '#':continue
		if line[0] == '>':
			a=line.strip().split()
			if( len(a) < 2 ):
				print Usage
				exit(0)
			name=a[0][1:]
			chains.append(a[1][0])
			if len(seq)>0:
				seqs.append(seq)
			seq=''
		else:
			seq+=line.strip().upper()
	seqs.append(seq)
	return chains,seqs

def prepare_model(chains,seqs):
	n=1#line number
	temp_map={'A':A_temp,'U':U_temp,'C':C_temp,'G':G_temp,}
	number_map={'A':22,'U':20,'C':20,'G':23,}
	out=''
	for chain,seq in zip(chains,seqs):
		rsn=0# residue seq no
		for k,i in enumerate(seq):
			rsn+=1
			xx=[]
			for j in range(0,number_map.get(i,0)):
				xx.append(n)
				xx.append(chain)
				xx.append(rsn)
				n+=1
			temp=temp_map.get(i,None)
			if temp != None:
				out+=temp%tuple(xx)
		out+='TER   %5d        %c %c%4d                      \n'%(n,i,chain,rsn)
		n+=1
	return out

def format_pdb(fp,num=5):
	chains,seqs=readfasta(fp)
	out=''
	for i in xrange(num):
		out+='MODEL       %2d                                              \n'%(i+1)
		out+=prepare_model(chains,seqs)
		out+='ENDMDL                                                      \n'
	out+='END                                                        \n'
	print out,

if __name__ == '__main__':
	if( len(sys.argv) < 2 ):
		print Usage
		exit(0)
	elif (len(sys.argv) > 2):
		format_pdb(sys.argv[1],int(sys.argv[2]))
	else:
		format_pdb(sys.argv[1])

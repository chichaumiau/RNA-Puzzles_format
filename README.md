# RNA-Puzzzles_format
The submission format of RNA-Puzzles
___
rna_puzzles_format.py is a script used to generate a standard formatted PDB file for given sequences in a fasta file. 
Please refer to '2gdi.fa' for an example. 
To use:
*$python rna_puzzles_format.py [fasta file] [number of model(optional)] >out.pdb*
___
format_check.py is a script used to check the format of a submitted PDB file by referring to a standard PDB file of good format. If more than 1 error exist, a 'xx.format_check.txt' file is generated to include all the error reports. 
To use:
*$python format_check.py [To_be_checked.pdb] [Reference.pdb]*


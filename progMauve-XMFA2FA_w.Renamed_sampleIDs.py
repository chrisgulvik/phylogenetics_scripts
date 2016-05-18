#!/usr/bin/env python

# parses progressiveMauve XMFA outfile for sample (filename) IDs
# outputs a renamed FastA file
# Usage:  python script.py Input.xmfa Input.fa Output.fa

import os.path
import re
import sys

# make sure input file is from progressiveMauve
with open(sys.argv[1], 'rU') as XMFA:
	version = XMFA.readline().rstrip('\n').split(' ')[1]
	if version != 'Mauve1':
		print '\n\tERROR: input XMFA lacks Mauve1 header'
		print '\tis this an output file from progressiveMauve?\n'
		sys.exit()

# get filenames (in order)
	filename_list = []  #order is essential
	for line in XMFA:
		if re.match('\#Sequence[0-9]{0,}File', line):  #beginning of the string
			file_w_ext = os.path.basename(line.rstrip('\n'))
			filename = '.'.join(file_w_ext.split('.')[:-1])
			print '\t{}'.format(filename)
			filename_list.append(filename)
print '\n\tFound {} samples\n'.format(str(len(filename_list)))

# replace numerical sample IDs with filenames in FastA output
with open(sys.argv[2], 'rU') as FASTA:
	out_fasta = []
	for line in FASTA:
		if line.startswith('>'):
			sampleID_header = line.lstrip('>')
			try:
				int(sampleID_header)
				filename_header = (filename_list[(int(sampleID_header) - 1)])
				line = line.replace(sampleID_header, filename_header)
				out_fasta.append(line + '\n')
			except:  #some headers are NODE sequences -- not numerical sample IDs with filenames
				out_fasta.append(line)
		elif re.match(r'[ATCGatcg-]', line):  #linewrap is okay; allow gaps
			out_fasta.append(line)
		else:
			 print '\tERROR: unrecognized FastA format\n'
			 sys.exit()
	out_fasta_string = ''.join(out_fasta)  #convert list to string

# output reformatted file
with open (sys.argv[3], 'w') as fasta_renamed:
	fasta_renamed.write(out_fasta_string)

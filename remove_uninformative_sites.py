#!/usr/bin/env python


import argparse
import glob
import itertools
import os
import sys
from Bio import AlignIO
from Bio import SeqIO
from Bio.Align import AlignInfo
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import generic_dna
from Bio.Nexus import Nexus
from Bio.Seq import Seq

def parseArgs():
	parser = argparse.ArgumentParser(description='Converts an aligned FastA file into an aligned Nexus with only variant sites (especially useful for BEAST)')
	parser.add_argument('-i', '--infile', help='aligned FastA input file', required=True)
	parser.add_argument('-o', '--outfile', help='aligned variant-only Nexus output file [cwd-of-input/basename.informative.nex]')
	return parser.parse_args()

def get_variant_sites(infile):
	with open(infile, 'r') as input_handle:
		aln = AlignIO.read(input_handle, 'fasta')
		alignment = MultipleSeqAlignment([rec.upper() for rec in aln], annotations=aln.annotations)
		summary_align = AlignInfo.SummaryInfo(alignment)
		sequence = (alignment[0].seq)  #arbitrarily chose 1st seq to compare all others to

		pssm_sequence = summary_align.pos_specific_score_matrix(sequence)#, chars_to_ignore=['N'])
		#'pssm_sequence' is an obj containing line-by-line [site_identity, -, A, C, G, T]
		count = 0
		num_invariants = 0
		invariant_position_index = []

		for site in pssm_sequence.pssm:
			print site
			if site[0] != '-' and site[1]['-'] == 0:
				#site is a tuple of site_identity and dict of corresponding PSSM keys=nucleotides, vals=scores
				gap_score = site[1]['-']
				A_score   = site[1]['A']
				C_score   = site[1]['C']
				G_score   = site[1]['G']
				T_score   = site[1]['T']
				N_score   = site[1]['N']
				site_scores = [gap_score, A_score, C_score, G_score, T_score, N_score]  #same order as pssm obj 'pssm_sequence'
				binary_index = []
				for score in site_scores:
					if score > 0:
						binary_index.append(1)  #simplify with binary (present)
					else:
						binary_index.append(0)  #simplify with binary (absent)
				if sum(binary_index[1:len(binary_index)]) > 1:
					pass
				elif sum(binary_index[1:len(binary_index)]) == 1:
					num_invariants += 1
					invariant_position_index.append(count)
				else:
					sys.exit('ERROR_1')
				count += 1
			elif site[0] == '-' or site[1]['-'] > 0:
				num_invariants += 1
				invariant_position_index.append(count)
				count += 1
			else:
				sys.exit('ERROR_2')

	variant_position_index = []
	for i in range(0, len(sequence)):
		if i not in invariant_position_index:
			variant_position_index.append(i)

	print 'Discarded %s invariant sites' % str(len(invariant_position_index))
	return variant_position_index, alignment

def ranges(i):
	for (a, b) in itertools.groupby(enumerate(i), lambda (x, y): y - x):
		b = list(b)
		yield b[0][1], b[-1][1]

def write_variant_sites(alignment, var_sites, outfile):
	nex_aligns = []  #Bio.Nexus.Nexus.Nexus objects
	blocks = list(ranges(var_sites))  #tuples of positions
	for i in blocks:
		alignment_iteration = MultipleSeqAlignment(alignment[:, i[0]:i[1] + 1], alphabet=generic_dna).format('nexus')
		nex_aligns.append(('site %s to %s' % (str(i[0]), str(i[1] + 1)), Nexus.Nexus(alignment_iteration)))

	combined = Nexus.combine(nex_aligns)
	with open(outfile, 'w') as out:
		combined.write_nexus_data(out)
	print 'Converted %s informative sites into nexus alignment' % str(len(blocks))

def main():
	opts = parseArgs()
	infile = opts.infile
	outfile = opts.outfile

	if outfile is None:
		base = os.path.splitext(os.path.basename(infile))[0]
		outfile = base + '.informative.nex'

	var_sites, alignment = get_variant_sites(infile)
	write_variant_sites(alignment, var_sites, outfile)

if __name__ == '__main__':
	main()

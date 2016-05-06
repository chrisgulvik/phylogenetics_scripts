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
	parser = argparse.ArgumentParser(description='Converts an aligned multi-FastA file into a Nexus file with only variant sites.\
		Gap sites, where no reads mapped, are removed (not SNPs), but\
		N sites, where reads mapped but failed a quality filter, are maintained as ambiguous.')
	parser.add_argument('-i', '--infile', help='aligned FastA input file', required=True)
	parser.add_argument('-o', '--outfile', help='aligned variant-only Nexus output file [cwd-of-infile/basename.informative.nex]')
	return parser.parse_args()

def get_variant_sites(infile):
	with open(infile, 'r') as input_handle:
		aln = AlignIO.read(input_handle, 'fasta')
		alignment = MultipleSeqAlignment([rec.upper() for rec in aln],
			annotations=aln.annotations)
		summary_align = AlignInfo.SummaryInfo(alignment)
		sequence = alignment[0].seq  #arbitrarily chose 1st seq to compare all others to

		pssm_sequence = summary_align.pos_specific_score_matrix(sequence)#, chars_to_ignore=['N'])
		#'pssm_sequence' is an obj containing line-by-line [site_identity, -, A, C, G, T]
		site_index = 0
		num_invariants = 0
		variant_position_index = []

		num_records = len(alignment)
		print 'Number of records in multifasta: {}'.format(num_records)
		for site in pssm_sequence.pssm:
			sys.stdout.write('.')
			if site[0] != '-' and site[1]['-'] == 0:
			#site is a tuple of site_identity and dict of corresponding PSSM keys=nucleotides, vals=scores
				site_pssm_freqs = []
				for nucleotide in site[1]:
					site_pssm_freqs.append(int(str(site[1][nucleotide]).split('.')[0]))
				
				if num_records in site_pssm_freqs:
					num_invariants += 1
				elif num_records not in site_pssm_freqs:
					variant_position_index.append(site_index)
				else:
					sys.exit('ERROR_1')
				site_index += 1

			elif site[0] == '-' or site[1]['-'] > 0:  #skip sites with any number of gaps
				num_invariants += 1
				site_index += 1
			else:
				sys.exit('ERROR_2')

	print '\nFound {} variant sites without gaps'.format(str(len(variant_position_index)))
	print 'Discarded {} invariant and gap-containing sites'.format(
		str(num_invariants))
	return variant_position_index, alignment

def ranges(i):
	for (a, b) in itertools.groupby(enumerate(i), lambda (x, y): y-x):
		b = list(b)
		yield b[0][1], b[-1][1]

def write_variant_sites(alignment, var_sites, outfile):
	nex_aligns = []  #Bio.Nexus.Nexus.Nexus objects
	blocks = list(ranges(var_sites))  #tuples of positions
	for i in blocks:
		alignment_iteration = MultipleSeqAlignment(alignment[:, i[0]:i[1]+1],
			alphabet=generic_dna).format('nexus')
		# if i[0] == i[1]:
		# 	nex_aligns.append(('site {}'.format(str(i[1] + 1)),
		# 		Nexus.Nexus(alignment_iteration)))
		# else:
		# 	nex_aligns.append(('site {} to {}'.format(str(i[0]), str(i[1] + 1)),
		# 		Nexus.Nexus(alignment_iteration)))
		nex_aligns.append(('site {} to {}'.format(str(i[0]), str(i[1]+1)),
			Nexus.Nexus(alignment_iteration)))

	combined = Nexus.combine(nex_aligns)
	with open(outfile, 'w') as out:
		combined.write_nexus_data(out)
	print 'Converted {} informative sites without gaps into nexus alignment'.format(str(len(blocks)))

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

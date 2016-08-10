#!/usr/bin/env python


import os
import argparse
from Bio import AlignIO, SeqIO

def parseArgs():
	parser = argparse.ArgumentParser(usage='%(prog)s [-h] [-o FILE] [-l \'STR\'] FILE FILE [FILE ...]',
		description='Concatenates aligned multi-FastA files into a single multi-FastA. Sequence identifiers are extracted from the shortest prefix split by a period in each input filename. Sample names are parsed from the first input file\'s deflines, where the prefix is extracted from a split on whitespace and again the prefix is taken from a split at underscore. For example, seqid=<seqid>.all.fa.aln and sampleid=><sampleid>_00165 Ribosome-binding ATPase YchF')
	parser.add_argument('file', metavar='FILE', nargs=1, type=str,
		help=argparse.SUPPRESS)
	parser.add_argument('files', metavar='FILE', nargs='+', type=str,
		help='input FastA files')
	parser.add_argument('-l', '--linkseq', metavar='\'STR\'', default='',
		help='sequence to add between each sequence during concatenation, '
		'optionally flanked by single quotes  [None]')
	parser.add_argument('-o', '--outfile', metavar='FILE',
		default='concat_merged.mfasta', help='output file '
		'[./concat_merged.mfasta]')
	return parser.parse_args()

def main():
	opts = parseArgs()
	infiles = opts.files + opts.file
	l = opts.linkseq.lstrip('\'').rstrip('\'')

	# parse sequence IDs from all input filenames
	seq_ids = []
	for file in infiles:
		seq_ids.append(os.path.basename(os.path.abspath(file)).split('.')[0])
	print '{} sequence identifiers:\n{}\n'.format(len(seq_ids),
												  ' '.join(seq_ids))

	# parse sample IDs from the first input file's deflines
	sample_ids = []
	f = SeqIO.parse(opts.file[0], 'fasta')
	for rec in f:
		sample_ids.append(rec.id.split('_')[0])  #handle locus tags as IDs
	print '{} sample identifiers:\n{}\n'.format(len(sample_ids),
												' '.join(sample_ids))

	# concatenate sequences with the same sample ID
	merged = {}
	for file in infiles:
		d = {}
		for record in AlignIO.parse(file, 'fasta'):
			for r in record:
				d[str(r.id.split('_')[0])] = str(r.seq)
		for sample_id in sample_ids:
			if sample_id in merged:
				merged[sample_id] = str(d[sample_id]) + l + merged[sample_id]
			else:
				merged[sample_id] = str(d[sample_id])

	# write the file
	with open(opts.outfile, 'w') as o:
		for k, v in merged.iteritems():
			o.write('>'+k+'\n'+v+'\n')

if __name__ == '__main__':
	main()

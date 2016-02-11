#!/usr/bin/env python
from __future__ import division

import argparse
import fileinput
import getpass
import glob
import logging
import os
import pwd
import re
import socket
import subprocess
import sys
import tempfile
import threading
from Bio import AlignIO
from decimal import Decimal, ROUND_UP


def parseArgs():
	parser = argparse.ArgumentParser(
		description='creates a Maximum Likelihood tree from an aligned FastA file')
	parser.add_argument('-f', '--fasta', required=True, help='input aligned FastA file')
	parser.add_argument('-e', '--rmext', required=False, help='extension of input FastA file to remove, e.g., -e .fas.aln [.*]')
	parser.add_argument('-r', '--raxml', required=False, default='raxmlHPC-PTHREADS-AVX', help='available raxml binary with PTHREADS [raxmlHPC-PTHREADS-AVX]')
	parser.add_argument('-o', '--outpath', required=False, default='ml_out', help='output folder [`cwd`/ml_out]')
	parser.add_argument('-b', '--boots', required=False, type=int, default='250', help='number of bootstraps [250]')
	parser.add_argument('-T', '--cpus', required=False, type=int, default='1', help='number of CPUs [1]')
	parser.add_argument('-k', '--keep', required=False, default=None, action='store_true', help='keep all intermediate files')
	args = parser.parse_args()
	return args

class LogPipe(threading.Thread):
	def __init__(self, level):
		''' sets up the object with a logger and starts up the thread '''
		threading.Thread.__init__(self)
		self.daemon = False
		self.level = level
		self.fd_read, self.fd_write = os.pipe()
		self.pipe_reader = os.fdopen(self.fd_read)
		self.start()

	def fileno(self):
		''' returns the write file descriptor of the pipe '''
		return self.fd_write

	def run(self):
		''' runs the thread and logs everything '''
		for line in iter(self.pipe_reader.readline, ''):
			# logging.debug(line.strip('\n'))
			logging.log(self.level, line.strip('\n'))
		self.pipe_reader.close()

	def close(self):
		''' closes the write end of the pipe '''
		os.close(self.fd_write)


def fa2phy(fasta, outphy):
	''' converts FastA format into phylip format '''
	with open(fasta) as handle:
		records = AlignIO.parse(handle, 'fasta')
		with open(outphy, 'w') as out:
			AlignIO.write(records, out, 'phylip-sequential')

def dependency(dep):
	''' checks for binary availability '''
	for path in os.environ.get('PATH', '').split(':'):
		if os.path.exists(os.path.join(path, dep)) and \
		not os.path.isdir(os.path.join(path, dep)):
			return os.path.join(path, dep)
	return None

def syscall(syscmd, outpath):
	''' executes a system call and optionally logs stdout and stderr '''
	if outpath is 'dump':
		with open(os.devnull) as dump:
			returncode = subprocess.call(syscmd.split(), stdout=dump, stderr=dump, shell=False)
			if returncode != 0:
				logging.error('failed system call ' + syscmd)
				sys.exit('ERROR: failed system call\ncheck the logfile')
	else:
		o = LogPipe(logging.DEBUG)
		e = LogPipe(logging.ERROR)
		s = subprocess.Popen(syscmd.split(), stdout=o, stderr=e)
		if s.wait() != 0:  #return code
			logging.error('failed system call ' + syscmd)
			o.close()
			e.close()
			sys.exit('ERROR: failed system call\ncheck the logfile')	
		o.close()
		e.close()

def run_PhyML(phylip, base, boots, outpath, freqs, kappa, alpha, search, xopts):
	''' runs PhyML '''
	cmd = 'phyml -i %s --datatype nt --bootstrap %s -f %s --ts/tv %s --alpha %s --search %s %s' % (phylip, boots, freqs, kappa, alpha, search, xopts)
	syscall(cmd, outpath)

def fill_empty_var(var):
	''' sets variable to e (estimate in phyml) if empty '''
	if var is None:
		logging.notice('unable to parse %s; PhyML will estimate this in subsequent runs' % var)
		return 'e'
	else:
		return var

def get_PhyML_stats(value=None):
	''' parses PhyML stats '''
	if os.path.isfile('init_phyml_stats.txt') and os.path.getsize('init_phyml_stats.txt') > 0:
		for phymlstats in open('init_phyml_stats.txt', 'r'):
			if phymlstats.startswith('. Log-likelihood:'):
				loglike = phymlstats.strip().split()[-1]
				logging.info('calculated %s log-likelihood from PhyML statistics' % loglike)
			elif phymlstats.startswith('  - Gamma shape parameter:'):
				alpha = phymlstats.strip().split()[-1]
				logging.info('calculated %s gamma shape from PhyML statistics' % alpha)
			elif phymlstats.startswith('  - f(A)='):
				fA = phymlstats.strip().split()[-1]
			elif phymlstats.startswith('  - f(C)='):
				fC = phymlstats.strip().split()[-1]
			elif phymlstats.startswith('  - f(G)='):
				fG = phymlstats.strip().split()[-1]
			elif phymlstats.startswith('  - f(T)='):
				fT = phymlstats.strip().split()[-1]
				freqs = '%s,%s,%s,%s' % (fA, fC, fG, fT)  # Order (A,C,G,T) matters for PhyML
				logging.info('calculated %s frequencies of A,C,G,T from PhyML statistics' % freqs)
			elif phymlstats.startswith('  A <-> C'):
				A2C = phymlstats.strip().split()[-1]
			elif phymlstats.startswith('  A <-> G'):
				A2G = phymlstats.strip().split()[-1]
			elif phymlstats.startswith('  A <-> T'):
				A2T = phymlstats.strip().split()[-1]
			elif phymlstats.startswith('  C <-> G'):
				C2G = phymlstats.strip().split()[-1]
			elif phymlstats.startswith('  C <-> T'):
				C2T = phymlstats.strip().split()[-1]
			elif phymlstats.startswith('  G <-> T'):
				G2T = phymlstats.strip().split()[-1]
				Transitions = (Decimal(A2G) + Decimal(C2T)) # A<->G && C<->T
				logging.info('calculated %s transitions from PhyML statistics' % Transitions)
				Transversions = (Decimal(A2C) + Decimal(A2T) + Decimal(C2G) + Decimal(G2T))  # All (4) others
				logging.info('calculated %s transversions from PhyML statistics' % Transversions)
				TsTv = (Transitions / Transversions).quantize(Decimal('.000001'), rounding=ROUND_UP)
				logging.info('calculated %s Ts:Tv from PhyML statistics' % TsTv)
		for var in [alpha, freqs, TsTv]:
			var = fill_empty_var(var)
	else:
		logging.notice('unable to open phyml stats file to optimize parameters')
		logging.notice('gamma shape, Ts:Tv, and f(A),f(T),f(C),f(G) values will be estimated in subsequent runs')
		alpha = 'e'
		TsTv = 'e'
		freqs = 'm'
	return alpha, TsTv, freqs

def run_RAxML(raxml, cpus, xopts, outpath, base):
	''' runs RAxML '''
	cmd = '%s -T %s %s' % (raxml, cpus, xopts)
	syscall(cmd, outpath)

def main():
	opts = parseArgs()
	cpus = opts.cpus
	fasta = opts.fasta
	outpath = opts.outpath
	raxml = opts.raxml

	# Setup logging
	if not os.path.exists(outpath): 
		os.mkdir(outpath)
	logging.basicConfig(filename=os.path.join(outpath, 'make_ML_tree.log'),
		format='%(asctime)s: %(levelname)s: %(message)s', level=logging.DEBUG,
		filemode='w', datefmt='%a %d-%b-%Y %H:%M:%S')
	logging.info('user: %s' % pwd.getpwuid(os.getuid()).pw_name)
	logging.info('release: %s' % os.uname()[3])
	logging.info('shell env: %s' % pwd.getpwuid(os.getuid()).pw_shell)
	logging.info('cwd: %s' % pwd.getpwuid(os.getuid()).pw_dir)
	logging.info('python version: %s' % sys.version)

	# Check dependencies
	RAxMLbin = dependency(raxml)
	PhyMLbin = dependency('phyml')
	if RAxMLbin:
		logging.info('found %s' % RAxMLbin)
		logging.info(subprocess.check_output('%s -v | grep "RAxML version"' % raxml, shell=True).rstrip())
	else:
		print '\t%s not found' % raxml
		sys.exit(1)
	if PhyMLbin:
		logging.info('found %s' % PhyMLbin)
		logging.info(subprocess.check_output('phyml --version | grep "PhyML version"', shell=True).rstrip())
	else:
		print '\tphyml not found'
		sys.exit(1)

	# I/O handling
	path = os.path.dirname(os.path.abspath(fasta))  #handles relational input
	if opts.rmext:
		base = os.path.basename(fasta).rstrip(opts.rmext)
	else:
		base = os.path.basename(fasta).rsplit('.', 1)[0]

	# Make phylip file
	fa2phy(fasta, base + '.phy')
	phylip = os.path.join(path, base + '.phy')

	# Calculate parameters and make initial tree
	run_PhyML(phylip, base, '0', outpath, 'm', 'e', 'e',
		'SPR', '-m GTR -o lr --quiet --sequential')
	if os.path.isfile(base + '.phy_phyml_stats.txt') and os.path.getsize(base + '.phy_phyml_stats.txt') > 0:
		os.rename(base + '.phy_phyml_stats.txt', 'init_phyml_stats.txt')
		os.rename(base + '.phy_phyml_tree.txt', 'init_tree.nwk')
	elif os.path.isfile(base + '.phy_phyml_stats') and os.path.getsize(base + '.phy_phyml_stats') > 0:
		os.rename(base + '.phy_phyml_stats', 'init_phyml_stats.txt')
		os.rename(base + '.phy_phyml_tree', 'init_tree.nwk')
	if not os.path.isfile('init_tree.nwk'):
		logging.error('init_tree.nwk absent')
		sys.exit('ERROR: init_tree.nwk absent')
	(alpha, kappa, freqs) = get_PhyML_stats()

	# Find ML topology
	find_ML_topology = '-s %s -n topology -f d -j -m GTRGAMMA -t init_tree.nwk' % phylip
	run_RAxML(raxml, cpus, find_ML_topology, outpath, base)
	
	# Optimize branch lengths of ML topology
	run_PhyML(phylip, base, '0', outpath, freqs, kappa, alpha,
		'BEST', '-u RAxML_result.topology -m GTR -o lr --quiet --sequential')
	# some PhyML versions (e.g., 20160116) append file extension to output files (.txt)
	try:
		os.rename(base + '.phy_phyml_tree', base + '.ML_topol.nwk')
	# some PhyML versions (e.g., 20131022) lack a file extension to output files (.txt)
	except OSError as error:
		if error.errno == 2:
			os.rename(base + '.phy_phyml_tree.txt', base + '.ML_topol.nwk')
		else:
			raise

	# Generate 250 boostrap trees
	generate_bootstrap_trees = '-n boots -s %s -N %s -f a -m GTRGAMMA -p 65432 -x 54321' % (phylip, opts.boots)
	run_RAxML(raxml, cpus, generate_bootstrap_trees, outpath, base)
	
	# Project boostrap values onto ML topology
	project_boots_onto_ML_topology = '-n final -t %s.ML_topol.nwk -z RAxML_bootstrap.boots -f b -m GTRGAMMA' % base
	run_RAxML(raxml, cpus, project_boots_onto_ML_topology, outpath, base)
	os.rename('RAxML_bipartitions.final', os.path.join(outpath, base + '.bootstrapped_ML.nwk'))
	
	# Cleanup the mess I made by default
	if opts.keep is None:
		os.remove(base + '.phy')
		os.remove('RAxML_flagCheck')
		try:
			os.remove(base + '.phy_phyml_stats.txt')
		except OSError as error:
			if error.errno == 2:
				os.remove(base + '.phy_phyml_stats')
			else:
				raise
		for s in ['tree.nwk', 'phyml_stats.txt']:
			os.remove('init_' + s)
		# RAxML topology files
		for junk in glob.glob('RAxML_checkpoint.topology.*'):
			os.remove(junk)
		for s in ['bestTree', 'info', 'log', 'result']:
			os.remove('RAxML_' + s + '.topology')
		os.remove(base + '.ML_topol.nwk')
		# RAxML boots files
		for s in ['bestTree', 'bipartitions', 'bipartitionsBranchLabels', 'bootstrap', 'info']:
			os.remove('RAxML_' + s + '.boots')
		# RAxML final files
		os.remove('RAxML_info' + '.final')
		os.remove('RAxML_bipartitionsBranchLabels' + '.final')
		logging.info('removed all intermediate files')
	logging.info('completed')
	sys.exit()

if __name__ == '__main__':
	main()

#!/usr/bin/env python
# -*- coding: utf-8 -*- 
import os
import sys
import argparse
from intervaltree import Interval, IntervalTree

def parse_snp(snp_file):
	snp = {}
	for line in open(snp_file):
		chr_name = line.split()[0]
		if not snp.has_key(chr_name):
			snp[chr_name] = IntervalTree()
		# 0-based
		beg = int(line.split()[1])
		end = int(line.split()[2])
		if beg == end:
			#insertion
			end += 0.1
		# not including upper limit
		snp[chr_name].add(Interval(beg, end))
#		snp[chr_name][beg:end] = line.split()[3]
	return snp

def process_tl(snp, tl_dir, out_dir, annotation):
	tl_dir_list = os.listdir(tl_dir)
	for tl_file in tl_dir_list:
		cur_chr = 'chr' + tl_file.split('_')[0]
		if cur_chr == 'chr23':
			cur_chr = 'chrX'
		elif cur_chr == 'chr24':
			cur_chr = 'chrY'
		if not cur_chr in snp.keys():
			continue
		out_file = os.path.join(out_dir, tl_file.split('_')[0] + '_tl_' + annotation + '_annotation.txt')
		fout = open(out_file, 'w')		
		for line in open(os.path.join(tl_dir, tl_file)):
			# 0-based
			# C coordinate
			pos = int(line.split()[1])
			strand = line.split()[2]
			if strand == '+':
				if len(snp[cur_chr][pos]) > 0 or len(snp[cur_chr][pos+1]) > 0:
					fout.write(str(pos) + '\t+\n')
#					fout.write(str(pos) + '\t')
#					for elem in snp[cur_chr][pos]:
#						fout.write(elem.data + '\t')
#					for elem in snp[cur_chr][pos+1]:
#						fout.write(elem.data + '\t')
#					fout.write('\n')
				else:
					fout.write(str(pos) + '\t-\n')
			elif strand == '-':
				if len(snp[cur_chr][pos]) > 0 or len(snp[cur_chr][pos-1]) > 0:
					fout.write(str(pos) + '\t+\n')
#					fout.write(str(pos) + '\t')
#					for elem in snp[cur_chr][pos]:
#						fout.write(elem.data + '\t')
#					for elem in snp[cur_chr][pos-1]:
#						fout.write(elem.data + '\t')
#					fout.write('\n')
				else:
					fout.write(str(pos) + '\t-\n')
			else:
				print 'unknown strad', strand
				fout.write(str(pos) + '\n')
		fout.close()

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print 'Usage:', sys.argv[0], '-s <snp annotation file> -t <traffic lights directory> -a <annotation> -o <output directory>'
		exit()

	parser = argparse.ArgumentParser(prog = sys.argv[0], description='GC and CpG content count')
	parser.add_argument('-s', '--snp_file', help='snp annotation file directory', required=True)
	parser.add_argument('-t', '--tl_dir', help='traffic lights directory', required=True)
	parser.add_argument('-o', '--out_dir', help='output directory', required=True)
	parser.add_argument('-a', '--ann', help='annotation', required=True)
	args = parser.parse_args()
	snp_file = args.snp_file
	tl_dir = args.tl_dir
	out_dir = args.out_dir
	annotation = args.ann

	if not os.path.isfile(snp_file):
		print >> sys.stderr, 'Not a file ' + snp_file
		sys.exit(1)

	if not os.path.isdir(tl_dir):
		print >> sys.stderr, 'Not a directory ' + tl_dir
		sys.exit(1)

	if not os.path.exists(out_dir):
		os.mkdir(out_dir)

	snp = parse_snp(snp_file)
	process_tl(snp, tl_dir, out_dir, annotation)
	

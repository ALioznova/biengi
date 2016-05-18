#!/usr/bin/env python
# -*- coding: utf-8 -*- 
import os
import sys
import argparse
from intervaltree import Interval, IntervalTree

def parse_enhancer(enhancer_file):
	enhancer = {}
	fin = open(enhancer_file)
	fin.readline()
	for line in fin:
		coords_info = line.split()[0]
		usage = sum([int(elem) for elem in line.split()[1:]])
		chr_name = coords_info.split(':')[0]
		if not enhancer.has_key(chr_name):
			enhancer[chr_name] = IntervalTree()
		# 0-based
		beg = int(coords_info.split(':')[1].split('-')[0])
		end = int(coords_info.split(':')[1].split('-')[1])
		enhancer[chr_name][beg:end] = str(usage)
	fin.close()
	return enhancer

def process_tl(enhancer, tl_dir, out_dir, out_file_name_ending):
	tl_dir_list = os.listdir(tl_dir)
	for tl_file in tl_dir_list:
		cur_chr = 'chr' + tl_file.split('_')[0]
		if cur_chr == 'chr23':
			cur_chr = 'chrX'
		elif cur_chr == 'chr24':
			cur_chr = 'chrY'
		out_file = os.path.join(out_dir, tl_file.split('_')[0] + out_file_name_ending)
		fout = open(out_file, 'w')		
		for line in open(os.path.join(tl_dir, tl_file)):
			# 0-based!
			# C coordinate
			pos = int(line.split()[1])
			strand = line.split()[2]
			if strand == '+':
				if not enhancer.has_key(cur_chr):
					fout.write(str(pos) + '\t-\n')
					continue
				if len(enhancer[cur_chr][pos]) > 0 or len(enhancer[cur_chr][pos+1]) > 0:
					fout.write(str(pos) + '\t')
					for elem in enhancer[cur_chr][pos]:
						fout.write(elem.data + '\t')
					for elem in enhancer[cur_chr][pos+1]:
						fout.write(elem.data + '\t')
					fout.write('\n')
				else:
					fout.write(str(pos) + '\t-\n')
			elif strand == '-':
				if not enhancer.has_key(cur_chr):
					fout.write(str(pos) + '\t-\n')
					continue
				if len(enhancer[cur_chr][pos]) > 0 or len(enhancer[cur_chr][pos-1]) > 0:
					fout.write(str(pos) + '\t')
					for elem in enhancer[cur_chr][pos]:
						fout.write(elem.data + '\t')
					for elem in enhancer[cur_chr][pos-1]:
						fout.write(elem.data + '\t')
					fout.write('\n')
				else:
					fout.write(str(pos) + '\t-\n')
			else:
				print 'unknown strad', strand
				fout.write(str(pos) + '\n')
		fout.close()

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print 'Usage:', sys.argv[0], '-e <enhancer annotation directory> -t <traffic lights directory> -o <output directory>'
		exit()

	parser = argparse.ArgumentParser(prog = sys.argv[0], description='GC and CpG content count')
	parser.add_argument('-e', '--enhancer_dir', help='enhancer annotation directory', required=True)
	parser.add_argument('-t', '--tl_dir', help='traffic lights directory', required=True)
	parser.add_argument('-o', '--out_dir', help='output directory', required=True)
	args = parser.parse_args()
	enhancer_dir = args.enhancer_dir
	tl_dir = args.tl_dir
	out_dir = args.out_dir

	if not os.path.isdir(enhancer_dir):
		print >> sys.stderr, 'Not a dir ' + enhancer_dir
		sys.exit(1)

	if not os.path.isdir(tl_dir):
		print >> sys.stderr, 'Not a directory ' + tl_dir
		sys.exit(1)

	if not os.path.exists(out_dir):
		os.mkdir(out_dir)

	for enhancer_file in os.listdir(enhancer_dir):
		enhancer = parse_enhancer(os.path.join(enhancer_dir, enhancer_file))
		out_file_name_ending = '_tl_' + enhancer_file[len('F5.hg38.enhancers_'):-len('.txt')] + '_annotation.txt'
		process_tl(enhancer, tl_dir, out_dir, out_file_name_ending)


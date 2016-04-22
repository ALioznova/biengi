#!/usr/bin/env python
# -*- coding: utf-8 -*- 
import os
import sys
import argparse
import bx.bbi.bigwig_file

def process_tl(ann_file, tl_dir, out_dir):
	bwh = bx.bbi.bigwig_file.BigWigFile(open(ann_file, "rb"))
	window = 5
	tl_dir_list = os.listdir(tl_dir)
	for tl_file in tl_dir_list:
		cur_chr = 'chr' + tl_file.split('_')[0]
		if cur_chr == 'chr23':
			cur_chr = 'chrX'
		elif cur_chr == 'chr24':
			cur_chr = 'chrY'
		out_file = os.path.join(out_dir, tl_file.split('_')[0] + '_tl_phyloP7way_annotation.txt')
		fout = open(out_file, 'w')		
		for line in open(os.path.join(tl_dir, tl_file)):
			# 0-based
			# C coordinate
			pos = int(line.split()[1])
			strand = line.split()[2]
			if strand == '+':
				data = bwh.get_as_array(cur_chr, pos-window, pos+1+window)
				if data != None:
					fout.write(str(pos))
					for elem in data:
						fout.write('\t' + str(elem))
					fout.write('\n')
				else:
					fout.write(str(pos) + '\t-\n')
			elif strand == '-':
				data = bwh.get_as_array(cur_chr, pos-1-window, pos+window)
				if data != None:
					fout.write(str(pos))
					for elem in data:
						fout.write('\t' + str(elem))
					fout.write('\n')
				else:
					fout.write(str(pos) + '\t-\n')
			else:
				print 'unknown strad', strand
				fout.write(str(pos) + '\n')
		fout.close()

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print 'Usage:', sys.argv[0], '-a <bigwig annotation file> -t <traffic lights directory> -o <output directory>'
		exit()
	
	parser = argparse.ArgumentParser(prog = sys.argv[0], description='Annotation from BigWig')
	parser.add_argument('-a', '--annotation', help='annotation file', required=True)
	parser.add_argument('-t', '--tl_dir', help='traffic lights directory', required=True)
	parser.add_argument('-o', '--out_dir', help='output directory', required=True)
	args = parser.parse_args()
	ann_file = args.annotation
	tl_dir = args.tl_dir
	out_dir = args.out_dir

	if not os.path.isfile(ann_file):
		print >> sys.stderr, 'Not a file ' + ann_file
		sys.exit(1)

	if not os.path.isdir(tl_dir):
		print >> sys.stderr, 'Not a directory ' + tl_dir
		sys.exit(1)

	if not os.path.exists(out_dir):
		os.mkdir(out_dir)

	process_tl(ann_file, tl_dir, out_dir)


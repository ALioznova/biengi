#!/usr/bin/env python
# -*- coding: utf-8 -*- 
import os
import sys
import argparse
from sets import Set

class tl_record:
	def __init__(self,  tl_line, tl_line_num, tl_id):
		self.tl_id = tl_id
		self.line_num = tl_line_num # 1-based
		self.gene_info = tl_line.split()[0]
		self.gene_name = self.gene_info.split(';')[0]
		self.pos = int(tl_line.split()[1])
		self.strand = tl_line.split()[2]
		self.annotation = {}
		for elem in tl_line.split()[3].split(','):
			if len(elem.split('=')) == 1:
				self.annotation[elem] = True
			else:
				self.annotation[elem.split('=')[0]] = float(elem.split('=')[1])
		calculated_data = [elem.split(')')[0].split() for elem in tl_line.split('(')[1:]]
		calculated_data = [[int(elem[0]), float(elem[1]), float(elem[2]), float(elem[3]), float(elem[4])] for elem in calculated_data]
		self.extra = calculated_data
		selected_data = calculated_data[0]
		self.num = selected_data[0]
		self.corr = selected_data[1]
		self.p_corr = selected_data[2]
		self.p_corr_fdr = selected_data[3]
		self.cause = selected_data[4]

def get_records(tl_dir):
	tl_records_scope = {}
	tl_dir_list = os.listdir(tl_dir)
	for tl_name in tl_dir_list:
		print 'processing', os.path.basename(tl_name)
		tl_id = len(tl_records_scope)
		tl_line_num = 1
		for line in open(os.path.join(tl_dir, tl_name)):
			tl_records_scope[tl_id] = tl_record(line, tl_line_num, tl_id)
			tl_id += 1
			tl_line_num += 1
	return tl_records_scope

def split_for_background(tl_records_scope):
	p_fdr_threshold = 0.2
	tl_rec = {}
	bg_rec = {}
	for (key, tl) in tl_records_scope.iteritems():
		if tl.p_corr_fdr < p_fdr_threshold:
			tl_rec[key] = tl
		else:
			bg_rec[key] = tl
	return (bg_rec, tl_rec)

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print 'Usage:', sys.argv[0], '-t <traffic lights directory> -o <output directory>'
		exit()
	
	parser = argparse.ArgumentParser(prog = sys.argv[0], description='Retrieve TL and BG genes')
	parser.add_argument('-t', '--tl_dir', help='traffic lights directory', required=True)
	parser.add_argument('-o', '--out_dir', help='output directory', required=True)
	args = parser.parse_args()
	tl_dir = args.tl_dir
	out_dir = args.out_dir

	if not os.path.isdir(tl_dir):
		print >> sys.stderr, 'Not a directory ' + tl_dir
		sys.exit(1)

	if not os.path.exists(out_dir):
		os.mkdir(out_dir)

	tl_records_scope = get_records(tl_dir)
	(bg, tl) = split_for_background(tl_records_scope)

	print len(bg), len(tl)
	
	bg_genes = Set()
	for elem in bg.itervalues():
		bg_genes.add(elem.gene_name)
	tl_genes = Set()
	for elem in tl.itervalues():
		tl_genes.add(elem.gene_name)

	fout_bg = open(os.path.join(out_dir, 'bg_genes.txt'), 'w')
	for elem in bg_genes:
		fout_bg.write(elem + '\n')
	fout_bg.close()

	fout_tl = open(os.path.join(out_dir, 'tl_genes.txt'), 'w')
	for elem in tl_genes:
		fout_tl.write(elem + '\n')
	fout_tl.close()


#!/usr/bin/env python
# -*- coding: utf-8 -*- 
import os
import sys
import argparse
import numpy

class Pos_rec:
	def __init__(self, description):
		(chr_name, begin1, end1, begin2, end2) = description.split(':')
		self.chr_name = chr_name[:-1]
		self.strand = chr_name[-1]
		self.beg1 = int(begin1)
		self.end1 = int(end1)
		self.beg2 = int(begin2)
		self.end2 = int(end2)

def read_data(in_fn):
	in_f = open(in_fn)
	header_line = in_f.readline()
	tumor_setd2_broken_num = int((header_line.split()[1]).split(':')[1])
	tumor_num = int((header_line.split()[2]).split(':')[1])
	normal_num = int((header_line.split()[3]).split(':')[1])
	pos = []
	tumor_setd2_broken = []
	tumor = []
	normal = []
	for line in in_f:
		pos.append(Pos_rec(line.split()[0]))
		tumor_setd2_broken.append(float(line.split()[1]))
		tumor.append(float(line.split()[2]))
		normal.append(float(line.split()[3]))
	in_f.close()
	return (pos, tumor_setd2_broken_num, tumor_setd2_broken, tumor_num, tumor, normal_num, normal)

class Genes:
	def __init__(self, name, begin, end, exons):
		self.name = name
		self.begin = begin
		self.end = end
		self.exons = exons

def parse_genes_file(genes_file):
	gf = open(genes_file)
	gf.readline()
	genes_data = []
	for line in gf:
		if line.split('\t')[2] == 'gene':
			name = line.split('\t')[1]
			exons = line.split('\t')[14]
			locus = line.split('\t')[16]
#			genes_data.append(Genes(name, begin, end, exons))
	return genes_data

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print 'Usage:', sys.argv[0], '-d <data directory> -g <gene annotation file>'
		exit()

	parser = argparse.ArgumentParser(prog = sys.argv[0], description='Analyze')
	parser.add_argument('-d', '--data_dir', help='data directory', required=True)
	parser.add_argument('-g', '--genes', help='gene annotation file', required=True)
	args = parser.parse_args()
	data_dir = args.data_dir
	genes_fn = args.genes
	
	if not os.path.isfile(genes_fn):
		print 'No such file', genes_fn
		sys.exit(1)
	
	parse_genes_file(genes_fn)

	data_dir_list = [os.path.join(data_dir, d) for d in os.listdir(data_dir)]
	for d in data_dir_list:
		print 'processing', os.path.basename(d)
		in_fn = os.path.join(d, os.path.basename(d) + '__t_setd2__t__n.txt')
		(pos, tumor_setd2_broken_num, tumor_setd2_broken, tumor_num, tumor, normal_num, normal) = read_data(in_fn)

		gene_expression_fn = None
		for elem in os.listdir(d):
			if 'gene_expression' in elem:
				gene_expression_fn = os.path.join(d, elem)
		if not (gene_expression_fn):
			continue
		
		
		

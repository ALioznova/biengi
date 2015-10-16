#!/usr/bin/env python
# -*- coding: utf-8 -*- 
import os
import sys
import argparse
import numpy
from sets import Set
from intervaltree import Interval, IntervalTree

class Pos_rec:
	def __init__(self, description):
		(chr_name, begin1, end1, begin2, end2, strand) = description.split(':')
		self.chr_name = chr_name
		self.beg1 = int(begin1)
		self.end1 = int(end1)
		self.beg2 = int(begin2)
		self.end2 = int(end2)
		if strand == '+':
			self.strand = True
		else:
			self.strand = False

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
	def __init__(self, name, locus, exons):
		self.name = name
		self.locus = locus
		self.exons = exons

class Locus:
	def __init__(self, description):
		(chr_name, coords, strand) = description.split(':')
		self.chr_name = chr_name
		self.beg = int(coords.split('-')[0])
		self.end = int(coords.split('-')[1])
		if strand == '+':
			self.strand = True
		else:
			self.strand = False
		
class Exons:
	def __init__(self, description):
		(chr_name, coords, strand) = description.split(':')
		self.chr_name = chr_name
		if strand == '+':
			self.strand = True
		else:
			self.strand = False
		self.begs = [int(c.split('-')[0]) for c in coords.split(',')]
		self.ends = [int(c.split('-')[1]) for c in coords.split(',')]

def parse_genes_file(genes_file):
	# data source: https://tcga-data.nci.nih.gov/docs/GAF/GAF.hg19.June2011.bundle/outputs/TCGA.hg19.June2011.gaf
	gf = open(genes_file)
	gf.readline()
	genes_data = {}
	exon_dict = {}
	interval_trees = {}
	for line in gf:
		if line.split('\t')[2] == 'gene':
			name = line.split('\t')[1]
			exons = Exons(line.split('\t')[14])
			locus = [Locus(l) for l in line.split('\t')[16].split(';')]
			# filling genes_data: gene_name -> [Genes]
			if len(name.split('|')) > 2: # gene has multiple locuses
				if not genes_data.has_key('|'.join(name.split('|')[:-1])):
					genes_data['|'.join(name.split('|')[:-1])] = []
				genes_data['|'.join(name.split('|')[:-1])].append(Genes(name, locus, exons))
			else:
				genes_data[name] = [Genes(name, locus, exons)] # the only line with gene description
			if not exon_dict.has_key((exons.chr_name, exons.strand)):
				exon_dict[(exons.chr_name, exons.strand)] = {}
			if not interval_trees.has_key((exons.chr_name, exons.strand)):
				interval_trees[(exons.chr_name, exons.strand)] = IntervalTree()
			for i in xrange(len(exons.begs)):
				# filling exon_dict: (chr, strand) -> (exon_start, exon_end) -> Set(gene_name)
				if exons.begs[i] == exons.ends[i]:
					continue
				if not exon_dict[(exons.chr_name, exons.strand)].has_key((exons.begs[i], exons.ends[i])):
					exon_dict[(exons.chr_name, exons.strand)][(exons.begs[i], exons.ends[i])] = Set()
				exon_dict[(exons.chr_name, exons.strand)][(exons.begs[i], exons.ends[i])].add('|'.join(name.split('|')[:2]))
				# filling interval_trees for annotated exons:
				interval_trees[(exons.chr_name, exons.strand)][exons.begs[i] : exons.ends[i]] = '|'.join(name.split('|')[:2])

	return (genes_data, exon_dict, interval_trees)

def exon_to_gene_names(exon_dict, interval_trees, pos):
	exon_names = {}
	for p in pos:
		if not exon_names.has_key((p.chr_name, p.strand)):
			exon_names[(p.chr_name, p.strand)] = {}
		if exon_dict[(p.chr_name, p.strand)].has_key((p.end1, p.beg2)):
			exon_names[(p.chr_name, p.strand)][(p.end1, p.beg2)] = exon_dict[(p.chr_name, p.strand)][(p.end1, p.beg2)]
		else:
			left_genes = Set([lg[2] for lg in interval_trees[(p.chr_name, p.strand)][p.end1 + 0.1]])
			right_genes = Set([rg[2] for rg in interval_trees[(p.chr_name, p.strand)][p.beg2 - 0.1]])
			if len(left_genes) != 0 or len(right_genes) != 0:
				if len(left_genes.intersection(right_genes)) != 0:
					exon_names[(p.chr_name, p.strand)][(p.end1, p.beg2)] = left_genes.intersection(right_genes)
	return exon_names

def get_dist_exon_gene_beg(exon_to_genes, genes_data, pos_rec):
	def get_overlap(a, b):
		return max(0, min(a[1], b[1]) - max(a[0], b[0]))

	exon_chr = pos_rec.chr_name
	exon_strand = pos_rec.strand
	exon_beg = pos_rec.end1
	exon_end = pos_rec.beg2
	if not exon_to_genes[(exon_chr, exon_strand)].has_key((exon_beg, exon_end)):
		return []
	exon_genes = list(exon_to_genes[(exon_chr, exon_strand)][(exon_beg, exon_end)])
	dist= []
	for exon_gene in exon_genes:
		for gene in genes_data[exon_gene]:
			if gene.exons.chr_name != exon_chr:
				continue
			for i in xrange(len(gene.exons.begs)):
				# if there is no intersection with exon_beg, exon_end => it is other exon's locus
				ge_beg = gene.exons.begs[i]
				ge_end = gene.exons.ends[i]
				if get_overlap((ge_beg, ge_end), (exon_beg, exon_end)) != 0:
					for loc in gene.locus:
						if loc.chr_name != gene.exons.chr_name:
							continue
						if get_overlap((loc.beg, loc.end), (ge_beg, ge_end)) != 0:
							if loc.strand == exon_strand:
								if exon_beg - loc.beg < 0 or exon_beg - loc.beg > loc.end - loc.beg:
									continue
								dist.append(exon_beg - loc.beg)
							else:
								if loc.end - exon_end < 0 or loc.end - exon_end > loc.end - loc.beg:
									continue
								dist.append(loc.end - exon_end)
	return dist

def get_gene_dist(exon_to_genes, genes_data, pos):
	gene_beg_dist = []
	for i in xrange(len(pos)):
		gene_beg_dist.append(get_dist_exon_gene_beg(exon_to_genes, genes_data, pos[i]))
	return gene_beg_dist

class Gene_expr:
	def __init__(self, name, tumor, tumor_num, tumor_setd2, tumor_setd2_num, normal, normal_num):
		self.name = name
		self.tumor = tumor
		self.tumor_num = tumor_num
		self.tumor_setd2 = tumor_setd2
		self.tumor_setd2_num = tumor_setd2_num
		self.normal = normal
		self.normal_num = normal_num

def process_gene_expression_file(gene_expression_fn, sample_classification_fn):
	gene_expression = {}
	fin = open(gene_expression_fn)
	header = fin.readline()
	fin.readline()
	for line in fin:
		pass
	return gene_expression

def get_gene_expression(gene_expression_fn, exon_to_genes, pos, sample_classification_fn):
	gene_expression = process_gene_expression_file(gene_expression_fn, sample_classification_fn)
	return

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
	
	(genes_data, exon_dict, interval_trees) = parse_genes_file(genes_fn)
	print 'Total genes numer', len(genes_data)

	data_dir_list = [os.path.join(data_dir, d) for d in os.listdir(data_dir)]
	for d in data_dir_list:
		print 'processing', os.path.basename(d)
		in_fn = os.path.join(d, os.path.basename(d) + '__t_setd2__t__n.txt')
		(pos, tumor_setd2_broken_num, tumor_setd2_broken, tumor_num, tumor, normal_num, normal) = read_data(in_fn)
		exon_to_genes = exon_to_gene_names(exon_dict, interval_trees, pos)

		gene_beg_dist = get_gene_dist(exon_to_genes, genes_data, pos)		

		gene_expression_fn = None
		for elem in os.listdir(d):
			if 'gene_expression' in elem:
				gene_expression_fn = os.path.join(d, elem)
		if not (gene_expression_fn):
			continue
		sample_classification_fn = os.path.join(d, os.path.basename(d) + '_sample_classification.txt') # generated on step 2
		gene_expression = get_gene_expression(gene_expression_fn, exon_to_genes, pos, sample_classification_fn)


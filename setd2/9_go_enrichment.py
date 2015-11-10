#!/usr/bin/env python
# -*- coding: utf-8 -*- 
import os
import sys
import argparse
import subprocess
import numpy
from enum import Enum, unique
from sets import Set

class Exon:
	def __init__(self, exon_name, psi_norma, psi_tumor, psi_tumor_setd2):
		self.name = exon_name
		self.psi_tumor_setd2 = psi_tumor_setd2
		self.psi_tumor = psi_tumor
		self.psi_norma = psi_norma

class Gene:
	def __init__(self, name):
		self.name = name
		self.exons = []

	def add_exon(self, exon):
		self.exons.append(exon)

def parse_psi_file(fin):
	f = open(fin)
	header = f.readline()
	exons_lst = []
	for line in f:
		(name, psi_tumor_setd2, psi_tumor, psi_norma) = line.split('\t')
		psi_tumor_setd2 = float(psi_tumor_setd2)
		psi_tumor = float(psi_tumor)
		psi_norma = float(psi_norma)
		e = Exon(name, psi_norma, psi_tumor, psi_tumor_setd2)
		exons_lst.append(e)
	f.close()
	return exons_lst

def parse_expr_and_dist_file(fin, exons):
	f = open(fin)
	header = f.readline()
	genes = {}
	line_num = 0
	for line in f:
		(name, _, _, _, _, _, gene_name) = line.split('\t')
		gene_name = gene_name.strip()
		if gene_name != '' and len(gene_name.split(',')) == 1: # ignore exons corresponding to multiple genes
			if not genes.has_key(gene_name):
				genes[gene_name] = Gene(gene_name)
			if exons[line_num].name == name:
				genes[gene_name].add_exon(exons[line_num])
		line_num += 1
	f.close()
	return genes

def prepare_gene_list_for_go_enrichment(genes, threshold, wt_setd2_fn, norma_wt_fn, norma_setd2_fn):
	wt_setd2_f = open(wt_setd2_fn, 'w')
	norma_setd2_f = open(norma_setd2_fn, 'w')
	norma_wt_f = open(norma_wt_fn, 'w')
	for (gene_name, g) in genes.iteritems():
		name = gene_name.split('|')[1] # entrez number only
		norma_setd2 = [abs(e.psi_norma - e.psi_tumor_setd2) for e in g.exons]
		norma_setd2 = [e for e in norma_setd2 if ~numpy.isnan(e)]
		if len(norma_setd2) > 0:
			if max(norma_setd2) > threshold:
				norma_setd2_f.write(name + '\t' + str(max(norma_setd2)) + '\n')
		norma_wt = [abs(e.psi_norma - e.psi_tumor) for e in g.exons]
		norma_wt = [e for e in norma_wt if ~numpy.isnan(e)]
		if len(norma_wt) > 0:
			if max(norma_wt) > threshold:
				norma_wt_f.write(name + '\t' + str(max(norma_wt)) + '\n')
		wt_setd2 = [abs(e.psi_tumor - e.psi_tumor_setd2) for e in g.exons]
		wt_setd2 = [e for e in wt_setd2 if ~numpy.isnan(e)]
		if len(wt_setd2) > 0:
			if max(wt_setd2) > threshold:
				wt_setd2_f.write(name + '\t' + str(max(wt_setd2)) + '\n')
	wt_setd2_f.close()
	norma_setd2_f.close()
	norma_wt_f.close()

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print 'Usage:', sys.argv[0], '-d <data directory>'
		exit()

	parser = argparse.ArgumentParser(prog = sys.argv[0], description='Analyze')
	parser.add_argument('-d', '--data_dir', help='data directory', required=True)
	args = parser.parse_args()
	data_dir = args.data_dir

	data_dir_list = [os.path.join(data_dir, d) for d in os.listdir(data_dir)]
	for d in data_dir_list:
		print 'processing', os.path.basename(d)
		expr_dist_fn = os.path.join(d, os.path.basename(d) + '_expression_and_dist.txt')
		psi_fn = os.path.join(d, os.path.basename(d) + '_PSI_average.txt')
		exons = parse_psi_file(psi_fn)
		genes = parse_expr_and_dist_file(expr_dist_fn, exons)
		wt_setd2_fn = os.path.join(d, os.path.basename(d) + '_go_wt_setd2.txt')
		norma_wt_fn = os.path.join(d, os.path.basename(d) + '_go_norma_wt.txt')
		norma_setd2_fn = os.path.join(d, os.path.basename(d) + '_go_norma_setd2.txt')
		threshold = 0.2
		prepare_gene_list_for_go_enrichment(genes, threshold, wt_setd2_fn, norma_wt_fn, norma_setd2_fn)


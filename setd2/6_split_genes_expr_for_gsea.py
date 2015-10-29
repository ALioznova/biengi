#!/usr/bin/env python
# -*- coding: utf-8 -*- 
import os
import sys
import argparse
import subprocess
import numpy
from enum import Enum, unique
from sets import Set

def is_tool(name):
	try:
		devnull = open(os.devnull)
		subprocess.Popen([name], stdout=devnull, stderr=devnull).communicate()
	except OSError as e:
		if e.errno == os.errno.ENOENT:
			return False
	return True

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
	exons = []
	for line in f:
		(name, psi_tumor_setd2, psi_tumor, psi_norma) = line.split('\t')
		psi_tumor_setd2 = float(psi_tumor_setd2)
		psi_tumor = float(psi_tumor)
		psi_norma = float(psi_norma)
		exons.append(Exon(name, psi_norma, psi_tumor, psi_tumor_setd2))
	f.close()
	return exons

def parse_expr_and_dist_file(fin, exons):
	f = open(fin)
	header = f.readline()
	genes = {}
	line_num = 0
	for line in f:
		(name, _, _, _, _, _, gene_name) = line.split('\t')
		gene_name = gene_name.strip()
		if gene_name != '':
			if not genes.has_key(gene_name):
				genes[gene_name] = Gene(gene_name)
			if exons[line_num].name == name:
				genes[gene_name].add_exon(exons[line_num])
		line_num += 1
	f.close()
	return genes

def output_info_for_gsea(data_fn, phenotype_fn, genes):
	can_run_gsea = False
	for (gene_name, g) in genes.iteritems():
		tumor_norm = []
		tumor_setd2_norm = []
		for e in g.exons:
			tumor_norm.append(abs(e.psi_tumor_setd2 - e.psi_norma))
			tumor_setd2_norm.append(abs(e.psi_tumor - e.psi_norma))
		if len(tumor_norm) == 0 and len(tumor_setd2_norm) == 0:
			continue
		tumor_norm = numpy.nanmax(tumor_norm)
		tumor_setd2_norm = numpy.nanmax(tumor_setd2_norm)
		if ~numpy.isnan(tumor_norm) and ~numpy.isnan(tumor_setd2_norm):
			can_run_gsea = True
			break
	if not can_run_gsea:
		return False
	#http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats
	fout = open(data_fn, 'w')
	fout.write('Name\tDescription\ttumor_norm_delta_psi\ttumor_setd2_norm_delta_psi\n')
	for (gene_name, exons) in genes.iteritems():
		fout.write(gene_name.split('|')[1] + '\tna') # entrez number only
		tumor_norm = []
		tumor_setd2_norm = []
		tumor_norm_max = None
		tumor_setd2_norm_max = None
		for e in g.exons:
			tumor_norm.append(abs(e.psi_tumor_setd2 - e.psi_norma))
			tumor_setd2_norm.append(abs(e.psi_tumor - e.psi_norma))
		if len(tumor_norm) != 0 and len(tumor_setd2_norm) != 0:
			tumor_norm_max = numpy.nanmax(tumor_norm)
			tumor_setd2_norm_max = numpy.nanmax(tumor_setd2_norm)
		fout.write('\t')
		if tumor_norm_max and ~numpy.isnan(tumor_norm_max):
			fout.write(str(tumor_norm_max))			
		fout.write('\t')
		if tumor_setd2_norm_max and ~numpy.isnan(tumor_setd2_norm_max):
			fout.write(str(tumor_setd2_norm_max))			
		fout.write('\n')
	fout.close()
	fout = open(phenotype_fn, 'w')
	fout.write('2 2 1\n')
	fout.write('# tumor_norm_delta_psi tumor_setd2_norm_delta_psi\n')
	fout.write('tumor_norm_delta_psi tumor_setd2_norm_delta_psi\n')
	fout.close()
	return True

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print 'Usage:', sys.argv[0], '-d <data directory> -gmx <gene set file> -gsea <gsea path>'
		exit()

	parser = argparse.ArgumentParser(prog = sys.argv[0], description='Analyze')
	parser.add_argument('-d', '--data_dir', help='data directory', required=True)
	parser.add_argument('--gmx', help='gene set file', required=True)
	parser.add_argument('--gsea', help='gsea path', required=True)
	args = parser.parse_args()
	data_dir = args.data_dir
	gmx = args.gmx
	gsea = args.gsea

	if not is_tool(gsea):
		print >> sys.stderr, "Not a tool", gsea
		sys.exit(1)
	
	data_dir_list = [os.path.join(data_dir, d) for d in os.listdir(data_dir)]
	for d in data_dir_list:
		print 'processing', os.path.basename(d)

		expr_dist_fn = os.path.join(d, os.path.basename(d) + '_expression_and_dist.txt')
		psi_fn = os.path.join(d, os.path.basename(d) + '_PSI_average.txt')
		exons = parse_psi_file(psi_fn)
		genes = parse_expr_and_dist_file(expr_dist_fn, exons)
		data_fn = os.path.join(d, os.path.basename(d) + '_gsea_data.txt')
		phenotype_fn = os.path.join(d, os.path.basename(d) + '_gsea_phenotype.txt')
		can_run_gsea = output_info_for_gsea(data_fn, phenotype_fn, genes)
		if not can_run_gsea:
			continue
		out_dir = os.path.join(d, os.path.basename(d) + '_gsea')
		#http://www.broadinstitute.org/gsea/doc/GSEAUserGuideTEXT.htm#_Run_GSEA_Page
		#http://www.broadinstitute.org/gsea/doc/GSEAUserGuideTEXT.htm#_Metrics_for_Ranking
		subprocess.call(['java', '-cp', gsea, '-Xmx512m', 'xtools.gsea.Gsea', '-res', data_fn, '-cls', phenotype_fn, '-gmx', gmx, '-out', out_dir, '-nperm', '10', '-collapse', 'false', '-metric', 'Diff_of_Classes'])
		#http://www.broadinstitute.org/gsea/doc/GSEAUserGuideTEXT.htm#_Interpreting_GSEA_Results



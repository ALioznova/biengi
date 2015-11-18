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

def output_info_for_gsea(data_fn, phenotype_fn, genes):
	can_run_gsea = False
	for (gene_name, g) in genes.iteritems():
		t_setd2_n = [abs(e.psi_tumor_setd2 - e.psi_norma) for e in g.exons]
		t_setd2_n = [e for e in t_setd2_n if ~numpy.isnan(e)]
		t_n = [abs(e.psi_tumor - e.psi_norma) for e in g.exons]
		t_n = [e for e in t_n if ~numpy.isnan(e)]
		if len(t_setd2_n) > 0 and len(t_n) > 0:
			can_run_gsea = True
			break
	if not can_run_gsea:
		return False
	#http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats
	fout = open(data_fn, 'w')
	fout.write('Name\tDescription\ttumor_norm_delta_psi\ttumor_setd2_norm_delta_psi\n')
	for (gene_name, g) in genes.iteritems():
		new_line = ''
		new_line += (gene_name.split('|')[1] + '\tna') # entrez number only
		t_setd2_n = [abs(e.psi_tumor_setd2 - e.psi_norma) for e in g.exons]
		t_setd2_n = [e for e in t_setd2_n if ~numpy.isnan(e)]
		t_n = [abs(e.psi_tumor - e.psi_norma) for e in g.exons]
		t_n = [e for e in t_n if ~numpy.isnan(e)]
		if len(t_n) > 0:
			tumor_norm_max = max(t_n)
		else:
			tumor_norm_max = None
		if len(t_setd2_n) > 0:
			tumor_setd2_norm_max = max(t_setd2_n)
		else:
			tumor_setd2_norm_max = None
		new_line += ('\t')
		if tumor_norm_max:
			new_line += (str(tumor_norm_max))
		new_line += ('\t')
		if tumor_setd2_norm_max:
			new_line += (str(tumor_setd2_norm_max))
		if (not tumor_norm_max) or (not tumor_setd2_norm_max):
			new_line = None
		if new_line:
			fout.write(new_line + '\n')
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
	parser.add_argument('-o', '--out_dir', help='output directory', required=True)
	parser.add_argument('--gmx', help='gene set file', required=True)
	parser.add_argument('--gsea', help='gsea path', required=True)
	args = parser.parse_args()
	data_dir = args.data_dir
	gmx = args.gmx
	gsea = args.gsea
	out_d = args.out_dir

	if not is_tool(gsea):
		print >> sys.stderr, "Not a tool", gsea
		sys.exit(1)

	if not os.path.exists(os.path.join(out_d, 'gsea')):
		os.mkdir(os.path.join(out_d, 'gsea'))
	out_d = os.path.join(out_d, 'gsea')

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
		out_dir = os.path.join(out_d, os.path.basename(d) + '_gsea')
		#http://www.broadinstitute.org/gsea/doc/GSEAUserGuideTEXT.htm#_Run_GSEA_Page
		#http://www.broadinstitute.org/gsea/doc/GSEAUserGuideTEXT.htm#_Metrics_for_Ranking
		subprocess.call(['java', '-cp', gsea, '-Xmx2000m', 'xtools.gsea.Gsea', '-res', data_fn, '-cls', phenotype_fn, '-gmx', gmx, '-out', out_dir, '-nperm', '10', '-collapse', 'false', '-metric', 'Diff_of_Classes'])
		#http://www.broadinstitute.org/gsea/doc/GSEAUserGuideTEXT.htm#_Interpreting_GSEA_Results



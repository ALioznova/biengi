#!/usr/bin/env python
# -*- coding: utf-8 -*- 
import os
import sys
import argparse
import subprocess
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

class Sample_type(Enum):
	unknown = -1
	norma = 0
	tumor_wild_type = 1
	tumor_setd2 = 2

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
		gene_expression_fn = None
		phenotype_fn = None
		for elem in os.listdir(d):
			if 'genes_expression.txt' in elem:
				gene_expression_fn = os.path.join(d, elem)
			if 'phenotype.cls' in elem:
				phenotype_fn = os.path.join(d, elem)
		if gene_expression_fn and phenotype_fn:
			#http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats
			phenotype_f = open(phenotype_fn)
			phenotype_f.readline()
			phenotype_f.readline()
			phenotype_data_str = phenotype_f.readline().split()
			phenotype_data = []
			for elem in phenotype_data_str:
				if elem == 'tumor_wild_type':
					phenotype_data.append(Sample_type.tumor_wild_type)
				elif elem == 'tumor_setd2':
					phenotype_data.append(Sample_type.tumor_setd2)
				elif elem == 'norma':
					phenotype_data.append(Sample_type.norma)
				else:
					phenotype_data.append(Sample_type.unknown)
			phenotype_f.close()
			phen_tumor_wt_vs_setd2 = []
			for i in xrange(len(phenotype_data)):
				if phenotype_data[i] == Sample_type.tumor_wild_type or phenotype_data[i] == Sample_type.tumor_setd2:
					phen_tumor_wt_vs_setd2.append(phenotype_data[i])
			if len(Set(phen_tumor_wt_vs_setd2)) != 2:
				continue
			cls = os.path.join(d, os.path.basename(d) + '_phenotype_tumor_wt_vs_setd2.cls')
			fout_phen = open(cls, 'w')
			fout_phen.write(str(len(phen_tumor_wt_vs_setd2)) + ' ' + str(len(Set(phen_tumor_wt_vs_setd2))) + ' 1\n')
			fout_phen.write('#type1 type2\n')
			types_line = ''
			for st in phen_tumor_wt_vs_setd2:
				if st == Sample_type.norma:
					types_line += 'norma '
				elif st == Sample_type.tumor_wild_type:
					types_line += 'tumor_wild_type '
				elif st == Sample_type.tumor_setd2:
					types_line += 'tumor_setd2 '
				else:
					types_line += 'unknown '
			fout_phen.write(types_line[:-1] + '\n')
			fout_phen.close()
			gene_expression_f = open(gene_expression_fn)
			res = os.path.join(d, os.path.basename(d) + '_genes_expression_tumor_wt_vs_setd2.txt')
			fout_expr = open(res, 'w')
			sample_names = gene_expression_f.readline().split()[2:]
			fout_expr.write('Name\tDescription')
			for i in xrange(len(sample_names)):
				if phenotype_data[i] == Sample_type.tumor_wild_type or phenotype_data[i] == Sample_type.tumor_setd2:
					fout_expr.write('\t' + sample_names[i])
			fout_expr.write('\n')
			for line in gene_expression_f:
				gene_name = line.split()[0]
				fout_expr.write(gene_name + '\tna')
				expr_data = line.split()[2:]
				for i in xrange(len(phenotype_data)):
					if phenotype_data[i] == Sample_type.tumor_wild_type or phenotype_data[i] == Sample_type.tumor_setd2:
						fout_expr.write('\t' + expr_data[i])
				fout_expr.write('\n')
			gene_expression_f.close()
			fout_expr.close()

			out_dir = os.path.join(d, os.path.basename(d) + '_gsea')

			#http://www.broadinstitute.org/gsea/doc/GSEAUserGuideTEXT.htm#_Run_GSEA_Page
			#http://www.broadinstitute.org/gsea/doc/GSEAUserGuideTEXT.htm#_Metrics_for_Ranking
			subprocess.call(['java', '-cp', gsea, '-Xmx512m', 'xtools.gsea.Gsea', '-res', res, '-cls', cls, '-gmx', gmx, '-out', out_dir, '-nperm', '10', '-collapse', 'false', '-metric', 'Diff_of_Classes'])


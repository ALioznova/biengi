#!/usr/bin/env python
# -*- coding: utf-8 -*- 
import os
import sys
import argparse
import numpy
from sets import Set
from intervaltree import Interval, IntervalTree
from enum import Enum, unique
import random

class Pos_rec:
	def __init__(self, description):
		(chr_name, begin1, end1, begin2, end2, strand) = description.split(':')
		self.desc = description
		self.chr_name = chr_name
		self.beg1 = int(begin1)
		self.end1 = int(end1)
		self.beg2 = int(begin2)
		self.end2 = int(end2)
		if strand == '+':
			self.strand = True
		else:
			self.strand = False

class Categorized_data_frame:
	def __init__(self, data, samples_num):
		self.data = data
		self.samples_num = samples_num
		self.current_arr_to_call_mean = []
		
	def count_mean(self):
		if len(self.current_arr_to_call_mean) != 0:
			self.data.append(numpy.nanmean(self.current_arr_to_call_mean))
		else:
			self.data.append(float('nan'))
		self.current_arr_to_call_mean = []

class Gene_expr:
	def __init__(self, data, samples_num):
		self.data = data
		self.samples_num = samples_num
		self.current_arr_to_call_mean = []
		
	def count_mean(self):
		if len(self.current_arr_to_call_mean) != 0:
			self.data = numpy.nanmean(self.current_arr_to_call_mean)
		else:
			self.data = float('nan')
		self.current_arr_to_call_mean = []

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

class OrderedEnum(Enum):
	 def __ge__(self, other):
		 if self.__class__ is other.__class__:
			 return self.value >= other.value
		 return NotImplemented
	 def __gt__(self, other):
		 if self.__class__ is other.__class__:
			 return self.value > other.value
		 return NotImplemented
	 def __le__(self, other):
		 if self.__class__ is other.__class__:
			 return self.value <= other.value
		 return NotImplemented
	 def __lt__(self, other):
		 if self.__class__ is other.__class__:
			 return self.value < other.value
		 return NotImplemented

@unique
class Impact(OrderedEnum):
	high = 3
	moderate = 2
	low = 1
	modifier = 0
	no = -1
	unknown = -2

class Sample_type(Enum):
	unknown = -1
	norma = 0
	tumor_wild_type = 1
	tumor_unclassified = 2
	tumor_mutant = 3

def read_psi_average_data(in_fn):
	in_f = open(in_fn)
	header_line = in_f.readline()
	tumor_mutant_num = int((header_line.split()[1]).split(':')[1])
	tumor_wild_type_num = int((header_line.split()[2]).split(':')[1])
	normal_num = int((header_line.split()[3]).split(':')[1])
	categorized_psi = {Sample_type.norma : Categorized_data_frame([], normal_num), Sample_type.tumor_wild_type : Categorized_data_frame([], tumor_wild_type_num), Sample_type.tumor_mutant : Categorized_data_frame([], tumor_mutant_num)}
	pos = []
	for line in in_f:
		categorized_psi[Sample_type.tumor_mutant].data.append(float(line.split()[1]))
		categorized_psi[Sample_type.tumor_wild_type].data.append(float(line.split()[2]))
		categorized_psi[Sample_type.norma].data.append(float(line.split()[3]))
		pos.append(Pos_rec(line.split()[0]))
	in_f.close()
	return (pos, categorized_psi)

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
	gf.close()
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
		return ([], [], [])
	exon_genes = list(exon_to_genes[(exon_chr, exon_strand)][(exon_beg, exon_end)])
	dist_bp = []
	dist_perc = []
	gene_name = []
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
							if loc.strand == True and exon_strand == True:
								if exon_beg < loc.beg or exon_beg > loc.end:
									continue
								dist_bp.append(exon_beg - loc.beg)
								dist_perc.append(float(exon_beg - loc.beg)/(loc.end - loc.beg))
								gene_name.append(gene.name)
							elif loc.strand == False and exon_strand == False:
								if loc.end < exon_end or exon_end < loc.beg:
									continue
								dist_bp.append(loc.end - exon_end)
								dist_perc.append(float(loc.end - exon_end)/(loc.end - loc.beg))
								gene_name.append(gene.name)
							else:
								continue
	return (dist_bp, dist_perc, gene_name)

def get_gene_dist(exon_to_genes, genes_data, pos):
	gene_beg_dist_bp = []
	gene_beg_dist_perc = []
	gene_names = []
	for i in xrange(len(pos)):
		(dist_bp, dist_perc, gene_name) = get_dist_exon_gene_beg(exon_to_genes, genes_data, pos[i])
		dist_bp = list(Set(dist_bp))
		dist_perc = list(Set(dist_perc))
		gene_name = list(Set(gene_name))
		gene_beg_dist_bp.append(dist_bp)
		gene_beg_dist_perc.append(dist_perc)
		gene_names.append(gene_name)
	return (gene_beg_dist_bp, gene_beg_dist_perc, gene_names)

def get_sample_classification(sample_names, samples_classification_fn):
	sample_classification = {}
	gene_mutation_rates = {}
	fin = open(samples_classification_fn)
	for line in fin:
		s_name = line.split()[0]
		sample_type = (line.split()[1]).strip()
		if sample_type == 'Impact.high':
			gene_mutation_rates[s_name] = Impact.high
		elif sample_type == 'Impact.moderate':
			gene_mutation_rates[s_name] = Impact.moderate
		elif sample_type == 'Impact.low':
			gene_mutation_rates[s_name] = Impact.low
		elif sample_type == 'Impact.modifier':
			gene_mutation_rates[s_name] = Impact.modifier
		elif sample_type == 'Impact.no':
			gene_mutation_rates[s_name] = Impact.no
		elif sample_type == 'Impact.unknown':
			gene_mutation_rates[s_name] = Impact.unknown
		else:
			gene_mutation_rates[s_name] = Impact.unknown
			print 'Unknown type', sample_type
	fin.close()
	# barcodes https://wiki.nci.nih.gov/display/TCGA/TCGA+Barcode
	# code tables: https://tcga-data.nci.nih.gov/datareports/codeTablesReport.htm?codeTable=sample%20type
	sample_type = [int(s.split('-')[3][:2]) for s in sample_names]
	for i in xrange(len(sample_type)):
		if (sample_type[i] >= 1) and (sample_type[i] <= 9): # tumor
			if gene_mutation_rates.has_key(sample_names[i][:15]):
				if gene_mutation_rates[sample_names[i][:15]] == Impact.high:
					sample_classification[sample_names[i]] = Sample_type.tumor_mutant
				elif gene_mutation_rates[sample_names[i][:15]] == Impact.no:
					sample_classification[sample_names[i]] = Sample_type.tumor_wild_type
				else:
					sample_classification[sample_names[i]] = Sample_type.tumor_unclassified
			else:
				sample_classification[sample_names[i]] = Sample_type.tumor_unclassified
		elif (sample_type[i] >= 10) and (sample_type[i] <= 19): # norma
			sample_classification[sample_names[i]] = Sample_type.norma
		else:
			sample_classification[sample_names[i]] = Sample_type.unknown
	return sample_classification

def process_gene_expression_file(gene_expression_fn, samples_classification_fn, filter_for_equal_size = False):
	fin = open(gene_expression_fn)
	sample_names = fin.readline().split('\t')[1:]
	sample_classification = get_sample_classification(sample_names, samples_classification_fn)
	if filter_for_equal_size:
		sample_type = []
		for name in sample_names:
			sample_type.append(sample_classification[name])
		indexes = {Sample_type.norma : [], Sample_type.tumor_wild_type : [], Sample_type.tumor_mutant : []}
		for i in xrange(len(sample_type)):
			if sample_type[i] in indexes.keys():
				indexes[sample_type[i]].append(i)
		equal_size = min(len(indexes[Sample_type.tumor_wild_type]), len(indexes[Sample_type.tumor_mutant]))
		random.shuffle(indexes[Sample_type.tumor_wild_type])
		random.shuffle(indexes[Sample_type.tumor_mutant])
		indexes[Sample_type.tumor_wild_type] = indexes[Sample_type.tumor_wild_type][:equal_size]
		indexes[Sample_type.tumor_mutant] = indexes[Sample_type.tumor_mutant][:equal_size]
	fin.readline()
	gene_expression = {}
	for line in fin:
		gene_name = line.split('\t')[0]
		if gene_name.endswith('_calculated'):
			gene_name = gene_name[:-len('_calculated')]
		data = line.split('\t')[1:]
		categorized_expr = {Sample_type.norma : Gene_expr(float('nan'), 0), Sample_type.tumor_wild_type : Gene_expr(float('nan'), 0), Sample_type.tumor_mutant : Gene_expr(float('nan'), 0)}
		for i in xrange(len(data)):
			if i % 3 == 2 and categorized_expr.has_key(sample_classification[sample_names[i]]):
				cur_type = sample_classification[sample_names[i]]
				if filter_for_equal_size:
					for key in indexes.iterkeys():
						if i in indexes[key]:
							categorized_expr[cur_type].current_arr_to_call_mean.append(float(data[i]))
							categorized_expr[cur_type].samples_num += 1
				else:
					categorized_expr[cur_type].current_arr_to_call_mean.append(float(data[i]))
					categorized_expr[cur_type].samples_num += 1
		for (cur_type, cur_expr_data) in categorized_expr.iteritems():
			cur_expr_data.count_mean()
		gene_expression[gene_name] = categorized_expr
	fin.close()
	return gene_expression

def get_exons_gene_expression(gene_expression_fn, exon_to_genes, pos, samples_classification_fn):
	gene_expression_data = process_gene_expression_file(gene_expression_fn, samples_classification_fn, True)
	categorized_expr = {Sample_type.norma : [], Sample_type.tumor_wild_type : [], Sample_type.tumor_mutant : []}
	for p in pos:
		for sample_type, expr_data in categorized_expr.iteritems():
			categorized_expr[sample_type].append([])
		if exon_to_genes[(p.chr_name, p.strand)].has_key((p.end1, p.beg2)):
			for gene in list(exon_to_genes[(p.chr_name, p.strand)][(p.end1, p.beg2)]):
				if not gene_expression_data.has_key(gene):
					continue
				for sample_type, expr_data in categorized_expr.iteritems():
					categorized_expr[sample_type][-1].append(gene_expression_data[gene][sample_type].data)
	return categorized_expr

def output_expression_and_dist_file(out_fn, pos, gene_beg_dist_bp, gene_beg_dist_perc, gene_expression_fn, categorized_expr, gene_names):
	fout = open(out_fn, 'w')
	fout.write('Pos:\tDist_to_gene_start_in_bp:\tDist_to_gene_start_in_percent_of_gene_len:\tGene_expression_in_tumor_mutant:\tGene_expression_in_tumor_wild_type:\tGene_expression_in_normal:\tGene_name:\n')
	for i in xrange(len(pos)):
		new_line = pos[i].desc + '\t'
		for dst in gene_beg_dist_bp[i]:
			new_line += (str(dst) + ',')
		if len(gene_beg_dist_bp[i]) != 0:
			new_line = new_line[:-1]
		new_line += '\t'
		for dst in gene_beg_dist_perc[i]:
			new_line += (str(dst) + ',')
		if len(gene_beg_dist_perc[i]) != 0:
			new_line = new_line[:-1]
		if gene_expression_fn:
			new_line += '\t'
			for expr in categorized_expr[Sample_type.tumor_mutant][i]:
				new_line += (str(expr) + ',')
			if len(categorized_expr[Sample_type.tumor_mutant][i]) != 0:
				new_line = new_line[:-1]
			new_line += '\t'
			for expr in categorized_expr[Sample_type.tumor_wild_type][i]:
				new_line += (str(expr) + ',')
			if len(categorized_expr[Sample_type.tumor_wild_type][i]) != 0:
				new_line = new_line[:-1]
			new_line += '\t'
			for expr in categorized_expr[Sample_type.norma][i]:
				new_line += (str(expr) + ',')
			if len(categorized_expr[Sample_type.norma][i]) != 0:
				new_line = new_line[:-1]
		else:
			new_line += '\t\t\t'
		new_line += '\t'
		for name in gene_names[i]:
			new_line += (name + ',')
		if len(gene_names[i]) != 0:
			new_line = new_line[:-1]
		new_line += '\n'
		fout.write(new_line)
	fout.close()

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print 'Usage:', sys.argv[0], '-d <data directory> -c <computations directory> -g <gene annotation file> -m <mutant gene name>'
		exit()

	parser = argparse.ArgumentParser(prog = sys.argv[0], description='Find distances from exons to genes and corresponding gene expression')
	parser.add_argument('-d', '--data_dir', help='data directory', required=True)
	parser.add_argument('-c', '--comp_dir', help='computations directory', required=True)
	parser.add_argument('-g', '--genes', help='gene annotation file', required=True)
	parser.add_argument('-m', '--mut_gene', help='mutatn gene name', required=True)
	args = parser.parse_args()
	data_dir = args.data_dir
	comp_dir = args.comp_dir
	genes_fn = args.genes
	mutant_gene = args.mut_gene
	
	if not os.path.isfile(genes_fn):
		print 'No such file', genes_fn
		sys.exit(1)
	
	(genes_data, exon_dict, interval_trees) = parse_genes_file(genes_fn)
	print 'Total genes numer', len(genes_data)

	random.seed(1)

	data_dir_list = [os.path.join(data_dir, d) for d in os.listdir(data_dir)]
	for d in data_dir_list:
		print 'processing', os.path.basename(d)
		d_comp = os.path.join(comp_dir, os.path.basename(d))
		in_fn = os.path.join(os.path.join(d_comp, mutant_gene), os.path.basename(d) + '_PSI_averaged_by_' + mutant_gene + '.txt')

		(pos, categorized_psi) = read_psi_average_data(in_fn)
		exon_to_genes = exon_to_gene_names(exon_dict, interval_trees, pos)
		(gene_beg_dist_bp, gene_beg_dist_perc, gene_names) = get_gene_dist(exon_to_genes, genes_data, pos)

		gene_expression_fn = None
		categorized_expr = None
		for elem in os.listdir(d):
			if 'gene_expression' in elem:
				gene_expression_fn = os.path.join(d, elem)
		if gene_expression_fn:
			samples_classification_fn = os.path.join(os.path.join(d_comp, mutant_gene), os.path.basename(d) + '_' + mutant_gene + '_mutation_impact_for_samples.txt')
			categorized_expr = get_exons_gene_expression(gene_expression_fn, exon_to_genes, pos, samples_classification_fn)
		out_fn = os.path.join(os.path.join(d_comp, mutant_gene), os.path.basename(d) + '_expression_and_dist_split_by_' + mutant_gene + '.txt')
		output_expression_and_dist_file(out_fn, pos, gene_beg_dist_bp, gene_beg_dist_perc, gene_expression_fn, categorized_expr, gene_names)


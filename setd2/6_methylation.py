#!/usr/bin/env python
# -*- coding: utf-8 -*- 
import os
import sys
import argparse
import numpy
from sets import Set
from enum import Enum, unique
from intervaltree import Interval, IntervalTree
import random

class Pos_rec:
	def __init__(self, description):
		(chr_name, begin1, end1, begin2, end2, strand) = description.split(':')
		self.desc = description
		if chr_name.startswith('chr'):
			self.chr_name = chr_name[len('chr'):]
		else:
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

def build_exon_interval_tree(pos):
	interval_trees = {}
	for p in pos:
		if not interval_trees.has_key(p.chr_name):
			interval_trees[p.chr_name] = IntervalTree()
		interval_trees[p.chr_name][p.end1 : p.beg2] = {Sample_type.norma: [], Sample_type.tumor_wild_type: [], Sample_type.tumor_mutant : []} # doesn't include right border
	return interval_trees

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

def process_methylation_data(methylation_fn, samples_classification_fn, interval_trees, filter_for_equal_size = False):
	methylation_f = open(methylation_fn)
	sample_names = methylation_f.readline().strip().split('\t')[1:]
	methylation_f.readline()
	sample_classification = get_sample_classification(list(Set(sample_names)), samples_classification_fn)
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
	sample_numbers = {Sample_type.tumor_mutant: 0, Sample_type.tumor_wild_type: 0, Sample_type.norma: 0}
	for (s_name, s_type) in sample_classification.iteritems():
		if sample_numbers.has_key(s_type):
			sample_numbers[Sample_type.tumor_wild_type] = equal_size
			sample_numbers[Sample_type.tumor_mutant] = equal_size
	if filter_for_equal_size:
		sample_numbers[s_type] += 1
	for line in methylation_f:
		cur_meth = {Sample_type.norma: [], Sample_type.tumor_wild_type: [], Sample_type.tumor_mutant : []}
		data = line.strip().split('\t')[1:]
		cur_chr = data[2]
		cur_pos = int(data[3])
		for i in xrange(len(data)):
			if i % 4 == 0: # Beta_value
				meth_val = data[i]
				if meth_val != 'NA' and cur_meth.has_key(sample_classification[sample_names[i]]):
					if filter_for_equal_size:
						for key in indexes.iterkeys():
							if i in indexes[key]:
								cur_meth[sample_classification[sample_names[i]]].append(float(meth_val))
					else:
						cur_meth[sample_classification[sample_names[i]]].append(float(meth_val))
			elif i % 4 == 2: # chr
				if data[i] != cur_chr:
					print 'wrong chr, expectation', cur_chr, 'reality', data[i]
			elif i % 4 == 3: # pos
				if int(data[i]) != cur_pos:
					print 'wrong pos, expectation', cur_pos, 'reality', int(data[i])
		for (key, value) in cur_meth.iteritems():
			if len(value) != 0:
				cur_meth[key] = numpy.mean(value)
			else:
				cur_meth[key] = float('nan')
		if interval_trees.has_key(cur_chr):
			for (cur_type, cur_type_meth) in cur_meth.iteritems():
				for interval in interval_trees[cur_chr][cur_pos]:
					if ~numpy.isnan(cur_type_meth):
						interval.data[cur_type].append(cur_type_meth)
	methylation_f.close()
	for (chr_name, intervals) in interval_trees.iteritems():
		for interval in intervals:
			for (cur_type, cur_type_meth) in interval.data.iteritems():
				if len(cur_type_meth) == 0:
					interval.data[cur_type] = float('nan')
				else:
					interval.data[cur_type] = numpy.mean(cur_type_meth)
	return (interval_trees, sample_numbers)

def output_methylation_info(output_fn, interval_trees, sample_numbers):
	outf = open(output_fn, 'w')
	outf.write('Exon_pos\tTumor_wild_type_methylation:' + str(sample_numbers[Sample_type.tumor_wild_type]) + '\tTumor_mutant_methylation:' + str(sample_numbers[Sample_type.tumor_mutant]) + '\tNorma_methylation:' + str(sample_numbers[Sample_type.norma]) + '\n')
	for (chr_name, intervals) in interval_trees.iteritems():
		for interval in intervals:
			if not (numpy.isnan(interval.data[Sample_type.tumor_wild_type]) and numpy.isnan(interval.data[Sample_type.tumor_mutant]) and numpy.isnan(interval.data[Sample_type.norma])):
				outf.write('chr' + chr_name + ':' + str(interval.begin) + ':' + str(interval.end) + '\t' + str(interval.data[Sample_type.tumor_wild_type]) + '\t' + str(interval.data[Sample_type.tumor_mutant]) + '\t' + str(interval.data[Sample_type.norma]) + '\n')
	outf.close()

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print 'Usage:', sys.argv[0], '-d <data directory> -c <computations directory> -m <mutant gene name>'
		exit()

	parser = argparse.ArgumentParser(prog = sys.argv[0], description='Categorize methylation to normal, tumor with mut_gene mutation and tumor wild type for every exon')
	parser.add_argument('-d', '--data_dir', help='data directory', required=True)
	parser.add_argument('-c', '--comp_dir', help='computations directory', required=True)
	parser.add_argument('-m', '--mut_gene', help='mutatn gene name', required=True)
	args = parser.parse_args()
	data_dir = args.data_dir
	comp_dir = args.comp_dir
	mutant_gene = args.mut_gene
	
	data_dir_list = [os.path.join(data_dir, d) for d in os.listdir(data_dir)]
	for d in data_dir_list:
		print 'processing', os.path.basename(d)
		d_comp = os.path.join(comp_dir, os.path.basename(d))
		in_fn = os.path.join(os.path.join(d_comp, mutant_gene), os.path.basename(d) + '_PSI_averaged_by_' + mutant_gene + '.txt')
		(pos, _) = read_psi_average_data(in_fn)
		interval_trees = build_exon_interval_tree(pos)
		methylation_fn = os.path.join(d, os.path.basename(d) + '.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt')
		samples_classification_fn = os.path.join(os.path.join(d_comp, mutant_gene), os.path.basename(d) + '_' + mutant_gene + '_mutation_impact_for_samples.txt')
		(interval_trees, sample_numbers) = process_methylation_data(methylation_fn, samples_classification_fn, interval_trees, True)
		output_fn = os.path.join(os.path.join(d_comp, mutant_gene), os.path.basename(d) + '_exon_methylation_averaged_by_' + mutant_gene + '.txt')
		output_methylation_info(output_fn, interval_trees, sample_numbers)


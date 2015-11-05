#!/usr/bin/env python
# -*- coding: utf-8 -*- 
import os
import sys
import argparse
import numpy
from sets import Set
from enum import Enum, unique
from intervaltree import Interval, IntervalTree

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

def read_psi_average_data(in_fn):
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

def build_exon_interval_tree(pos):
	interval_trees = {}
	for p in pos:
		if not interval_trees.has_key(p.chr_name):
			interval_trees[p.chr_name] = IntervalTree()
		interval_trees[p.chr_name][p.end1 : p.beg2] = {Sample_type.norma: [], Sample_type.tumor_wild_type: [], Sample_type.tumor_setd2 : []} # doesn't include right border
	return interval_trees

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
	tumor_setd2 = 3

def get_sample_classification(sample_names, samples_classification_fn):
	sample_classification = {}
	setd2_mutation_rates = {}
	fin = open(samples_classification_fn)
	for line in fin:
		s_name = line.split()[0]
		sample_type = (line.split()[1]).strip()
		if sample_type == 'Impact.high':
			setd2_mutation_rates[s_name] = Impact.high
		elif sample_type == 'Impact.moderate':
			setd2_mutation_rates[s_name] = Impact.moderate
		elif sample_type == 'Impact.low':
			setd2_mutation_rates[s_name] = Impact.low
		elif sample_type == 'Impact.modifier':
			setd2_mutation_rates[s_name] = Impact.modifier
		elif sample_type == 'Impact.no':
			setd2_mutation_rates[s_name] = Impact.no
		elif sample_type == 'Impact.unknown':
			setd2_mutation_rates[s_name] = Impact.unknown
		else:
			setd2_mutation_rates[s_name] = Impact.unknown
			print 'Unknown type', sample_type
	fin.close()
	# barcodes https://wiki.nci.nih.gov/display/TCGA/TCGA+Barcode
	# code tables: https://tcga-data.nci.nih.gov/datareports/codeTablesReport.htm?codeTable=sample%20type
	sample_type = [int(s.split('-')[3][:2]) for s in sample_names]
	for i in xrange(len(sample_type)):
		if (sample_type[i] >= 1) and (sample_type[i] <= 9): # tumor
			if setd2_mutation_rates.has_key(sample_names[i][:15]):
				if setd2_mutation_rates[sample_names[i][:15]] == Impact.high:
					sample_classification[sample_names[i]] = Sample_type.tumor_setd2
				elif setd2_mutation_rates[sample_names[i][:15]] == Impact.no:
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

def read_methylation_data(methylation_f, sample_names, sample_classification, interval_trees):
	for line in methylation_f:
		cur_meth = {Sample_type.norma: [], Sample_type.tumor_wild_type: [], Sample_type.tumor_setd2 : []}
		data = line.strip().split('\t')[1:]
		cur_chr = data[2]
		cur_pos = int(data[3])
		for i in xrange(len(data)):
			if i % 4 == 0: # Beta_value
				meth_val = data[i]
				if meth_val != 'NA' and cur_meth.has_key(sample_classification[sample_names[i]]):
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
			for (key, value) in cur_meth.iteritems():
				for interval in interval_trees[cur_chr][cur_pos]:
					interval[key].data.append(value)


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
##
		if os.path.basename(d) != 'OV':
			continue
##
		print 'processing', os.path.basename(d)
		in_fn = os.path.join(d, os.path.basename(d) + '_PSI_average.txt')
		(pos, tumor_setd2_broken_num, tumor_setd2_broken, tumor_num, tumor, normal_num, normal) = read_psi_average_data(in_fn)
		interval_trees = build_exon_interval_tree(pos)
		methylation_fn = os.path.join(d, os.path.basename(d) + '.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt')
		samples_classification_fn = os.path.join(d, os.path.basename(d) + '_setd2_mutation_impact_for_samples.txt')

		methylation_f = open(methylation_fn)
		sample_names = methylation_f.readline().strip().split('\t')[1:]
		methylation_f.readline()
		sample_classification = get_sample_classification(list(Set(sample_names)), samples_classification_fn)
		read_methylation_data(methylation_f, sample_names, sample_classification, interval_trees)
		methylation_f.close()


#!/usr/bin/env python
# -*- coding: utf-8 -*- 
import os
import sys
import argparse
import pylab
import numpy
from enum import Enum, unique
import random

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

class Categorized_data_frame:
	def __init__(self, data, samples_num):
		self.data = data
		self.samples_num = samples_num

class Data_frame:
	def __init__(self, label, color, y):
		self.label = label
		self.color = color
		self.y = y
		self.num = len(y)

	def delete_elems_by_index(self, will_del_index):
		y_new = []
		for i in xrange(self.num):
			if not will_del_index[i]:
				y_new.append(self.y[i])
		self.y = y_new
		self.num = len(y_new)

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
			if i % 3 == 2 and indexes.has_key(sample_type[i]):
				indexes[sample_type[i]].append(i)
		if len(indexes[Sample_type.norma]) > 0:
			equal_size = min(len(indexes[Sample_type.norma]), len(indexes[Sample_type.tumor_wild_type]), len(indexes[Sample_type.tumor_mutant]))
			random.shuffle(indexes[Sample_type.norma])
			indexes[Sample_type.norma] = indexes[Sample_type.norma][:equal_size]
		else:
			equal_size = min(len(indexes[Sample_type.tumor_wild_type]), len(indexes[Sample_type.tumor_mutant]))
		random.shuffle(indexes[Sample_type.tumor_wild_type])
		random.shuffle(indexes[Sample_type.tumor_mutant])
		indexes[Sample_type.tumor_wild_type] = indexes[Sample_type.tumor_wild_type][:equal_size]
		indexes[Sample_type.tumor_mutant] = indexes[Sample_type.tumor_mutant][:equal_size]
	fin.readline()
	gene_expr = {Sample_type.norma : [], Sample_type.tumor_wild_type : [], Sample_type.tumor_mutant : []}
	for line in fin:
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
		for s_type in gene_expr.iterkeys():
			gene_expr[s_type].append(categorized_expr[s_type].data)
	fin.close()
	return gene_expr

def filter_data(data_frames):
	df_filtered = []
	for df in data_frames:
		only_nan = True
		for elem in df.y:
			if ~numpy.isnan(elem):
				only_nan = False
				break
		if not only_nan:
			df_filtered.append(df)
	data_frames = df_filtered
	will_del_index = []
	if len(data_frames) > 0:
		for i in xrange(data_frames[0].num):
			will_del = False
			for df in data_frames:
				if numpy.isnan(df.y[i]):
					will_del = True
			will_del_index.append(will_del)
	for df in data_frames:
		df.delete_elems_by_index(will_del_index)
	return data_frames

def draw_hist(data_frames, pic_path, plot_label, nbins):
	data = []
	colors = []
	labels = []
	max_limit = float("-inf")
	for df in data_frames:
		data.append(df.y)
		max_limit = max(max_limit, max(df.y))
		colors.append(df.color)
		labels.append(df.label)
	fig_hist = pylab.figure()
	n, bins, patches = pylab.hist(data, bins=(10 ** numpy.linspace(0.000000001, numpy.log10(max_limit), nbins)), normed=0, histtype='bar', color=colors, label=labels, log=True)
	pylab.gca().set_xscale("log")
	pylab.title(plot_label)
	pylab.legend(loc='best')
	fig_hist.savefig(pic_path)
	pylab.close(fig_hist)

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print 'Usage:', sys.argv[0], '-d <data directory> -c <computations directory> -p <pictures directory>  -m <mutant gene name>'
		exit()

	parser = argparse.ArgumentParser(prog = sys.argv[0], description='Build histogram for gene expression, data is split to tumor_mutant, tumor_wild_type and norma')
	parser.add_argument('-d', '--data_dir', help='data directory', required=True)
	parser.add_argument('-c', '--comp_dir', help='computations directory', required=True)
	parser.add_argument('-p', '--pics_dir', help='pictures directory', required=True)
	parser.add_argument('-m', '--mut_gene', help='mutatn gene name', required=True)
	args = parser.parse_args()
	data_dir = args.data_dir
	comp_dir = args.comp_dir
	pics_dir = args.pics_dir
	mutant_gene = args.mut_gene

	random.seed(4)

	gene_expr_dir = os.path.join(pics_dir, 'gene_expression_hist_for_' + mutant_gene)

	if not os.path.exists(gene_expr_dir):
		os.mkdir(gene_expr_dir)

	data_dir_list = [os.path.join(data_dir, d) for d in os.listdir(data_dir)]
	for d in data_dir_list:
		print 'processing', os.path.basename(d)
		d_comp = os.path.join(comp_dir, os.path.basename(d))
	
		gene_expression_fn = None
		for elem in os.listdir(d):
			if 'gene_expression' in elem:
				gene_expression_fn = os.path.join(d, elem)
		if gene_expression_fn:
			samples_classification_fn = os.path.join(os.path.join(d_comp, mutant_gene), os.path.basename(d) + '_' + mutant_gene + '_mutation_impact_for_samples.txt')
			gene_expr = process_gene_expression_file(gene_expression_fn, samples_classification_fn, True)
			df = [el for el in [Data_frame('Tumor_wt', 'blue', gene_expr[Sample_type.tumor_wild_type]), Data_frame('Normal', 'chartreuse', gene_expr[Sample_type.norma]), Data_frame('Tumor_mutant', 'red', gene_expr[Sample_type.tumor_mutant])] if el.num > 0]
			df = filter_data(df)
			if len(df) > 0:
				gene_expr_path = os.path.join(gene_expr_dir, os.path.basename(d) + '_gene_expr_hist.png')
				draw_hist(df, gene_expr_path, mutant_gene + '. Gene expression hist for ' + os.path.basename(d), 20)



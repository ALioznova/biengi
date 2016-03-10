#!/usr/bin/env python
# -*- coding: utf-8 -*- 
import os
import sys
import argparse
import pylab
import numpy
from enum import Enum

class Sample_type(Enum):
	norma = 0
	tumor_wild_type = 1
	tumor_mutant = 2

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

def read_expr_dist_data(in_fn):
	in_f = open(in_fn)
	header_line = in_f.readline()
	dist_bp = []
	dist_perc = []
	categorized_expr = {Sample_type.norma : Categorized_data_frame([], None), Sample_type.tumor_wild_type : Categorized_data_frame([], None), Sample_type.tumor_mutant : Categorized_data_frame([], None)}
	for line in in_f:
		dist_bp.append([int(el) for el in line.split('\t')[1].split(',') if (el != '' and el != 'nan')])
		dist_perc.append([float(el) for el in line.split('\t')[2].split(',') if (el != '' and el != 'nan')])
		categorized_expr[Sample_type.tumor_mutant].data.append([float(el) for el in line.split('\t')[3].split(',') if el.strip() != '' and el.strip() != 'nan'])
		categorized_expr[Sample_type.tumor_wild_type].data.append([float(el) for el in line.split('\t')[4].split(',') if el.strip() != '' and el.strip() != 'nan'])
		categorized_expr[Sample_type.norma].data.append([float(el) for el in line.split('\t')[5].strip().split(',') if el.strip() != '' and el.strip() != 'nan'])
	in_f.close()
	return (dist_bp, dist_perc, categorized_expr)

def read_exon_expr_data(in_fn):
	in_f = open(in_fn)
	header_line = in_f.readline()
	categorized_expr = {Sample_type.norma : Categorized_data_frame([], None), Sample_type.tumor_wild_type : Categorized_data_frame([], None), Sample_type.tumor_mutant : Categorized_data_frame([], None)}
	for line in in_f:
		categorized_expr[Sample_type.tumor_wild_type].data.append([float(el) for el in line.split('\t')[1].split(',') if el.strip() != '' and el.strip() != 'nan'])
		categorized_expr[Sample_type.tumor_mutant].data.append([float(el) for el in line.split('\t')[2].split(',') if el.strip() != '' and el.strip() != 'nan'])
		categorized_expr[Sample_type.norma].data.append([float(el) for el in line.split('\t')[3].strip().split(',') if el.strip() != '' and el.strip() != 'nan'])
	in_f.close()
	return categorized_expr

def prepare_data_for_plot(y_arr_of_arr):
	y = []
	for i in xrange(len(y_arr_of_arr)):
		if len(y_arr_of_arr[i]) != 1: # ignore exons corresponding to multiple genes
			y.append(float('nan'))
			continue
		y.append(y_arr_of_arr[i][0])
	return y

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
	for i in xrange(data_frames[0].num):
		will_del = False
		for df in data_frames:
			if numpy.isnan(df.y[i]):
				will_del = True
		will_del_index.append(will_del)
	for df in data_frames:
		df.delete_elems_by_index(will_del_index)

def draw_hist(data_frames, pic_path, plot_label, nbins):
	data = []
	colors = []
	labels = []
	for df in data_frames:
		data.append(df.y)
		colors.append(df.color)
		labels.append(df.label)
	fig_hist = pylab.figure()
	n, bins, patches = pylab.hist(data, bins=nbins, normed=0, histtype='bar', color=colors, label=labels, log=True)
	pylab.title(plot_label)
	pylab.legend(loc='best')
	fig_hist.savefig(pic_path)
	pylab.close(fig_hist)

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print 'Usage:', sys.argv[0], '-c <computations directory> -p <pictures directory>  -m <mutant gene name>'
		exit()

	parser = argparse.ArgumentParser(prog = sys.argv[0], description='Build histogram for expression (exons and genes), data is split to tumor_mutant, tumor_wild_type and norma')
	parser.add_argument('-c', '--comp_dir', help='computations directory', required=True)
	parser.add_argument('-p', '--pics_dir', help='pictures directory', required=True)
	parser.add_argument('-m', '--mut_gene', help='mutatn gene name', required=True)
	args = parser.parse_args()
	comp_dir = args.comp_dir
	pics_dir = args.pics_dir
	mutant_gene = args.mut_gene

	exon_expr_dir = os.path.join(pics_dir, 'gene_expression_hist_for_' + mutant_gene)
	gene_expr_dir = os.path.join(pics_dir, 'exon_expression_hist_for_' + mutant_gene)

	if not os.path.exists(exon_expr_dir):
		os.mkdir(exon_expr_dir)
	if not os.path.exists(gene_expr_dir):
		os.mkdir(gene_expr_dir)

	data_dir_list = [os.path.join(comp_dir, d) for d in os.listdir(comp_dir)]
	for d in data_dir_list:
		print 'processing', os.path.basename(d)
		expr_dist_fn = os.path.join(os.path.join(d, mutant_gene), os.path.basename(d) + '_expression_and_dist_split_by_' + mutant_gene + '.txt')
		(_, _, gene_expr) = read_expr_dist_data(expr_dist_fn)
		# check for gene expression presence
		no_gene_expr = True
		for ge in gene_expr.itervalues():
			for value in ge.data:
				if len(value) > 0:
					no_gene_expr = False
					break
		if no_gene_expr:
			gene_expr = None

		exon_expr_fn = os.path.join(os.path.join(d, mutant_gene), os.path.basename(d) + '_exon_expresion_split_by_' + mutant_gene + '.txt')
		if not os.path.isfile(exon_expr_fn):
			exon_expr = None
		else:
			exon_expr = read_exon_expr_data(exon_expr_fn)
			# check for exon expression presence
			no_exon_expr = True
			for ge in exon_expr.itervalues():
				for value in ge.data:
					if len(value) > 0:
						no_exon_expr = False
						break
			if no_exon_expr:
				exon_expr = None

		# gene expression
		if gene_expr:
			df = [el for el in [Data_frame('Tumor_wt', 'blue', prepare_data_for_plot(gene_expr[Sample_type.tumor_wild_type].data)), Data_frame('Normal', 'chartreuse', prepare_data_for_plot(gene_expr[Sample_type.norma].data)), Data_frame('Tumor_mutant', 'red', prepare_data_for_plot(gene_expr[Sample_type.tumor_mutant].data))] if el.num > 0]
			filter_data(df)
			if len(df) > 0:
				gene_expr_path = os.path.join(gene_expr_dir, os.path.basename(d) + '_gene_expr_hist.png')
				draw_hist(df, gene_expr_path, mutant_gene + '. Gene expression hist for ' + os.path.basename(d), 20)

		# exon expression
		if exon_expr:
			df = [el for el in [Data_frame('Tumor_wt', 'blue', prepare_data_for_plot(exon_expr[Sample_type.tumor_wild_type].data)), Data_frame('Normal', 'chartreuse', prepare_data_for_plot(exon_expr[Sample_type.norma].data)), Data_frame('Tumor_mutant', 'red', prepare_data_for_plot(exon_expr[Sample_type.tumor_mutant].data))] if el.num > 0]
			filter_data(df)
			if len(df) > 0:
				exon_expr_path = os.path.join(exon_expr_dir, os.path.basename(d) + '_exon_expr_hist.png')
				draw_hist(df, exon_expr_path, mutant_gene + '. Exon expression hist for ' + os.path.basename(d), 20)



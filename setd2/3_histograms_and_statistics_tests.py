#!/usr/bin/env python
# -*- coding: utf-8 -*- 
import os
import sys
import argparse
import numpy
import pylab
from scipy import stats
from scipy.stats import mstats
from enum import Enum

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

class Categorized_data_frame:
	def __init__(self, data, samples_num):
		self.data = data
		self.samples_num = samples_num

class Sample_type(Enum):
	norma = 0
	tumor_wild_type = 1
	tumor_mutant = 2

class Difference_type(Enum):
	tumor_wild_type_minus_norma = 1
	tumor_mutant_minus_norma = 2

def read_psi_average_data(in_fn):
	in_f = open(in_fn)
	header_line = in_f.readline()
	tumor_mutant_num = int((header_line.split()[1]).split(':')[1])
	tumor_wild_type_num = int((header_line.split()[2]).split(':')[1])
	normal_num = int((header_line.split()[3]).split(':')[1])
	categorized_psi = {Sample_type.norma : Categorized_data_frame([], normal_num), Sample_type.tumor_wild_type : Categorized_data_frame([], tumor_wild_type_num), Sample_type.tumor_mutant : Categorized_data_frame([], tumor_mutant_num)}
	for line in in_f:
		categorized_psi[Sample_type.tumor_mutant].data.append(float(line.split()[1]))
		categorized_psi[Sample_type.tumor_wild_type].data.append(float(line.split()[2]))
		categorized_psi[Sample_type.norma].data.append(float(line.split()[3]))
	in_f.close()
	return categorized_psi

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
	n, bins, patches = pylab.hist(data, bins=nbins, normed=0, histtype='bar', color=colors, label=labels)
	pylab.title(plot_label)
	pylab.legend(loc='best')
	fig_hist.savefig(pic_path)
	pylab.close(fig_hist)

def get_difference_data(categorized_psi):
	data = {}
	if not categorized_psi.has_key(Sample_type.norma):
		return None
	if categorized_psi.has_key(Sample_type.tumor_wild_type):
		tumor_wild_type_norma = [categorized_psi[Sample_type.tumor_wild_type].data[i] - categorized_psi[Sample_type.norma].data[i] for i in xrange(len(categorized_psi[Sample_type.norma].data)) if (~numpy.isnan(categorized_psi[Sample_type.norma].data[i]) and ~numpy.isnan(categorized_psi[Sample_type.tumor_wild_type].data[i]) and ~numpy.isnan(categorized_psi[Sample_type.tumor_mutant].data[i]))]
		if len(tumor_wild_type_norma) > 0:
			data[Difference_type.tumor_wild_type_minus_norma] = Data_frame('wt - norma', 'blue', tumor_wild_type_norma)
	if categorized_psi.has_key(Sample_type.tumor_mutant):
		tumor_mutant_norma = [categorized_psi[Sample_type.tumor_mutant].data[i] - categorized_psi[Sample_type.norma].data[i] for i in xrange(len(categorized_psi[Sample_type.norma].data)) if (~numpy.isnan(categorized_psi[Sample_type.norma].data[i]) and ~numpy.isnan(categorized_psi[Sample_type.tumor_wild_type].data[i]) and ~numpy.isnan(categorized_psi[Sample_type.tumor_mutant].data[i]))]
		if len(tumor_mutant_norma) > 0:
			data[Difference_type.tumor_mutant_minus_norma] = Data_frame('mut - norma', 'red', tumor_mutant_norma)
	return data

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print 'Usage:', sys.argv[0], '-d <data directory> -p <pictures directory>  -m <mutant gene name>'
		exit()

	parser = argparse.ArgumentParser(prog = sys.argv[0], description='Build histogram of average PSI data and delta PSI tumor and norma')
	parser.add_argument('-d', '--data_dir', help='data directory', required=True)
	parser.add_argument('-p', '--pics_dir', help='pictures directory', required=True)
	parser.add_argument('-m', '--mut_gene', help='mutatn gene name', required=True)
	args = parser.parse_args()
	data_dir = args.data_dir
	pics_dir = args.pics_dir
	mutant_gene = args.mut_gene.strip()

	hist_dir = os.path.join(pics_dir, 'histogram_of_PSI_for_' + mutant_gene)
	delta_hist_dir = os.path.join(pics_dir, 'histogram_of_delta_PSI_tumor_norm_for_' + mutant_gene)

	if not os.path.exists(hist_dir):
		os.mkdir(hist_dir)
	if not os.path.exists(delta_hist_dir):
		os.mkdir(delta_hist_dir)

	data_dir_list = [os.path.join(data_dir, d) for d in os.listdir(data_dir)]
	for d in data_dir_list:
		print 'processing', os.path.basename(d)
		in_fn = os.path.join(d, os.path.basename(d) + '_PSI_averaged_by_' + mutant_gene + '.txt')
		categorized_psi = read_psi_average_data(in_fn)

		df = [el for el in [Data_frame('Tumor_wild_type, ' + str(categorized_psi[Sample_type.tumor_wild_type].samples_num) + ' samples', 'blue', categorized_psi[Sample_type.tumor_wild_type].data), Data_frame('Normal, ' + str(categorized_psi[Sample_type.norma].samples_num) + ' samples', 'chartreuse', categorized_psi[Sample_type.norma].data), Data_frame('Tumor_mutant, ' + str(categorized_psi[Sample_type.tumor_mutant].samples_num) + ' samples', 'red', categorized_psi[Sample_type.tumor_mutant].data)] if el.num > 0]
		filter_data(df)
		hist_path = os.path.join(hist_dir, os.path.basename(d) + '_hist.png')
		draw_hist(df, hist_path, 'PSI for ' + os.path.basename(d), 15)

		diff_df = get_difference_data(categorized_psi)
		if len(diff_df) > 0:
			delta_path = os.path.join(delta_hist_dir, os.path.basename(d) + '_delta_hist.png')
			draw_hist(diff_df.values(), delta_path, 'Delta PSI for ' + os.path.basename(d), 25)

"""		U, pval = stats.mannwhitneyu(tumor_setd2_broken_filtered, tumor_filtered)
		pval *= 2 # two-sided hypothesis
		print 'mann-whitneyu for PSI pval', pval
		if pval < 0.05:
			print('Reject NULL hypothesis - Significant differences exist between groups.')
		if pval > 0.05:
			print('Accept NULL hypothesis - No significant difference between groups.')

		if len(tumor_setd2_normal) == 0 or len(tumor_normal) == 0:
			print
			continue

		U, pval = stats.mannwhitneyu(tumor_setd2_normal, tumor_normal)
		pval *= 2 # two-sided hypothesis
		print 'mann-whitneyu for delta PSI pval', pval
		if pval < 0.05:
			print('Reject NULL hypothesis - Significant differences exist between groups.')
		if pval > 0.05:
			print('Accept NULL hypothesis - No significant difference between groups.')


		T, pval = stats.ranksums([abs(x) for x in tumor_setd2_normal], [abs(x) for x in tumor_normal])
		print 'ranksums for abs(delta PSI) pval', pval
		if pval < 0.05:
			print('Reject NULL hypothesis - Significant differences exist between groups.')
		if pval > 0.05:
			print('Accept NULL hypothesis - No significant difference between groups.')
		
		print
"""

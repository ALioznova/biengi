#!/usr/bin/env python
# -*- coding: utf-8 -*- 
import os
import sys
import argparse
import numpy
import pylab
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats
from scipy.stats import mstats
from enum import Enum

class Sample_type(Enum):
	norma = 0
	tumor_wild_type = 1
	tumor_mutant = 2

class Categorized_data_frame:
	def __init__(self, data, samples_num):
		self.data = data
		self.samples_num = samples_num

class Data_frame_xy:
	def __init__(self, label, color, (x, y)):
		self.label = label
		self.color = color
		self.x = x
		self.y = y
		self.num = len(x)

	def delete_elems_by_index(self, will_del_index):
		x_new = []
		y_new = []
		for i in xrange(self.num):
			if not will_del_index[i]:
				x_new.append(self.x[i])
				y_new.append(self.y[i])
		self.x = x_new
		self.y = y_new
		self.num = len(x_new)

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

def read_expr_dist_data(in_fn):
	in_f = open(in_fn)
	header_line = in_f.readline()
	dist_bp = []
	dist_perc = []
	categorized_expr = {Sample_type.norma : Categorized_data_frame([], None), Sample_type.tumor_wild_type : Categorized_data_frame([], None), Sample_type.tumor_mutant : Categorized_data_frame([], None)}
	for line in in_f:
		dist_bp.append([int(el) for el in line.split('\t')[1].split(',') if (el != '' and el != 'nan')])
		dist_perc.append([float(el) for el in line.split('\t')[2].split(',') if (el != '' and el != 'nan')])
		categorized_expr[Sample_type.tumor_mutant].data.append([float(el) for el in line.split('\t')[3].split(',') if el != '' and el != 'nan'])
		categorized_expr[Sample_type.tumor_wild_type].data.append([float(el) for el in line.split('\t')[4].split(',') if el != '' and el != 'nan'])
		categorized_expr[Sample_type.norma].data.append([float(el) for el in line.split('\t')[5].strip().split(',') if el != '' and el != 'nan'])
	in_f.close()
	return (dist_bp, dist_perc, categorized_expr)

def prepare_data_for_plot(x_arr_of_arr, y_arr):
	only_nan = True
	for elem in y_arr:
		if ~numpy.isnan(elem):
			only_nan = False
			break
	if only_nan:
		return ([],[])
	x = []
	y = []
	for i in xrange(len(y_arr)):
		if len(x_arr_of_arr[i]) != 1: # ignore exons corresponding to multiple genes
			continue
		x.append(x_arr_of_arr[i][0])
		y.append(y_arr[i])
	return (x, y) 

def prepare_data_for_plot_delta_dist(x_arr_of_arr, y_arr1, y_arr2):
	only_nan1 = True
	for elem in y_arr1:
		if ~numpy.isnan(elem):
			only_nan1 = False
			break
	only_nan2 = True
	for elem in y_arr2:
		if ~numpy.isnan(elem):
			only_nan2 = False
			break
	if only_nan1 or only_nan2:
		return ([],[])
	x = []
	y = []
	for i in xrange(len(y_arr1)):
		if len(x_arr_of_arr[i]) != 1: # ignore exons corresponding to multiple genes
			continue
		x.append(x_arr_of_arr[i][0])
		y.append(abs(y_arr1[i] - y_arr2[i]))
	return (x, y) 

def prepare_data_for_plot_delta_expr(x_arr_of_arr1, x_arr_of_arr2, x_arr_of_arr3, y_arr1, y_arr2, y_arr3):
	only_nan1 = True
	for elem in y_arr1:
		if ~numpy.isnan(elem):
			only_nan1 = False
			break
	only_nan2 = True
	for elem in y_arr2:
		if ~numpy.isnan(elem):
			only_nan2 = False
			break
	only_nan3 = True
	for elem in y_arr3:
		if ~numpy.isnan(elem):
			only_nan3 = False
			break
	if only_nan1 or only_nan2 or only_nan3:
		return ([],[])
	x = []
	y = []
	for i in xrange(len(y_arr1)):
		if len(x_arr_of_arr1[i]) != 1 or len(x_arr_of_arr2[i]) != 1  or len(x_arr_of_arr3[i]) != 1: # ignore exons corresponding to multiple genes
			continue
		if numpy.isnan(y_arr1[i]) or numpy.isnan(y_arr2[i]) or numpy.isnan(y_arr3[i]):
			continue
#		x.append(float(x_arr_of_arr1[i][0] - x_arr_of_arr2[i][0]))
		x.append(float(x_arr_of_arr1[i][0] + x_arr_of_arr2[i][0])/2.0)
		y.append(y_arr1[i] - y_arr2[i])
	return (x, y) 

def plot_points(fig_path, plot_label, xscale, xy_data_frames):
	will_del_index = []
	for i in xrange(xy_data_frames[0].num):
		will_del = False
		for df in xy_data_frames:
			if numpy.isnan(df.y[i]):
				will_del = True
		will_del_index.append(will_del)
	for df in xy_data_frames:
		df.delete_elems_by_index(will_del_index)

	fig = plt.figure()
	plt.xscale(xscale)
	for df in xy_data_frames:
		plt.plot(df.x, df.y, marker='.', color=df.color, ls='', label=df.label, alpha=.5, ms=6)
	plt.title(plot_label)
	plt.legend(loc='best')
	fig.savefig(fig_path)
	plt.close(fig)

def plot_bins(fig_path, plot_label, xscale, nbins_x, nbins_y, y_lim, xy_data_frames):
	def find_max_num_for_plot(bins_y, xy_data_frames, y_lim):
		max_num = 0
		for df in xy_data_frames:
			for elem in df.y_split:
				y_bins = numpy.linspace(y_lim[0], y_lim[1], num=nbins_y+1, endpoint=True)
				freq, _ = numpy.histogram(elem, bins = y_bins)
				if max(freq) > max_num:
					max_num = max(freq)
		return max_num

	will_del_index = []
	for i in xrange(xy_data_frames[0].num):
		will_del = False
		for df in xy_data_frames:
			if numpy.isnan(df.x[i]):
				will_del = True
			if numpy.isnan(df.y[i]):
				will_del = True
		will_del_index.append(will_del)
	for df in xy_data_frames:
		df.delete_elems_by_index(will_del_index)

	x_data = []
	for df in xy_data_frames:
		x_data = numpy.concatenate([x_data, df.x])
	bins_x = stats.mstats.mquantiles(x_data, numpy.linspace(0.0, 1.0, num=nbins_x+1, endpoint=True))
	freq, _ = numpy.histogram(x_data, bins = bins_x)
	for df in xy_data_frames:
		df.ind = numpy.digitize(df.x, bins_x)
		y_split = []
		for i in xrange(1, nbins_x + 1):
			y_split.append([df.y[j] for j in xrange(df.num) if df.ind[j] == i])
#			y_split.append([el for el in ([df.y[j] for j in xrange(df.num) if df.ind[j] == i]) if ~numpy.isnan(el)])
		df.y_split = y_split
	max_num = find_max_num_for_plot(nbins_y, xy_data_frames, y_lim)
	max_num *= 1.1

	plot_label = plot_label + ', x_max=' + str(max_num)
	fig_bin, axes = plt.subplots(nrows=1, ncols=nbins_x, sharey=True, figsize=(24,6))
	x_bin_num = 0
	for ax in axes.flat:
		y_data = []
		colors = []
		labels = []
		for df in xy_data_frames:
			y_data.append(df.y_split[x_bin_num])
			colors.append(df.color)
			labels.append(df.label)
		y_bins = numpy.linspace(y_lim[0], y_lim[1], num=nbins_y+1, endpoint=True)
		ax.hist(y_data, bins=y_bins, normed=0, histtype='bar', color=colors, label=labels, orientation="horizontal")
		ax.set_xscale(xscale)
		ax.set_xlim((0, max_num))
		ax.set_ylim(y_lim)
		ax.get_xaxis().set_ticks([])
		ax.set_xlabel("{0:.1f}-{1:.1f}".format(bins_x[x_bin_num], bins_x[x_bin_num+1]))
		x_bin_num += 1
	fig_bin.suptitle(plot_label, fontsize=14)
	fig_bin.subplots_adjust(bottom=0.25)
	patches = []
	for df in xy_data_frames:
		patches.append(mpatches.Patch(color=df.color, label=df.label))
	plt.figlegend(handles=patches, labels=labels, loc = (0.1, 0.05))
	fig_bin.savefig(fig_path, dpi=100)
	plt.close(fig_bin)


if __name__ == '__main__':
	if len(sys.argv) == 1:
		print 'Usage:', sys.argv[0], '-d <data directory> -p <pictures directory>  -m <mutant gene name>'
		exit()

	parser = argparse.ArgumentParser(prog = sys.argv[0], description='Build plots for psi vs distance to gene start (in bp and percents) and for psi vs gene expression, data is split to tumor_mutant, tumor_wild_type and norma')
	parser.add_argument('-d', '--data_dir', help='data directory', required=True)
	parser.add_argument('-p', '--pics_dir', help='pictures directory', required=True)
	parser.add_argument('-m', '--mut_gene', help='mutatn gene name', required=True)
	args = parser.parse_args()
	data_dir = args.data_dir
	pics_dir = args.pics_dir
	mutant_gene = args.mut_gene.strip()

	expr_dir = os.path.join(pics_dir, 'expression_vs_PSI_for_' + mutant_gene)
	expr_average_dir = os.path.join(pics_dir, 'average_expr_vs_delta_PSI_for_' + mutant_gene)
	dist_bp_dir = os.path.join(pics_dir, 'distance_bp_vs_PSI_for_' + mutant_gene)
	dist_bp_delta_dir = os.path.join(pics_dir, 'distance_bp_vs_delta_PSI_for_' + mutant_gene)
	dist_perc_dir = os.path.join(pics_dir, 'distance_perc_vs_PSI_for_' + mutant_gene)
	dist_perc_delta_dir = os.path.join(pics_dir, 'distance_perc_vs_delta_PSI_for_' + mutant_gene)

	if not os.path.exists(expr_dir):
		os.mkdir(expr_dir)
	if not os.path.exists(expr_average_dir):
		os.mkdir(expr_average_dir)
	if not os.path.exists(dist_bp_dir):
		os.mkdir(dist_bp_dir)
	if not os.path.exists(dist_bp_delta_dir):
		os.mkdir(dist_bp_delta_dir)
	if not os.path.exists(dist_perc_dir):
		os.mkdir(dist_perc_dir)
	if not os.path.exists(dist_perc_delta_dir):
		os.mkdir(dist_perc_delta_dir)


	data_dir_list = [os.path.join(data_dir, d) for d in os.listdir(data_dir)]
	for d in data_dir_list:
		print 'processing', os.path.basename(d)
		psi_fn = os.path.join(d, os.path.basename(d) + '_PSI_averaged_by_' + mutant_gene + '.txt')
		categorized_psi = read_psi_average_data(psi_fn)
		expr_dist_fn = os.path.join(d, os.path.basename(d) + '_expression_and_dist_split_by_' + mutant_gene + '.txt')
		(dist_bp, dist_perc, categorized_expr) = read_expr_dist_data(expr_dist_fn)

		nbins_x = 10
		nbins_y = 15
		y_lim = (0.0, 1.0)

		# expression
		df = [el for el in [Data_frame_xy('Tumor', 'blue', prepare_data_for_plot(categorized_expr[Sample_type.tumor_wild_type].data, categorized_psi[Sample_type.tumor_wild_type].data)), Data_frame_xy('Normal', 'chartreuse', prepare_data_for_plot(categorized_expr[Sample_type.norma].data, categorized_psi[Sample_type.norma].data)), Data_frame_xy('Tumor_setd2', 'red', prepare_data_for_plot(categorized_expr[Sample_type.tumor_mutant].data, categorized_psi[Sample_type.tumor_mutant].data))] if el.num > 0]
		if len(df) > 0:
			expr_path = os.path.join(expr_dir, os.path.basename(d) + '_psi_vs_expr.png')
			plot_points(expr_path, 'PSI vs expression, ' + os.path.basename(d), 'log', df)
			expr_bin_path = os.path.join(expr_dir, os.path.basename(d) + '_psi_vs_expr_bin.png')
			plot_bins(expr_bin_path, 'PSI vs expression for ' + os.path.basename(d), 'linear', nbins_x, nbins_y, y_lim, df)

		# average expr, delpa PSI
		df = [el for el in [Data_frame_xy('tumor-norm', 'blue', prepare_data_for_plot_delta_expr(categorized_expr[Sample_type.tumor_wild_type].data, categorized_expr[Sample_type.norma].data, categorized_expr[Sample_type.tumor_mutant].data, categorized_psi[Sample_type.tumor_wild_type].data, categorized_psi[Sample_type.norma].data, categorized_psi[Sample_type.tumor_mutant].data)), Data_frame_xy('tumor_setd2-norm', 'red', prepare_data_for_plot_delta_expr(categorized_expr[Sample_type.tumor_mutant].data, categorized_expr[Sample_type.norma].data, categorized_expr[Sample_type.tumor_wild_type].data, categorized_psi[Sample_type.tumor_mutant].data, categorized_psi[Sample_type.norma].data, categorized_psi[Sample_type.tumor_wild_type].data))] if el.num > 0]
		if len(df) > 0:
			dist_perc_delta_path = os.path.join(expr_average_dir, os.path.basename(d) + '_delta_psi_vs_av_expression.png')
			plot_points(dist_perc_delta_path, 'delta PSI vs average expression, ' + os.path.basename(d), 'linear', df)
			dist_perc_delta_bin_path = os.path.join(expr_average_dir, os.path.basename(d) + '_delta_psi_vs_av_expression_bin.png')
			plot_bins(dist_perc_delta_bin_path, 'delta PSI vs average expression for ' + os.path.basename(d), 'linear', nbins_x, nbins_y, (-1.0, 1.0), df)

		# distance in bp
		df = [el for el in [Data_frame_xy('Tumor', 'blue', prepare_data_for_plot(dist_bp, categorized_psi[Sample_type.tumor_wild_type].data)), Data_frame_xy('Normal', 'chartreuse', prepare_data_for_plot(dist_bp, categorized_psi[Sample_type.norma].data)), Data_frame_xy('Tumor_setd2', 'red', prepare_data_for_plot(dist_bp, categorized_psi[Sample_type.tumor_mutant].data))] if el.num > 0]
		if len(df) > 0:
			dist_bp_path = os.path.join(dist_bp_dir, os.path.basename(d) + '_psi_vs_dist_bp.png')
			plot_points(dist_bp_path, 'PSI vs distance in bp, ' + os.path.basename(d), 'log', df)
			dist_bp_bin_path = os.path.join(dist_bp_dir, os.path.basename(d) + '_psi_vs_dist_bp_bin.png')
			plot_bins(dist_bp_bin_path, 'PSI vs distance in bp for ' + os.path.basename(d), 'linear', nbins_x, nbins_y, y_lim, df)

		# distance in percent
		df = [el for el in [Data_frame_xy('Tumor', 'blue', prepare_data_for_plot(dist_perc, categorized_psi[Sample_type.tumor_wild_type].data)), Data_frame_xy('Normal', 'chartreuse', prepare_data_for_plot(dist_perc, categorized_psi[Sample_type.norma].data)), Data_frame_xy('Tumor_setd2', 'red', prepare_data_for_plot(dist_perc, categorized_psi[Sample_type.tumor_mutant].data))] if el.num > 0]
		if len(df) > 0:
			dist_perc_path = os.path.join(dist_perc_dir, os.path.basename(d) + '_psi_vs_dist_perc.png')
			plot_points(dist_perc_path, 'PSI vs distance in perc, ' + os.path.basename(d), 'linear', df)
			dist_perc_bin_path = os.path.join(dist_perc_dir, os.path.basename(d) + '_psi_vs_dist_perc_bin.png')
			plot_bins(dist_perc_bin_path, 'PSI vs distance in perc for ' + os.path.basename(d), 'linear', nbins_x, nbins_y, y_lim, df)

		# delta distance in bp, delpa PSI
		df = [el for el in [Data_frame_xy('abs(tumor-norm)', 'blue', prepare_data_for_plot_delta_dist(dist_bp, categorized_psi[Sample_type.tumor_wild_type].data, categorized_psi[Sample_type.norma].data)), Data_frame_xy('abs(tumor_setd2-norm)', 'red', prepare_data_for_plot_delta_dist(dist_bp, categorized_psi[Sample_type.tumor_mutant].data, categorized_psi[Sample_type.norma].data))] if el.num > 0]
		if len(df) > 0:
			dist_bp_delta_path = os.path.join(dist_bp_delta_dir, os.path.basename(d) + '_delta_psi_vs_dist_bp.png')
			plot_points(dist_bp_delta_path, 'delta PSI vs distance in bp, ' + os.path.basename(d), 'log', df)
			dist_bp_delta_bin_path = os.path.join(dist_bp_delta_dir, os.path.basename(d) + '_delta_psi_vs_dist_bp_bin.png')
			plot_bins(dist_bp_delta_bin_path, 'delta PSI vs distance in bp for ' + os.path.basename(d), 'linear', nbins_x, nbins_y, y_lim, df)

		# delta distance in perc, delpa PSI
		df = [el for el in [Data_frame_xy('abs(tumor-norm)', 'blue', prepare_data_for_plot_delta_dist(dist_perc, categorized_psi[Sample_type.tumor_wild_type].data, categorized_psi[Sample_type.norma].data)), Data_frame_xy('abs(tumor_setd2-norm)', 'red', prepare_data_for_plot_delta_dist(dist_perc, categorized_psi[Sample_type.tumor_mutant].data, categorized_psi[Sample_type.norma].data))] if el.num > 0]
		if len(df) > 0:
			dist_perc_delta_path = os.path.join(dist_perc_delta_dir, os.path.basename(d) + '_delta_psi_vs_dist_perc.png')
			plot_points(dist_perc_delta_path, 'delta PSI vs distance in perc, ' + os.path.basename(d), 'linear', df)
			dist_perc_delta_bin_path = os.path.join(dist_perc_delta_dir, os.path.basename(d) + '_delta_psi_vs_dist_perc_bin.png')
			plot_bins(dist_perc_delta_bin_path, 'delta PSI vs distance in perc for ' + os.path.basename(d), 'linear', nbins_x, nbins_y, y_lim, df)



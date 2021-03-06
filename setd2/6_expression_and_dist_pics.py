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
from scipy.stats import chi2

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
		self.y_split = None
		self.x_split = None

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
		x.append(numpy.mean([x_arr_of_arr1[i][0], x_arr_of_arr2[i][0], x_arr_of_arr3[i][0]]))
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
		if len(df.x) == 0:
			return

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
		if len(df.x) == 0:
			return

	x_data = []
	for df in xy_data_frames:
		x_data = numpy.concatenate([x_data, df.x])
	bins_x = stats.mstats.mquantiles(x_data, numpy.linspace(0.0, 1.0, num=nbins_x+1, endpoint=True))
	freq, _ = numpy.histogram(x_data, bins = bins_x)
	for df in xy_data_frames:
		df.ind = numpy.digitize(df.x, bins_x)
		y_split = []
		x_split = []
		for i in xrange(1, nbins_x + 1):
			if len([df.y[j] for j in xrange(df.num) if df.ind[j] == i]) == 0:
				return
			y_split.append([df.y[j] for j in xrange(df.num) if df.ind[j] == i])
			x_split.append([df.x[j] for j in xrange(df.num) if df.ind[j] == i])
		df.y_split = y_split
		df.x_split = x_split
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

def plot_std_and_ci_for_expr(data_frames, pic_path, plot_title, alpha):
	def ci_for_std(arr, alpha):
		# http://stats.stackexchange.com/questions/76444/why-is-chi-square-used-when-creating-a-confidence-interval-for-the-variance
		df = len(arr) - 1
		ci_lower = numpy.sqrt((df*((numpy.std(arr))**2))/(chi2.ppf(1.0 - alpha/2.0, df)))
		ci_upper = numpy.sqrt((df*((numpy.std(arr))**2))/(chi2.ppf(alpha/2.0, df)))
		return (ci_lower, ci_upper)
	class Std_data:
		def __init__(self, label, color, x, y, y_err_lower, y_err_upper):
			self.label = label
			self.color = color
			self.x = x
			self.y = y
			self.y_err_lower = y_err_lower
			self.y_err_upper = y_err_upper
	std_data = []
	for elem in df:
		std_data.append(Std_data(elem.label, elem.color, [], [], [], []))
		for bin_num in xrange(len(elem.y_split)):
			arr = elem.y_split[bin_num]
			std_data[-1].x.append(numpy.mean(elem.x_split[bin_num]))
			std_val = numpy.std(arr)
			std_data[-1].y.append(std_val)
			(ci_lower, ci_upper) = ci_for_std(arr, alpha)
			std_data[-1].y_err_lower.append(std_val - ci_lower)
			std_data[-1].y_err_upper.append(ci_upper - std_val)
		std_data[-1].asymmetric_error = [std_data[-1].y_err_lower, std_data[-1].y_err_upper]

	fig, ax = plt.subplots(nrows=1, sharex=True)
	for elem in std_data:
		ax.errorbar(elem.x, elem.y, yerr=elem.asymmetric_error, fmt='-o', color=elem.color, label=elem.label, capsize=7, elinewidth=2)
		ax.fill_between(elem.x, numpy.array(elem.y)-numpy.array(elem.asymmetric_error[0]), numpy.array(elem.y)+numpy.array(elem.asymmetric_error[1]), alpha=0.5, facecolor=elem.color)
	ax.legend(loc='best')
	ax.set_title(plot_title)
	ax.set_xscale('log')
	fig.savefig(pic_path)
	plt.close(fig)


if __name__ == '__main__':
	if len(sys.argv) == 1:
		print 'Usage:', sys.argv[0], '-c <computations directory> -p <pictures directory>  -m <mutant gene name>'
		exit()

	parser = argparse.ArgumentParser(prog = sys.argv[0], description='Build plots for psi vs distance to gene start (in bp and percents) and for psi vs gene expression, data is split to tumor_mutant, tumor_wild_type and norma')
	parser.add_argument('-c', '--comp_dir', help='computations directory', required=True)
	parser.add_argument('-p', '--pics_dir', help='pictures directory', required=True)
	parser.add_argument('-m', '--mut_gene', help='mutatn gene name', required=True)
	args = parser.parse_args()
	comp_dir = args.comp_dir
	pics_dir = args.pics_dir
	mutant_gene = args.mut_gene

	expr_dir = os.path.join(pics_dir, 'gene_expression_vs_PSI_for_' + mutant_gene)
	expr_average_dir = os.path.join(pics_dir, 'gene_expression_averaged_vs_delta_PSI_for_' + mutant_gene)
	exon_expr_dir = os.path.join(pics_dir, 'exon_expression_vs_PSI_for_' + mutant_gene)
	exon_expr_average_dir = os.path.join(pics_dir, 'exon_expression_averaged_vs_delta_PSI_for_' + mutant_gene)
	dist_bp_dir = os.path.join(pics_dir, 'distance_bp_vs_PSI_for_' + mutant_gene)
	dist_bp_delta_dir = os.path.join(pics_dir, 'distance_bp_vs_delta_PSI_for_' + mutant_gene)
	dist_perc_dir = os.path.join(pics_dir, 'distance_perc_vs_PSI_for_' + mutant_gene)
	dist_perc_delta_dir = os.path.join(pics_dir, 'distance_perc_vs_delta_PSI_for_' + mutant_gene)

	if not os.path.exists(expr_dir):
		os.mkdir(expr_dir)
	if not os.path.exists(expr_average_dir):
		os.mkdir(expr_average_dir)
	if not os.path.exists(exon_expr_dir):
		os.mkdir(exon_expr_dir)
	if not os.path.exists(exon_expr_average_dir):
		os.mkdir(exon_expr_average_dir)
	if not os.path.exists(dist_bp_dir):
		os.mkdir(dist_bp_dir)
	if not os.path.exists(dist_bp_delta_dir):
		os.mkdir(dist_bp_delta_dir)
	if not os.path.exists(dist_perc_dir):
		os.mkdir(dist_perc_dir)
	if not os.path.exists(dist_perc_delta_dir):
		os.mkdir(dist_perc_delta_dir)

	data_dir_list = [os.path.join(comp_dir, d) for d in os.listdir(comp_dir)]
	for d in data_dir_list:
		print 'processing', os.path.basename(d)
		psi_fn = os.path.join(os.path.join(d, mutant_gene), os.path.basename(d) + '_PSI_averaged_by_' + mutant_gene + '.txt')
		categorized_psi = read_psi_average_data(psi_fn)
		expr_dist_fn = os.path.join(os.path.join(d, mutant_gene), os.path.basename(d) + '_expression_and_dist_split_by_' + mutant_gene + '.txt')
		(dist_bp, dist_perc, categorized_expr) = read_expr_dist_data(expr_dist_fn)

		exon_expr_fn = os.path.join(os.path.join(d, mutant_gene), os.path.basename(d) + '_exon_expresion_split_by_' + mutant_gene + '.txt')
		if not os.path.isfile(exon_expr_fn):
			exon_expr = None
		else:
			exon_expr = read_exon_expr_data(exon_expr_fn)

		nbins_x = 10
		nbins_y = 15
		y_lim = (0.0, 1.0)

		# gene expression
		df = [el for el in [Data_frame_xy('Tumor_wt', 'blue', prepare_data_for_plot(categorized_expr[Sample_type.tumor_wild_type].data, categorized_psi[Sample_type.tumor_wild_type].data)), Data_frame_xy('Normal', 'chartreuse', prepare_data_for_plot(categorized_expr[Sample_type.norma].data, categorized_psi[Sample_type.norma].data)), Data_frame_xy('Tumor_mutant', 'red', prepare_data_for_plot(categorized_expr[Sample_type.tumor_mutant].data, categorized_psi[Sample_type.tumor_mutant].data))] if el.num > 0]
		if len(df) > 0:
			expr_path = os.path.join(expr_dir, os.path.basename(d) + '_psi_vs_expr.png')
			plot_points(expr_path, mutant_gene + '. PSI vs expression, ' + os.path.basename(d), 'log', df)
			expr_bin_path = os.path.join(expr_dir, os.path.basename(d) + '_psi_vs_expr_bin.png')
			plot_bins(expr_bin_path, mutant_gene + '. PSI vs expression for ' + os.path.basename(d), 'linear', nbins_x, nbins_y, y_lim, df)

		# average gene expr, delpa PSI
		df = [el for el in [Data_frame_xy('tumor_wt-norm', 'blue', prepare_data_for_plot_delta_expr(categorized_expr[Sample_type.tumor_wild_type].data, categorized_expr[Sample_type.norma].data, categorized_expr[Sample_type.tumor_mutant].data, categorized_psi[Sample_type.tumor_wild_type].data, categorized_psi[Sample_type.norma].data, categorized_psi[Sample_type.tumor_mutant].data)), Data_frame_xy('tumor_mut-norm', 'red', prepare_data_for_plot_delta_expr(categorized_expr[Sample_type.tumor_mutant].data, categorized_expr[Sample_type.norma].data, categorized_expr[Sample_type.tumor_wild_type].data, categorized_psi[Sample_type.tumor_mutant].data, categorized_psi[Sample_type.norma].data, categorized_psi[Sample_type.tumor_wild_type].data))] if el.num > 0]
		if len(df) > 0:
			expr_delta_path = os.path.join(expr_average_dir, os.path.basename(d) + '_delta_psi_vs_av_expression.png')
			plot_points(expr_delta_path, mutant_gene + '. delta PSI vs average expression, ' + os.path.basename(d), 'linear', df)
			expr_delta_bin_path = os.path.join(expr_average_dir, os.path.basename(d) + '_delta_psi_vs_av_expression_bin.png')
			plot_bins(expr_delta_bin_path, mutant_gene + '. delta PSI vs average expression for ' + os.path.basename(d), 'linear', nbins_x, nbins_y, (-1.0, 1.0), df)
			expr_std_delta_path = os.path.join(expr_average_dir, os.path.basename(d) + '_std_for_delta_psi_vs_average_expression.png')
			plot_std_and_ci_for_expr(df, expr_std_delta_path, mutant_gene + '. std for delta PSI vs average expression for ' + os.path.basename(d), 0.05)

		if exon_expr:
			# exon expression
			df = [el for el in [Data_frame_xy('Tumor_wt', 'blue', prepare_data_for_plot(exon_expr[Sample_type.tumor_wild_type].data, categorized_psi[Sample_type.tumor_wild_type].data)), Data_frame_xy('Normal', 'chartreuse', prepare_data_for_plot(exon_expr[Sample_type.norma].data, categorized_psi[Sample_type.norma].data)), Data_frame_xy('Tumor_mutant', 'red', prepare_data_for_plot(exon_expr[Sample_type.tumor_mutant].data, categorized_psi[Sample_type.tumor_mutant].data))] if el.num > 0]
			if len(df) > 0:
				exon_expr_path = os.path.join(exon_expr_dir, os.path.basename(d) + '_psi_vs_exon_expr.png')
				plot_points(exon_expr_path, mutant_gene + '. PSI vs exon expression, ' + os.path.basename(d), 'log', df)
				exon_expr_bin_path = os.path.join(exon_expr_dir, os.path.basename(d) + '_psi_vs_exon_expr_bin.png')
				plot_bins(exon_expr_bin_path, mutant_gene + '. PSI vs exon expression for ' + os.path.basename(d), 'linear', nbins_x, nbins_y, y_lim, df)

			# average gene expr, delpa PSI
			df = [el for el in [Data_frame_xy('tumor_wt-norm', 'blue', prepare_data_for_plot_delta_expr(exon_expr[Sample_type.tumor_wild_type].data, exon_expr[Sample_type.norma].data, exon_expr[Sample_type.tumor_mutant].data, categorized_psi[Sample_type.tumor_wild_type].data, categorized_psi[Sample_type.norma].data, categorized_psi[Sample_type.tumor_mutant].data)), Data_frame_xy('tumor_mut-norm', 'red', prepare_data_for_plot_delta_expr(exon_expr[Sample_type.tumor_mutant].data, exon_expr[Sample_type.norma].data, exon_expr[Sample_type.tumor_wild_type].data, categorized_psi[Sample_type.tumor_mutant].data, categorized_psi[Sample_type.norma].data, categorized_psi[Sample_type.tumor_wild_type].data))] if el.num > 0]
			if len(df) > 0:
				exon_expr_delta_path = os.path.join(exon_expr_average_dir, os.path.basename(d) + '_delta_psi_vs_av_exon_expression.png')
				plot_points(exon_expr_delta_path, mutant_gene + '. delta PSI vs average exon expression, ' + os.path.basename(d), 'linear', df)
				exon_expr_delta_bin_path = os.path.join(exon_expr_average_dir, os.path.basename(d) + '_delta_psi_vs_av_exon_expression_bin.png')
				plot_bins(exon_expr_delta_bin_path, mutant_gene + '. delta PSI vs average exon expression for ' + os.path.basename(d), 'linear', nbins_x, nbins_y, (-1.0, 1.0), df)
				exon_expr_std_delta_path = os.path.join(exon_expr_average_dir, os.path.basename(d) + '_std_for_delta_psi_vs_average_exon_expression.png')
				plot_std_and_ci_for_expr(df, exon_expr_std_delta_path, mutant_gene + '. std for delta PSI vs average exon expression for ' + os.path.basename(d), 0.05)

		# distance in bp
		df = [el for el in [Data_frame_xy('Tumor_wt', 'blue', prepare_data_for_plot(dist_bp, categorized_psi[Sample_type.tumor_wild_type].data)), Data_frame_xy('Normal', 'chartreuse', prepare_data_for_plot(dist_bp, categorized_psi[Sample_type.norma].data)), Data_frame_xy('Tumor_mut', 'red', prepare_data_for_plot(dist_bp, categorized_psi[Sample_type.tumor_mutant].data))] if el.num > 0]
		if len(df) > 0:
			dist_bp_path = os.path.join(dist_bp_dir, os.path.basename(d) + '_psi_vs_dist_bp.png')
			plot_points(dist_bp_path, mutant_gene + '. PSI vs distance in bp, ' + os.path.basename(d), 'log', df)
			dist_bp_bin_path = os.path.join(dist_bp_dir, os.path.basename(d) + '_psi_vs_dist_bp_bin.png')
			plot_bins(dist_bp_bin_path, mutant_gene + '. PSI vs distance in bp for ' + os.path.basename(d), 'linear', nbins_x, nbins_y, y_lim, df)

		# distance in percent
		df = [el for el in [Data_frame_xy('Tumor_wt', 'blue', prepare_data_for_plot(dist_perc, categorized_psi[Sample_type.tumor_wild_type].data)), Data_frame_xy('Normal', 'chartreuse', prepare_data_for_plot(dist_perc, categorized_psi[Sample_type.norma].data)), Data_frame_xy('Tumor_mut', 'red', prepare_data_for_plot(dist_perc, categorized_psi[Sample_type.tumor_mutant].data))] if el.num > 0]
		if len(df) > 0:
			dist_perc_path = os.path.join(dist_perc_dir, os.path.basename(d) + '_psi_vs_dist_perc.png')
			plot_points(dist_perc_path, mutant_gene + '. PSI vs distance in perc, ' + os.path.basename(d), 'linear', df)
			dist_perc_bin_path = os.path.join(dist_perc_dir, os.path.basename(d) + '_psi_vs_dist_perc_bin.png')
			plot_bins(dist_perc_bin_path, mutant_gene + '. PSI vs distance in perc for ' + os.path.basename(d), 'linear', nbins_x, nbins_y, y_lim, df)

		# delta distance in bp, delpa PSI
		df = [el for el in [Data_frame_xy('abs(tumor_wt-norm)', 'blue', prepare_data_for_plot_delta_dist(dist_bp, categorized_psi[Sample_type.tumor_wild_type].data, categorized_psi[Sample_type.norma].data)), Data_frame_xy('abs(tumor_mut-norm)', 'red', prepare_data_for_plot_delta_dist(dist_bp, categorized_psi[Sample_type.tumor_mutant].data, categorized_psi[Sample_type.norma].data))] if el.num > 0]
		if len(df) > 0:
			dist_bp_delta_path = os.path.join(dist_bp_delta_dir, os.path.basename(d) + '_delta_psi_vs_dist_bp.png')
			plot_points(dist_bp_delta_path, mutant_gene + '. delta PSI vs distance in bp, ' + os.path.basename(d), 'log', df)
			dist_bp_delta_bin_path = os.path.join(dist_bp_delta_dir, os.path.basename(d) + '_delta_psi_vs_dist_bp_bin.png')
			plot_bins(dist_bp_delta_bin_path, mutant_gene + '. delta PSI vs distance in bp for ' + os.path.basename(d), 'linear', nbins_x, nbins_y, y_lim, df)

		# delta distance in perc, delpa PSI
		df = [el for el in [Data_frame_xy('abs(tumor_wt-norm)', 'blue', prepare_data_for_plot_delta_dist(dist_perc, categorized_psi[Sample_type.tumor_wild_type].data, categorized_psi[Sample_type.norma].data)), Data_frame_xy('abs(tumor_mut-norm)', 'red', prepare_data_for_plot_delta_dist(dist_perc, categorized_psi[Sample_type.tumor_mutant].data, categorized_psi[Sample_type.norma].data))] if el.num > 0]
		if len(df) > 0:
			dist_perc_delta_path = os.path.join(dist_perc_delta_dir, os.path.basename(d) + '_delta_psi_vs_dist_perc.png')
			plot_points(dist_perc_delta_path, mutant_gene + '. delta PSI vs distance in perc, ' + os.path.basename(d), 'linear', df)
			dist_perc_delta_bin_path = os.path.join(dist_perc_delta_dir, os.path.basename(d) + '_delta_psi_vs_dist_perc_bin.png')
			plot_bins(dist_perc_delta_bin_path, mutant_gene + '. delta PSI vs distance in perc for ' + os.path.basename(d), 'linear', nbins_x, nbins_y, y_lim, df)



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

def read_psi_data(in_fn):
	in_f = open(in_fn)
	header_line = in_f.readline()
	tumor_setd2_broken_num = int((header_line.split()[1]).split(':')[1])
	tumor_num = int((header_line.split()[2]).split(':')[1])
	normal_num = int((header_line.split()[3]).split(':')[1])
	tumor_setd2_broken = []
	tumor = []
	normal = []
	for line in in_f:
		tumor_setd2_broken.append(float(line.split()[1]))
		tumor.append(float(line.split()[2]))
		normal.append(float(line.split()[3]))
	in_f.close()
	return (tumor_setd2_broken_num, tumor_setd2_broken, tumor_num, tumor, normal_num, normal)

def read_expr_dist_data(in_fn):
	in_f = open(in_fn)
	header_line = in_f.readline()
	dist_bp = []
	dist_perc = []
	tumor_setd2 = []
	tumor = []
	normal = []
	for line in in_f:
		dist_bp.append([int(el) for el in line.split('\t')[1].split(',') if (el != '' and el != 'nan')])
		dist_perc.append([float(el) for el in line.split('\t')[2].split(',') if (el != '' and el != 'nan')])
		tumor_setd2.append([float(el) for el in line.split('\t')[3].split(',') if el != '' and el != 'nan'])
		tumor.append([float(el) for el in line.split('\t')[4].split(',') if el != '' and el != 'nan'])
		normal.append([float(el) for el in line.split('\t')[5].strip().split(',') if el != '' and el != 'nan'])
	in_f.close()
	return (dist_bp, dist_perc, tumor_setd2, tumor, normal)

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
	for i in xrange(len(y_arr)): # ignore exons corresponding to multiple genes
		if len(x_arr_of_arr[i]) != 1:
			continue
		x.append(x_arr_of_arr[i][0])
		y.append(y_arr[i])
	return (x, y) 

def prepare_data_for_plot_delta(x_arr_of_arr, y_arr1, y_arr2):
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
	for i in xrange(len(y_arr1)): # ignore exons corresponding to multiple genes
		if len(x_arr_of_arr[i]) != 1:
			continue
		x.append(x_arr_of_arr[i][0])
		y.append(abs(y_arr1[i] - y_arr2[i]))
	return (x, y) 

def plot_points(fig_path, plot_label, xscale, x1, y1, label1, color1, x2, y2, label2, color2, x3=None, y3=None, label3=None, color3=None):
	fig = plt.figure()
	plt.xscale(xscale)
	plt.plot(x1, y1, marker='.', color=color1, ls='', label=label1, alpha=.5, ms=6)
	plt.plot(x2, y2, marker='.', color=color2, ls='', label=label2, alpha=.5, ms=6)
	if x3:
		plt.plot(x3, y3, marker='.', color=color3, ls='', label=label3, alpha=.5, ms=6)
	plt.title(plot_label)
	plt.legend(loc='best')
	if len(x1) != 0 or len(x2) != 0 or (x3 and len(x3) != 0):
		fig.savefig(fig_path)
	plt.close(fig)

def plot_bins(fig_path, plot_label, nbins_x, nbins_y, y_lim, x1, y1, label1, color1, x2, y2, label2, color2, x3=[], y3=[], label3=None, color3=None):
	def find_max_num_for_plot(bins_y, y1, y2, y3=None):
		max_num = 0
		for y in y1:
			freq, _ = numpy.histogram(y, bins = bins_y)
			for elem in freq:
				if elem > max_num:
					max_num = elem
		for y in y2:
			freq, _ = numpy.histogram(y, bins = bins_y)
			for elem in freq:
				if elem > max_num:
					max_num = elem
		if y3:
			for y in y3:
				freq, _ = numpy.histogram(y, bins = bins_y)
				for elem in freq:
					if elem > max_num:
						max_num = elem
		return max_num
	
	if len(x1) == 0 and len(x2) == 0 and len(x3) == 0:
		return
	data = numpy.concatenate([x1, x2, x3])
	bins_x = stats.mstats.mquantiles(data, numpy.linspace(0.0, 1.0, num=nbins_x+1, endpoint=True))
	freq, _ = numpy.histogram(data, bins = bins_x)
	if len(x1) > 0:
		ind1 = numpy.digitize(x1, bins_x)
	y1_split = []
	for i in xrange(1, nbins_x + 1):
		y1_split.append([el for el in ([y1[j] for j in xrange(len(x1)) if ind1[j] == i]) if ~numpy.isnan(el)])
	if len(x2) > 0:
		ind2 = numpy.digitize(x2, bins_x)
	y2_split = []
	for i in xrange(1, nbins_x + 1):
		y2_split.append([el for el in ([y2[j] for j in xrange(len(x2)) if ind2[j] == i]) if ~numpy.isnan(el)])
	if len(x3) > 0:
		ind3 = numpy.digitize(x3, bins_x)
	y3_split = []
	for i in xrange(1, nbins_x + 1):
		y3_split.append([el for el in ([y3[j] for j in xrange(len(x3)) if ind3[j] == i]) if ~numpy.isnan(el)])
	if color3:
		max_num = find_max_num_for_plot(nbins_y, y1_split, y2_split, y3_split)
	else:
		max_num = find_max_num_for_plot(nbins_y, y1_split, y2_split)
	fig_bin, axes = plt.subplots(nrows=1, ncols=nbins_x, sharey=True, figsize=(24,6))
	bin_num = 0
	for ax in axes.flat:
		if color3:
			ax.hist((y1_split[bin_num], y2_split[bin_num], y3_split[bin_num]), bins=nbins_y, normed=0, histtype='bar', color=[color1, color2, color3], label=[label1, label2, label3], orientation="horizontal")[0]
			ax.set_xlim((0, max_num*1.1))
		else:
			ax.hist((y1_split[bin_num], y2_split[bin_num]), bins=nbins_y, normed=0, histtype='bar', color=[color1, color2], label=[label1, label2], orientation="horizontal")[0]
			ax.set_xlim((0, max_num*1.1))
		ax.get_xaxis().set_ticks([])
		ax.set_xlabel("{0:.1f}".format(float(bins_x[bin_num] + bins_x[bin_num+1])/2))
		ax.set_ylim(y_lim)
		bin_num += 1
	fig_bin.suptitle(plot_label, fontsize=14)
	fig_bin.subplots_adjust(bottom=0.25)
	if color3:
		patches = [mpatches.Patch(color=color1, label=label1), mpatches.Patch(color=color2, label=label2), mpatches.Patch(color=color3, label=label3)]
		plt.figlegend(handles=patches, labels=[label1, label2, label3], loc = (0.1, 0.05))
	else:
		patches = [mpatches.Patch(color=color1, label=label1), mpatches.Patch(color=color2, label=label2)]
		plt.figlegend(handles=patches, labels=[label1, label2], loc = (0.1, 0.05))
	fig_bin.savefig(fig_path, dpi=100)
	plt.close(fig_bin)


if __name__ == '__main__':
	if len(sys.argv) == 1:
		print 'Usage:', sys.argv[0], '-d <data directory> -p <pictures directory>'
		exit()

	parser = argparse.ArgumentParser(prog = sys.argv[0], description='Analyze')
	parser.add_argument('-d', '--data_dir', help='data directory', required=True)
	parser.add_argument('-p', '--pics_dir', help='pictures directory', required=True)
	args = parser.parse_args()
	data_dir = args.data_dir
	pics_dir = args.pics_dir

	if not os.path.exists(os.path.join(pics_dir, 'expression')):
		os.mkdir(os.path.join(pics_dir, 'expression'))
	if not os.path.exists(os.path.join(pics_dir, 'distance_bp')):
		os.mkdir(os.path.join(pics_dir, 'distance_bp'))
	if not os.path.exists(os.path.join(pics_dir, 'distance_bp_delta')):
		os.mkdir(os.path.join(pics_dir, 'distance_bp_delta'))
	if not os.path.exists(os.path.join(pics_dir, 'distance_perc')):
		os.mkdir(os.path.join(pics_dir, 'distance_perc'))
	if not os.path.exists(os.path.join(pics_dir, 'distance_perc_delta')):
		os.mkdir(os.path.join(pics_dir, 'distance_perc_delta'))


	data_dir_list = [os.path.join(data_dir, d) for d in os.listdir(data_dir)]
	for d in data_dir_list:
		print 'processing', os.path.basename(d)
		psi_fn = os.path.join(d, os.path.basename(d) + '__t_setd2__t__n.txt')
		(tumor_setd2_broken_num, tumor_setd2_psi, tumor_num, tumor_psi, normal_num, normal_psi) = read_psi_data(psi_fn)
		expr_dist_fn = os.path.join(d, os.path.basename(d) + '_expression_and_dist.txt')
		(dist_bp, dist_perc, tumor_setd2_expr, tumor_expr, normal_expr) = read_expr_dist_data(expr_dist_fn)


		# expression
		(x1, y1) = prepare_data_for_plot(tumor_expr, tumor_psi)
		(x2, y2) = prepare_data_for_plot(normal_expr, normal_psi)
		(x3, y3) = prepare_data_for_plot(tumor_setd2_expr, tumor_setd2_psi)
		expr_path = os.path.join(os.path.join(pics_dir, 'expression'), os.path.basename(d) + '_expr.png')
		plot_points(expr_path, 'PSI vs expression, ' + os.path.basename(d), 'log', x1, y1, 'Tumor', 'blue', x2, y2, 'Normal', 'chartreuse', x3, y3, 'Tumor_setd2', 'red')
		expr_bin_path = os.path.join(os.path.join(pics_dir, 'expression'), os.path.basename(d) + '_expr_bin.png')
		plot_bins(expr_bin_path, 'PSI vs expression for ' + os.path.basename(d), 10, 15, [0,1], x3, y3, 'Tumor_setd2', 'red', x1, y1, 'Tumor', 'blue', x2, y2, 'Normal', 'chartreuse')

		# distance in bp
		(x1, y1) = prepare_data_for_plot(dist_bp, tumor_psi)
		(x2, y2) = prepare_data_for_plot(dist_bp, normal_psi)
		(x3, y3) = prepare_data_for_plot(dist_bp, tumor_setd2_psi)
		dist_bp_path = os.path.join(os.path.join(pics_dir, 'distance_bp'), os.path.basename(d) + '_dist_bp.png')
		plot_points(dist_bp_path, 'PSI vs distance in bp, ' + os.path.basename(d), 'log', x1, y1, 'Tumor', 'blue', x2, y2, 'Normal', 'chartreuse', x3, y3, 'Tumor_setd2', 'red')
		dist_bp_bin_path = os.path.join(os.path.join(pics_dir, 'distance_bp'), os.path.basename(d) + '_dist_bp_bin.png')
		plot_bins(dist_bp_bin_path, 'PSI vs distance in bp for ' + os.path.basename(d), 10, 15, [0,1], x3, y3, 'Tumor_setd2', 'red', x1, y1, 'Tumor', 'blue', x2, y2, 'Normal', 'chartreuse')

		# distance in percent
		(x1, y1) = prepare_data_for_plot(dist_perc, tumor_psi)
		(x2, y2) = prepare_data_for_plot(dist_perc, normal_psi)
		(x3, y3) = prepare_data_for_plot(dist_perc, tumor_setd2_psi)
		dist_perc_path = os.path.join(os.path.join(pics_dir, 'distance_perc'), os.path.basename(d) + '_dist_perc.png')
		plot_points(dist_perc_path, 'PSI vs distance in perc, ' + os.path.basename(d), 'linear', x1, y1, 'Tumor', 'blue', x2, y2, 'Normal', 'chartreuse', x3, y3, 'Tumor_setd2', 'red')
		dist_perc_bin_path = os.path.join(os.path.join(pics_dir, 'distance_perc'), os.path.basename(d) + '_dist_perc_bin.png')
		plot_bins(dist_perc_bin_path, 'PSI vs distance in perc for ' + os.path.basename(d), 10, 15, [0,1], x3, y3, 'Tumor_setd2', 'red', x1, y1, 'Tumor', 'blue', x2, y2, 'Normal', 'chartreuse')

		# distance in bp, delpa PSI
		(x1, y1) = prepare_data_for_plot_delta(dist_bp, tumor_psi, normal_psi)
		(x2, y2) = prepare_data_for_plot_delta(dist_bp, tumor_setd2_psi, normal_psi)
		dist_bp_delta_path = os.path.join(os.path.join(pics_dir, 'distance_bp_delta'), os.path.basename(d) + '_dist_bp_delta.png')
		plot_points(dist_bp_delta_path, 'delta PSI vs distance in bp, ' + os.path.basename(d), 'log', x1, y1, 'abs(tumor-norm)', 'blue', x2, y2, 'abs(tumor_setd2-norm)', 'red')
		dist_bp_delta_bin_path = os.path.join(os.path.join(pics_dir, 'distance_bp_delta'), os.path.basename(d) + '_dist_bp_delta_bin.png')
		plot_bins(dist_bp_delta_bin_path, 'delta PSI vs distance in bp for ' + os.path.basename(d), 10, 15, [0,1], x1, y1, 'abs(tumor-norm)', 'blue', x2, y2, 'abs(tumor_setd2-norm)', 'red')

		# distance in perc, delpa PSI
		(x1, y1) = prepare_data_for_plot_delta(dist_perc, tumor_psi, normal_psi)
		(x2, y2) = prepare_data_for_plot_delta(dist_perc, tumor_setd2_psi, normal_psi)
		dist_perc_delta_path = os.path.join(os.path.join(pics_dir, 'distance_perc_delta'), os.path.basename(d) + '_dist_perc_delta.png')
		plot_points(dist_perc_delta_path, 'delta PSI vs distance in perc, ' + os.path.basename(d), 'linear', x1, y1, 'abs(tumor-norm)', 'blue', x2, y2, 'abs(tumor_setd2-norm)', 'red')
		dist_perc_delta_bin_path = os.path.join(os.path.join(pics_dir, 'distance_perc_delta'), os.path.basename(d) + '_dist_perc_delta_bin.png')
		plot_bins(dist_perc_delta_bin_path, 'delta PSI vs distance in perc for ' + os.path.basename(d), 10, 15, [0,1], x1, y1, 'abs(tumor-norm)', 'blue', x2, y2, 'abs(tumor_setd2-norm)', 'red')


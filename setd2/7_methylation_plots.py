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
from enum import Enum, unique
from sets import Set

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

class Sample_type(Enum):
	unknown = -1
	norma = 0
	tumor_wild_type = 1
	tumor_unclassified = 2
	tumor_mutant = 3

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

class Categorized_data_frame:
	def __init__(self, data, samples_num):
		self.data = data
		self.samples_num = samples_num

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

def read_methylation_data(methylation_fn):
	methylation_data = {}
	inf = open(methylation_fn)
	inf.readline()
	for line in inf:
		(pos, tumor_wt, tumor_mut, norma) = line.strip().split()
		tumor_wt = float(tumor_wt)
		tumor_mut = float(tumor_mut)
		norma = float(norma)
		(chr_name, begin, end) = pos.split(':')
		if chr_name.startswith('chr'):
			chr_name = chr_name[len('chr'):]
		begin = int(begin)
		end = int(end)
		if not methylation_data.has_key(chr_name):
			methylation_data[chr_name] = {}
		methylation_data[chr_name][(begin, end)] = {Sample_type.norma: norma, Sample_type.tumor_wild_type: tumor_wt, Sample_type.tumor_mutant : tumor_mut}
	inf.close()
	return methylation_data

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
		x_split = []
		for i in xrange(1, nbins_x + 1):
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

def get_meth_array(methylation_data, pos, s_type):
	meth = []
	for i in xrange(len(pos)):
		if not methylation_data.has_key(pos[i].chr_name):
			meth.append(float('nan'))
			continue
		if not methylation_data[pos[i].chr_name].has_key((pos[i].end1, pos[i].beg2)):
			meth.append(float('nan'))
			continue
		meth.append(methylation_data[pos[i].chr_name][(pos[i].end1, pos[i].beg2)][s_type])
	return meth

def get_average_meth(categorized_meth):
	def only_nan(arr):
		only_nan = True
		for elem in arr:
			if ~numpy.isnan(elem):
				only_nan = False
				break
		return only_nan
	not_nan_arr = []
	for (s_type, meth) in categorized_meth.iteritems():
		if not only_nan(categorized_meth[s_type]):
			not_nan_arr.append(categorized_meth[s_type])
	if len(numpy.array(not_nan_arr)) == 0:
		return []
	average_meth = numpy.mean(numpy.array(not_nan_arr), axis = 0)
	return average_meth

def prepare_data_for_plot(x, y):
	def only_nan(arr):
		only_nan = True
		for elem in arr:
			if ~numpy.isnan(elem):
				only_nan = False
				break
		return only_nan

	if only_nan(x) or only_nan(y):
		return ([],[])
	else:
		return (x, y)

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print 'Usage:', sys.argv[0], '-c <computations directory> -p <pictures directory> -m <mutant gene name>'
		exit()

	parser = argparse.ArgumentParser(prog = sys.argv[0], description='Build plots for psi vs distance to gene start (in bp and percents) and for psi vs gene expression, data is split to tumor_mutant, tumor_wild_type and norma')
	parser.add_argument('-c', '--comp_dir', help='computations directory', required=True)
	parser.add_argument('-p', '--pics_dir', help='pictures directory', required=True)
	parser.add_argument('-m', '--mut_gene', help='mutatn gene name', required=True)
	args = parser.parse_args()
	comp_dir = args.comp_dir
	pics_dir = args.pics_dir
	mutant_gene = args.mut_gene

	meth_vs_psi_dir = os.path.join(pics_dir, 'methylation_vs_PSI_for_' + mutant_gene)
	av_meth_vs_delta_psi_dir = os.path.join(pics_dir, 'methylation_averaged_vs_delta_PSI_for_' + mutant_gene)
	delta_meth_vs_delta_psi_dir = os.path.join(pics_dir, 'methylation_delta_vs_delta_PSI_tumor_vs_norma_for_' + mutant_gene)
	delta_meth_vs_delta_psi_tumor_dir = os.path.join(pics_dir, 'methylation_delta_vs_delta_PSI_wt_vs_mut_for_' + mutant_gene)
	meth_wt_vs_mut_dir = os.path.join(pics_dir, 'methylation_wt_vs_setd2')

	if not os.path.exists(meth_vs_psi_dir):
		os.mkdir(meth_vs_psi_dir)
	if not os.path.exists(av_meth_vs_delta_psi_dir):
		os.mkdir(av_meth_vs_delta_psi_dir)
	if not os.path.exists(delta_meth_vs_delta_psi_dir):
		os.mkdir(delta_meth_vs_delta_psi_dir)
	if not os.path.exists(delta_meth_vs_delta_psi_tumor_dir):
		os.mkdir(delta_meth_vs_delta_psi_tumor_dir)
	if not os.path.exists(meth_wt_vs_mut_dir):
		os.mkdir(meth_wt_vs_mut_dir)

	data_dir_list = [os.path.join(comp_dir, d) for d in os.listdir(comp_dir)]
	for d in data_dir_list:
		print 'processing', os.path.basename(d)
		psi_fn = os.path.join(os.path.join(d, mutant_gene), os.path.basename(d) + '_PSI_averaged_by_' + mutant_gene + '.txt')
		(pos, categorized_psi) = read_psi_average_data(psi_fn)
		methylation_fn = os.path.join(os.path.join(d, mutant_gene), os.path.basename(d) + '_exon_methylation_averaged_by_' + mutant_gene + '.txt')
		methylation_data = read_methylation_data(methylation_fn)
		categorized_meth = {Sample_type.tumor_wild_type : get_meth_array(methylation_data, pos, Sample_type.tumor_wild_type), Sample_type.tumor_mutant : get_meth_array(methylation_data, pos, Sample_type.tumor_mutant), Sample_type.norma : get_meth_array(methylation_data, pos, Sample_type.norma)}
		average_meth = get_average_meth(categorized_meth)
		if len(average_meth) == 0:
			continue
		nbins_x = 10
		nbins_y = 15

		# PSI vs meth
		df = [el for el in [Data_frame_xy('Tumor_wt', 'blue', prepare_data_for_plot(categorized_meth[Sample_type.tumor_wild_type], categorized_psi[Sample_type.tumor_wild_type].data)), Data_frame_xy('Normal', 'chartreuse', prepare_data_for_plot(categorized_meth[Sample_type.norma], categorized_psi[Sample_type.norma].data)), Data_frame_xy('Tumor_mutant', 'red', prepare_data_for_plot(categorized_meth[Sample_type.tumor_mutant], categorized_psi[Sample_type.tumor_mutant].data))] if el.num > 0]
		if len(df) != 0:
			meth_path = os.path.join(meth_vs_psi_dir, os.path.basename(d) + '_psi_vs_meth.png')
			plot_points(meth_path, mutant_gene + '. PSI vs methylation, ' + os.path.basename(d), 'linear', df)
			y_lim = (0.0, 1.0)
			meth_bin_path = os.path.join(meth_vs_psi_dir, os.path.basename(d) + '_psi_vs_meth_bin.png')
			plot_bins(meth_bin_path, mutant_gene + 'PSI vs methylation for ' + os.path.basename(d), 'linear', nbins_x, nbins_y, y_lim, df)

		# delta PSI vs average meth
		df = [el for el in [Data_frame_xy('Tumor_wt - norma', 'blue', prepare_data_for_plot(average_meth, numpy.array(categorized_psi[Sample_type.tumor_wild_type].data) - numpy.array(categorized_psi[Sample_type.norma].data))), Data_frame_xy('Tumor_mutant - norma', 'red', prepare_data_for_plot(average_meth, numpy.array(categorized_psi[Sample_type.tumor_mutant].data) - numpy.array(categorized_psi[Sample_type.norma].data)))] if el.num > 0]
		if len(df) != 0:
			meth_path = os.path.join(av_meth_vs_delta_psi_dir, os.path.basename(d) + '_delta_psi_vs_average_meth.png')
			plot_points(meth_path, mutant_gene + '. delta PSI vs average methylation, ' + os.path.basename(d), 'linear', df)
			y_lim = (-1.0, 1.0)
			meth_bin_path = os.path.join(av_meth_vs_delta_psi_dir, os.path.basename(d) + '_delta_psi_vs_average_meth_bin.png')
			plot_bins(meth_bin_path, mutant_gene + '. delta PSI vs average methylation for ' + os.path.basename(d), 'linear', nbins_x, nbins_y, y_lim, df)

		# delta PSI vs delta meth
		df = [el for el in [Data_frame_xy('Tumor_wt - norma', 'blue', prepare_data_for_plot(numpy.array(categorized_meth[Sample_type.tumor_wild_type]) - numpy.array(categorized_meth[Sample_type.norma]), numpy.array(categorized_psi[Sample_type.tumor_wild_type].data) - numpy.array(categorized_psi[Sample_type.norma].data))), Data_frame_xy('Tumor_mutant - norma', 'red', prepare_data_for_plot(numpy.array(categorized_meth[Sample_type.tumor_mutant]) - numpy.array(categorized_meth[Sample_type.norma]), numpy.array(categorized_psi[Sample_type.tumor_mutant].data) - numpy.array(categorized_psi[Sample_type.norma].data)))] if el.num > 0]
		if len(df) != 0:
			meth_path = os.path.join(delta_meth_vs_delta_psi_dir, os.path.basename(d) + '_delta_psi_vs_delta_meth.png')
			plot_points(meth_path, mutant_gene + '. delta PSI vs delta methylation, ' + os.path.basename(d), 'linear', df)
			y_lim = (-1.0, 1.0)
			meth_bin_path = os.path.join(delta_meth_vs_delta_psi_dir, os.path.basename(d) + '_delta_psi_vs_delta_meth_bin.png')
			plot_bins(meth_bin_path, mutant_gene + '. delta PSI vs delta methylation for ' + os.path.basename(d), 'linear', nbins_x, nbins_y, y_lim, df)

		# delta PSI vs delta meth wt vs setd2
		df = [el for el in [Data_frame_xy('Tumor_wt - tumor_mut', 'magenta', prepare_data_for_plot(numpy.array(categorized_meth[Sample_type.tumor_wild_type]) - numpy.array(categorized_meth[Sample_type.tumor_mutant]), numpy.array(categorized_psi[Sample_type.tumor_wild_type].data) - numpy.array(categorized_psi[Sample_type.tumor_mutant].data)))] if el.num > 0]
		if len(df) != 0:
			meth_path = os.path.join(delta_meth_vs_delta_psi_tumor_dir, os.path.basename(d) + '_delta_psi_vs_delta_meth_tumor.png')
			plot_points(meth_path,  mutant_gene + '. delta PSI vs delta methylation, wt vs mutant, ' + os.path.basename(d), 'linear', df)
			y_lim = (-1.0, 1.0)
			meth_bin_path = os.path.join(delta_meth_vs_delta_psi_tumor_dir, os.path.basename(d) + '_delta_psi_vs_delta_meth_tumor_bin.png')
			plot_bins(meth_bin_path,  mutant_gene + '. delta PSI vs delta methylation, wt vs mutant for ' + os.path.basename(d), 'linear', nbins_x, nbins_y, y_lim, df)

		# methylation 
		df = [el for el in [Data_frame_xy('wt vs mut', 'magenta', prepare_data_for_plot(categorized_meth[Sample_type.tumor_wild_type], categorized_meth[Sample_type.tumor_mutant]))] if el.num > 0]
		if len(df) != 0:
			meth_path = os.path.join(meth_wt_vs_mut_dir, os.path.basename(d) + '_meth.png')
			plot_points(meth_path, mutant_gene + '. methylation vs methylation (wt vs mut), ' + os.path.basename(d), 'linear', df)




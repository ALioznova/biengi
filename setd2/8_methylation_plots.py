#!/usr/bin/env python
# -*- coding: utf-8 -*- 
import os
import sys
import argparse
import numpy
from sets import Set
from enum import Enum, unique
import pylab
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats
from scipy.stats import mstats

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

class Sample_type(Enum):
	unknown = -1
	norma = 0
	tumor_wild_type = 1
	tumor_unclassified = 2
	tumor_setd2 = 3

def read_methylation_data(methylation_fn):
	methylation_data = {}
	inf = open(methylation_fn)
	inf.readline()
	for line in inf:
		(pos, tumor_wt, tumor_setd2, norma) = line.strip().split()
		tumor_wt = float(tumor_wt)
		tumor_setd2 = float(tumor_setd2)
		norma = float(norma)
		(chr_name, begin, end) = pos.split(':')
		if chr_name.startswith('chr'):
			chr_name = chr_name[len('chr'):]
		begin = int(begin)
		end = int(end)
		if not methylation_data.has_key(chr_name):
			methylation_data[chr_name] = {}
		methylation_data[chr_name][(begin, end)] = {Sample_type.norma: norma, Sample_type.tumor_wild_type: tumor_wt, Sample_type.tumor_setd2 : tumor_setd2}
	inf.close()
	return methylation_data

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

def plot_points(fig_path, plot_label, xscale, xy_data_frames):
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

def prepare_data_for_plot(methylation_data, pos, tumor_setd2, tumor_wt, normal):
	def only_nan(arr):
		only_nan = True
		for elem in arr:
			if ~numpy.isnan(elem):
				only_nan = False
				break
		return only_nan
	
	tumor_wt_psi = []
	tumor_wt_meth = []
	tumor_setd2_psi = []
	tumor_setd2_meth = []
	normal_psi = []
	normal_meth = []
	for i in xrange(len(pos)):
		if not methylation_data.has_key(pos[i].chr_name):
			continue
		if not methylation_data[pos[i].chr_name].has_key((pos[i].end1, pos[i].beg2)):
			continue
		tumor_wt_psi.append(tumor_wt[i])
		tumor_wt_meth.append(methylation_data[pos[i].chr_name][(pos[i].end1, pos[i].beg2)][Sample_type.tumor_wild_type])
		tumor_setd2_psi.append(tumor_setd2[i])
		tumor_setd2_meth.append(methylation_data[pos[i].chr_name][(pos[i].end1, pos[i].beg2)][Sample_type.tumor_setd2])
		normal_psi.append(normal[i])
		normal_meth.append(methylation_data[pos[i].chr_name][(pos[i].end1, pos[i].beg2)][Sample_type.norma])
	df = []
	if not only_nan(tumor_wt_psi) and not only_nan(tumor_wt_meth):
		df.append(Data_frame_xy('tumor_wt', 'blue', (tumor_wt_meth, tumor_wt_psi)))
	if not only_nan(tumor_setd2_psi) and not only_nan(tumor_setd2_meth):
		df.append(Data_frame_xy('tumor_setd2', 'red', (tumor_setd2_meth, tumor_setd2_psi)))
	if not only_nan(normal_psi) and not only_nan(normal_meth):
		df.append(Data_frame_xy('normal', 'chartreuse', (normal_meth, normal_psi)))
	return df

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print 'Usage:', sys.argv[0], '-d <data directory>'
		exit()

	parser = argparse.ArgumentParser(prog = sys.argv[0], description='Analyze')
	parser.add_argument('-d', '--data_dir', help='data directory', required=True)
	parser.add_argument('-p', '--pics_dir', help='pictures directory', required=True)
	args = parser.parse_args()
	data_dir = args.data_dir
	pics_dir = args.pics_dir

	if not os.path.exists(os.path.join(pics_dir, 'methylation')):
		os.mkdir(os.path.join(pics_dir, 'methylation'))
	
	data_dir_list = [os.path.join(data_dir, d) for d in os.listdir(data_dir)]
	for d in data_dir_list:
		print 'processing', os.path.basename(d)
		in_fn = os.path.join(d, os.path.basename(d) + '_PSI_average.txt')
		(pos, tumor_setd2_broken_num, tumor_setd2_broken, tumor_num, tumor, normal_num, normal) = read_psi_average_data(in_fn)
		methylation_fn = os.path.join(d, os.path.basename(d) + '_exon_methylation.txt')
		methylation_data = read_methylation_data(methylation_fn)
		df = prepare_data_for_plot(methylation_data, pos, tumor_setd2_broken, tumor, normal)
		if len(df) != 0:

			meth_path = os.path.join(os.path.join(pics_dir, 'methylation'), os.path.basename(d) + '_meth.png')
			plot_points(meth_path, 'PSI vs methylation, ' + os.path.basename(d), 'linear', df)
			nbins_x = 10
			nbins_y = 15
			y_lim = (0.0, 1.0)
			meth_bin_path = os.path.join(os.path.join(pics_dir, 'methylation'), os.path.basename(d) + '_meth_bin.png')
			plot_bins(meth_bin_path, 'PSI vs methylation for ' + os.path.basename(d), 'linear', nbins_x, nbins_y, y_lim, df)


#!/usr/bin/env python
# -*- coding: utf-8 -*- 
import os
import sys
import argparse
import numpy
import pylab
from scipy import stats
from scipy.stats import mstats

def read_data(in_fn):
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

def draw_hist(data_frames, pic_path, plot_label):
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
	data = []
	colors = []
	labels = []
	for df in data_frames:
		data.append(df.y)
		colors.append(df.color)
		labels.append(df.label)
	fig_hist = pylab.figure()
	n, bins, patches = pylab.hist(data, bins=15, normed=0, histtype='bar', color=colors, label=labels)
	pylab.title(plot_label)
	pylab.legend(loc='best')
	fig_hist.savefig(pic_path)
	pylab.close(fig_hist)
	tumor = []
	tumor_setd2 = []
	normal = []
	for df in data_frames:
		if 'Tumor_setd2,' in df.label:
			tumor_setd2 = df.y
		if 'Tumor,' in df.label:
			tumor = df.y
		if 'Normal,' in df.label:
			normal = df.y
	return (tumor_setd2, tumor, normal)

def draw_delta_hist(tumor_setd2_broken, tumor, normal, pic_path):
	tumor_normal = []
	if len(tumor) != 0:
		for i in xrange(len(normal)):
			if  ~numpy.isnan(tumor[i]) and ~numpy.isnan(normal[i]):
				tumor_normal.append(tumor[i] - normal[i])
	tumor_setd2_normal = []
	if len(tumor_setd2_broken) != 0:
		for i in xrange(len(normal)):
			if  ~numpy.isnan(tumor_setd2_broken[i]) and ~numpy.isnan(normal[i]):
				tumor_setd2_normal.append(tumor_setd2_broken[i] - normal[i])
	if len(tumor_normal) == 0 and len(tumor_setd2_normal) == 0:
		return ([], [])
	data_normal = [tumor_setd2_normal, tumor_normal]
	fig_hist_delta = pylab.figure()
	n, bins, patches = pylab.hist(data_normal, bins=35, normed=0, histtype='bar', color=['crimson', 'blue'], label=['tsetd2-n', 't-n'])
	pylab.title('Delta PSI for ' + os.path.basename(d))
	pylab.legend(loc='best')
	fig_hist_delta.savefig(pic_path)
	pylab.close(fig_hist_delta)
	return (tumor_setd2_normal, tumor_normal)

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

	if not os.path.exists(os.path.join(pics_dir, 'hist')):
		os.mkdir(os.path.join(pics_dir, 'hist'))
	if not os.path.exists(os.path.join(pics_dir, 'delta_hist')):
		os.mkdir(os.path.join(pics_dir, 'delta_hist'))


	data_dir_list = [os.path.join(data_dir, d) for d in os.listdir(data_dir)]
	for d in data_dir_list:
		print 'processing', os.path.basename(d)
		in_fn = os.path.join(d, os.path.basename(d) + '_PSI_average.txt')
		(tumor_setd2_broken_num, tumor_setd2_broken, tumor_num, tumor, normal_num, normal) = read_data(in_fn)

		df = [el for el in [Data_frame('Tumor, ' + str(tumor_num) + ' samples', 'blue', tumor), Data_frame('Normal, ' + str(normal_num) + ' samples', 'chartreuse', normal), Data_frame('Tumor_setd2, ' + str(tumor_setd2_broken_num) + ' samples', 'red', tumor_setd2_broken)] if el.num > 0]
		hist_path = os.path.join(os.path.join(pics_dir, 'hist'), os.path.basename(d) + '_hist.png')
		(tumor_setd2_broken_filtered, tumor_filtered, normal_filtered) = draw_hist(df, hist_path, 'PSI for ' + os.path.basename(d))

		delta_path = os.path.join(os.path.join(pics_dir, 'delta_hist'), os.path.basename(d) + '_delta_hist.png')
		(tumor_setd2_normal, tumor_normal) = draw_delta_hist(tumor_setd2_broken_filtered, tumor_filtered, normal_filtered, delta_path)


		U, pval = stats.mannwhitneyu(tumor_setd2_broken_filtered, tumor_filtered)
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


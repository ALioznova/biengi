#!/usr/bin/env python
# -*- coding: utf-8 -*- 
import os
import sys
import argparse
import numpy
import pylab
from scipy import stats
from scipy.stats import mstats
import random

def read_data(in_fn):
	in_f = open(in_fn)
	header_line = in_f.readline()
	tumor_setd2_broken_num = int((header_line.split()[0]).split(':')[1])
	tumor_num = int((header_line.split()[1]).split(':')[1])
	normal_num = int((header_line.split()[2]).split(':')[1])
	tumor_setd2_broken = []
	tumor = []
	normal = []
	for line in in_f:
		tumor_setd2_broken.append(float(line.split()[0]))
		tumor.append(float(line.split()[1]))
		normal.append(float(line.split()[2]))
	in_f.close()
	return (tumor_setd2_broken_num, tumor_setd2_broken, tumor_num, tumor, normal_num, normal)

def draw_hist(tumor_setd2_broken_num, tumor_setd2_broken, tumor_num, tumor, normal_num, normal, pic_path):
	tumor_setd2_broken_filtered = [t for t in tumor_setd2_broken if ~numpy.isnan(t)]
	tumor_filtered = [n for n in tumor if ~numpy.isnan(n)]
	normal_filtered = [n for n in normal if ~numpy.isnan(n)]
	data = [tumor_setd2_broken_filtered, tumor_filtered, normal_filtered]

	fig_hist = pylab.figure()
	n, bins, patches = pylab.hist(data, bins=15, normed=1, histtype='bar', color=['crimson', 'blue', 'chartreuse'], label=['Tumor setd2 mutation impact high, ' + str(tumor_setd2_broken_num) + ' samples', 'Tumor no setd2 mutation, ' + str(tumor_num) + ' samples', 'Normal, ' + str(normal_num) + ' samples'])
	pylab.title('PSI for ' + os.path.basename(d))
	pylab.legend(loc='best')
	fig_hist.savefig(pic_path)
	pylab.close(fig_hist)
	return (tumor_setd2_broken_filtered, tumor_filtered, normal_filtered)

def draw_delta_hist(tumor_setd2_broken, tumor, normal, pic_path):
	tumor_normal = []
	for i in xrange(len(normal)):
		if  ~numpy.isnan(tumor[i]) and ~numpy.isnan(normal[i]):
			tumor_normal.append(tumor[i] - normal[i])
	tumor_setd2_normal = []
	for i in xrange(len(normal)):
		if  ~numpy.isnan(tumor_setd2_broken[i]) and ~numpy.isnan(normal[i]):
			tumor_setd2_normal.append(tumor_setd2_broken[i] - normal[i])
	if len(tumor_normal) == 0 and len(tumor_setd2_normal) == 0:
		return ([], [])
	data_normal = [tumor_setd2_normal, tumor_normal]
	fig_hist_delta = pylab.figure()
	n, bins, patches = pylab.hist(data_normal, bins=35, normed=1, histtype='bar', color=['crimson', 'blue'], label=['tsetd2-n', 't-n'])
	pylab.title('Delta PSI for ' + os.path.basename(d))
	pylab.legend(loc='best')
	fig_hist_delta.savefig(pic_path)
	pylab.close(fig_hist_delta)
	return (tumor_setd2_normal, tumor_normal)

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print "Usage:", sys.argv[0], "-d <data directory> -p <pictures directory>"
		exit()

	parser = argparse.ArgumentParser(prog = sys.argv[0], description='Analyze')
	parser.add_argument("-d", "--data_dir", help="data directory", required=True)
	parser.add_argument("-p", "--pics_dir", help="pictures directory", required=True)
	args = parser.parse_args()
	data_dir = args.data_dir
	pics_dir = args.pics_dir

	if not os.path.exists(os.path.join(pics_dir, 'hist')):
		os.mkdir(os.path.join(pics_dir, 'hist'))
	if not os.path.exists(os.path.join(pics_dir, 'delta_hist')):
		os.mkdir(os.path.join(pics_dir, 'delta_hist'))


	data_dir_list = [os.path.join(data_dir, d) for d in os.listdir(data_dir)]
	for d in data_dir_list:
		print "processing", os.path.basename(d)
		in_fn = os.path.join(d, os.path.basename(d) + '__t_setd2__t__n.txt')
		(tumor_setd2_broken_num, tumor_setd2_broken, tumor_num, tumor, normal_num, normal) = read_data(in_fn)
		hist_path = os.path.join(os.path.join(pics_dir, 'hist'), os.path.basename(d) + '_hist.png')
		(tumor_setd2_broken_filtered, tumor_filtered, normal_filtered) = draw_hist(tumor_setd2_broken_num, tumor_setd2_broken, tumor_num, tumor, normal_num, normal, hist_path)
		delta_path = os.path.join(os.path.join(pics_dir, 'delta_hist'), os.path.basename(d) + '_delta_hist.png')
		(tumor_setd2_normal, tumor_normal) = draw_delta_hist(tumor_setd2_broken, tumor, normal, delta_path)

		U, pval = stats.mannwhitneyu(tumor_setd2_broken_filtered, tumor_filtered)
		pval *= 2 # two-sided hypothesis
		print 'mann-whitneyu PSI pval', pval
		if pval < 0.05:
			print("Reject NULL hypothesis - Significant differences exist between groups.")
		if pval > 0.05:
			print("Accept NULL hypothesis - No significant difference between groups.")

		if len(tumor_setd2_normal) == 0 or len(tumor_normal) == 0:
			continue

		U, pval = stats.mannwhitneyu(tumor_setd2_normal, tumor_normal)
		pval *= 2 # two-sided hypothesis
		print 'mann-whitneyu delta PSI pval', pval
		if pval < 0.05:
			print("Reject NULL hypothesis - Significant differences exist between groups.")
		if pval > 0.05:
			print("Accept NULL hypothesis - No significant difference between groups.")


		T, pval = stats.ranksums([abs(x) for x in tumor_setd2_normal], [abs(x) for x in tumor_normal])
		print 'ranksums delta PSI pval', pval
		if pval < 0.05:
			print("Reject NULL hypothesis - Significant differences exist between groups.")
		if pval > 0.05:
			print("Accept NULL hypothesis - No significant difference between groups.")


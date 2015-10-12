#!/usr/bin/env python
# -*- coding: utf-8 -*- 
import os
import sys
import argparse
import numpy
import pylab

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

	data_dir_list = [os.path.join(data_dir, d) for d in os.listdir(data_dir)]
	for d in data_dir_list:
		print "processing", d
		tumor_setd2_broken_fn = os.path.join(d, os.path.basename(d) + '_tumor_setd2_broken_filtered.txt')
		tumor_fn = os.path.join(d, os.path.basename(d) + '_tumor_filtered.txt')
		normal_fn = os.path.join(d, os.path.basename(d) + '_normal_filtered.txt')
		
		tumor_setd2_broken_f = open(tumor_setd2_broken_fn)
		tumor_setd2_broken_num = int(tumor_setd2_broken_f.readline().split(':')[1])
		tumor_setd2_broken = []
		for line in tumor_setd2_broken_f:
			tumor_setd2_broken.append(float(line))
		tumor_setd2_broken_f.close()
		
		tumor_f = open(tumor_fn)
		tumor_num = int(tumor_f.readline().split(':')[1])
		tumor = []
		for line in tumor_f:
			tumor.append(float(line))
		tumor_f.close()
		
		normal_f = open(normal_fn)
		normal_num = int(normal_f.readline().split(':')[1])
		normal = []
		for line in normal_f:
			normal.append(float(line))
		normal_f.close()


		tumor_setd2_broken_filtered = [t for t in tumor_setd2_broken if ~numpy.isnan(t)]
		tumor_filtered = [n for n in tumor if ~numpy.isnan(n)]
		normal_filtered = [n for n in normal if ~numpy.isnan(n)]

		data = [tumor_setd2_broken_filtered, tumor_filtered, normal_filtered]

		fig_hist = pylab.figure()
		n, bins, patches = pylab.hist(data, bins=15, normed=1, histtype='bar', color=['crimson', 'blue', 'chartreuse'], label=['Tumor setd2 mutation impact high, ' + str(tumor_setd2_broken_num) + ' samples', 'Tumor no setd2 mutation, ' + str(tumor_num) + ' samples', 'Normal, ' + str(normal_num) + ' samples'])
		pylab.title('PSI for ' + os.path.basename(d))
		pylab.legend(loc='best')
		fig_hist.savefig(os.path.join(d, os.path.basename(d) + '_tsetd2_t_n_filtered.png'))
		fig_hist.savefig(pics_dir + 't_setd2_t_n_filtered/' + os.path.basename(d) + '_tsetd2_t_n_filtered.png')

		tumor_normal = []
		for i in xrange(len(normal)):
			if  ~numpy.isnan(tumor[i]) and ~numpy.isnan(normal[i]):
				tumor_normal.append(tumor[i] - normal[i])

		tumor_setd2_normal = []
		for i in xrange(len(normal)):
			if  ~numpy.isnan(tumor_setd2_broken[i]) and ~numpy.isnan(normal[i]):
				tumor_setd2_normal.append(tumor_setd2_broken[i] - normal[i])

		if len(tumor_normal) == 0 and len(tumor_setd2_normal) == 0:
			continue

		data_normal = [tumor_setd2_normal, tumor_normal]

		fig_hist_normal = pylab.figure()
		n, bins, patches = pylab.hist(data_normal, bins=35, normed=1, histtype='bar', color=['crimson', 'blue'], label=['tsetd2-n', 't-n'])
		pylab.title('Delta PSI for ' + os.path.basename(d))
		pylab.legend(loc='best')
		fig_hist_normal.savefig(os.path.join(d, os.path.basename(d) + '_minus_normal_filtered.png'))
		fig_hist_normal.savefig(pics_dir + 'minus_normal_filtered/' + os.path.basename(d) + '_minus_normal_filtered.png')



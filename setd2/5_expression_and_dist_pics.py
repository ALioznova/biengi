#!/usr/bin/env python
# -*- coding: utf-8 -*- 
import os
import sys
import argparse
import numpy
import pylab
import matplotlib.pyplot as plt
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
	dist = []
	tumor_setd2 = []
	tumor = []
	normal = []
	for line in in_f:
		dist.append([int(el) for el in line.split('\t')[1].split(',') if (el != '' and el != 'nan')])
		tumor_setd2.append([float(el) for el in line.split('\t')[2].split(',') if el != '' and el != 'nan'])
		tumor.append([float(el) for el in line.split('\t')[3].split(',') if el != '' and el != 'nan'])
		normal.append([float(el) for el in line.split('\t')[4].strip().split(',') if el != '' and el != 'nan'])
	in_f.close()
	return (dist, tumor_setd2, tumor, normal)

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
		if len(x_arr_of_arr[i]) != 1:
			continue
		x.append(x_arr_of_arr[i][0])
		y.append(y_arr[i])
	return (x, y) 

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
	if not os.path.exists(os.path.join(pics_dir, 'distance')):
		os.mkdir(os.path.join(pics_dir, 'distance'))


	data_dir_list = [os.path.join(data_dir, d) for d in os.listdir(data_dir)]
	for d in data_dir_list:
		print 'processing', os.path.basename(d)
		psi_fn = os.path.join(d, os.path.basename(d) + '__t_setd2__t__n.txt')
		(tumor_setd2_broken_num, tumor_setd2_psi, tumor_num, tumor_psi, normal_num, normal_psi) = read_psi_data(psi_fn)
		expr_dist_fn = os.path.join(d, os.path.basename(d) + '_expression_and_dist.txt')
		(dist, tumor_setd2_expr, tumor_expr, normal_expr) = read_expr_dist_data(expr_dist_fn)

		expr_path = os.path.join(os.path.join(pics_dir, 'expression'), os.path.basename(d) + '_expr.png')
		dist_path = os.path.join(os.path.join(pics_dir, 'distance'), os.path.basename(d) + '_dist.png')

		fig_dist = plt.figure()
		plt.xscale('log')
		(x1, y1) = prepare_data_for_plot(dist, tumor_setd2_psi)
		plt.plot(x1, y1, marker='.', color='r', ls='', label = 'tumor_setd2')
		(x2, y2) = prepare_data_for_plot(dist, tumor_psi)
		plt.plot(x2, y2, marker='.', color='b', ls='', label = 'tumor')
		(x3, y3) = prepare_data_for_plot(dist, normal_psi)
		plt.plot(x3, y3, marker='.', color='g', ls='', label = 'normal')
		plt.title('PSI vs distance, ' + os.path.basename(d))
		plt.legend(loc='best')
		fig_dist.savefig(dist_path)
		plt.close(fig_dist)

		fig_expr = plt.figure()
		plt.xscale('log')
		(x1, y1) = prepare_data_for_plot(tumor_setd2_expr, tumor_setd2_psi)
		plt.plot(x1, y1, marker='.', color='r', ls='', label = 'tumor_setd2')
		(x2, y2) = prepare_data_for_plot(tumor_expr, tumor_psi)
		plt.plot(x2, y2, marker='.', color='b', ls='', label = 'tumor')
		(x3, y3) = prepare_data_for_plot(normal_expr, normal_psi)
		plt.plot(x3, y3, marker='.', color='g', ls='', label = 'normal')
		plt.title('PSI vs expression, ' + os.path.basename(d))
		plt.legend(loc='best')
		fig_expr.savefig(expr_path)
		plt.close(fig_expr)



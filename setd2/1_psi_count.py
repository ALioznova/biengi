#!/usr/bin/env python
# -*- coding: utf-8 -*- 
import os
import sys
import argparse
from intervaltree import Interval, IntervalTree
import numpy

class data_frame:
	def __init__(self, begin, end, coverage):
		self.begin = begin
		self.end = end
		self.coverage = coverage

def read_data(in_f):
	samples_names = in_f.readline()
	samples_number = len(in_f.readline().split()) - 1
	data = {}
	for line in in_f:
		(junction_beg, junction_end) = (line.split()[0]).split(',')
		(beg_chr, beg_coord, beg_strand) = junction_beg.split(':')
		beg_coord = int(beg_coord)
		(end_chr, end_coord, end_strand) = junction_end.split(':')
		end_coord = int(end_coord)
		if beg_chr != end_chr:
			print 'CHR'
			continue
		if beg_strand != end_strand:
			print 'STRAND'
			continue
		coverage = [int(n) for n in line.split()[1:]]
		if not data.has_key(beg_chr + beg_strand):
			data[beg_chr + beg_strand] = []
		data[beg_chr + beg_strand].append(data_frame(beg_coord, end_coord, coverage))
	return (data, samples_number, samples_names)

def build_interval_trees(data):
	interval_trees = {}
	for (chr_name, chr_data) in data.iteritems():
		chr_interval_tree = IntervalTree()
		for i in xrange(len(chr_data)):
			chr_interval_tree[chr_data[i].begin : chr_data[i].end] = i
		interval_trees[chr_name] = chr_interval_tree
	return interval_trees

class three_intervals:
	def __init__(self, ab_index, bc_index, ac_index):
		self.ab = ab_index
		self.bc = bc_index
		self.ac = ac_index

def find_spliced_intervals(data):
	interval_trees = build_interval_trees(data)
	intervals = {}
	for (chr_name, chr_data) in data.iteritems():
		chr_tree = interval_trees[chr_name] 
		chr_intervals = []
		for i in xrange(len(chr_data)): # treat current interval as a long one
			left_set_all = chr_tree[chr_data[i].begin + 0.1]
			left_intervals = []
			for elem in left_set_all:
				if elem.begin == chr_data[i].begin and elem.end < chr_data[i].end:
					left_intervals.append(elem)
			right_set_all = chr_tree[chr_data[i].end - 0.1]
			right_intervals = []
			for elem in right_set_all:
				if elem.end == chr_data[i].end and elem.begin > chr_data[i].begin:
					right_intervals.append(elem)
			for left in left_intervals:
				for right in right_intervals:
					if left.end < right.begin:
						chr_intervals.append(three_intervals(left.data, right.data, i))
		intervals[chr_name] = chr_intervals
	return intervals

def calculate_psi(N_ab, N_bc, N_ac):
	# http://geuvadiswiki.crg.es/index.php/Percentage_Splicing_Index
	if ((N_ab + N_bc) * 0.5 <= 5) or (N_ac <= 5):
		return float('nan')
	return (float(N_ab + N_bc)/(N_ab + N_ac + N_bc))

def compute(in_fn, out_fn):
	in_f = open(in_fn)
	(data, samples_num, samples_names) = read_data(in_f)
	in_f.close()

	intervals = find_spliced_intervals(data)

	out_f = open(out_fn, 'w')
	out_f.write(samples_names)
	for (chr_name, chr_intervals) in intervals.iteritems():
		chr_data = data[chr_name]
		for elem in chr_intervals:
			new_line = ''
			new_line += (chr_name[:-1] + ':' + str(chr_data[elem.ab].begin) + ':' + str(chr_data[elem.ab].end) + ':' + str(chr_data[elem.bc].begin) + ':' + str(chr_data[elem.bc].end) + ':' + chr_name[-1])
			only_nan = True
			for i in xrange(samples_num):
				new_line += '\t'
				cur_psi = calculate_psi(chr_data[elem.ab].coverage[i], chr_data[elem.bc].coverage[i], chr_data[elem.ac].coverage[i])
				new_line += str(cur_psi)
				if not numpy.isnan(cur_psi):
					only_nan = False
			if not only_nan:
				out_f.write(new_line + '\n')
	out_f.close()

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print 'Usage:', sys.argv[0], '-d <data directory> -c <computations directory>'
		exit()

	parser = argparse.ArgumentParser(prog = sys.argv[0], description='PSI count')
	parser.add_argument('-d', '--data_dir', help='data directory', required=True)
	parser.add_argument('-c', '--comp_dir', help='computations directory', required=True)
	args = parser.parse_args()
	data_dir = args.data_dir
	comp_dir = args.comp_dir

	if not os.path.isdir(data_dir):
		print >> sys.stderr, 'Not a directory ' + data_dir
		sys.exit(1)

	if not os.path.exists(comp_dir):
		os.mkdir(comp_dir)

	data_dir_list = [os.path.join(data_dir, d) for d in os.listdir(data_dir)]
	for d in data_dir_list:
		print 'processing', os.path.basename(d)
		d_out = os.path.join(comp_dir, os.path.basename(d))
		if not os.path.exists(d_out):
			os.mkdir(d_out)

		for elem in os.listdir(d):
			if 'junction_quantification' in elem:
				in_f = os.path.join(d, elem)
		if not os.path.exists(in_f):
			print 'no such file:', in_f
			continue
		out_f = os.path.join(d_out, os.path.basename(d) + '_PSI.txt')
		compute(in_f, out_f)
		


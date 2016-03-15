#!/usr/bin/env python
# -*- coding: utf-8 -*- 
import os
import sys
import argparse
import numpy
import pylab
from sets import Set

def draw_hist(data, data_labels, colors, pic_path, plot_label, nbins):
	fig_hist = pylab.figure()
	n, bins, patches = pylab.hist(data, bins=nbins, normed=0, histtype='bar', color=colors, label=data_labels)
	pylab.title(plot_label)
	pylab.legend(loc='best')
	fig_hist.savefig(pic_path)
	pylab.close(fig_hist)

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print 'Usage:', sys.argv[0], '-d <data directory> -p <pictures directory>'
		exit()

	parser = argparse.ArgumentParser(prog = sys.argv[0], description='plot data split by annotation')
	parser.add_argument('-d', '--data_dir', help='data directory', required=True)
	parser.add_argument('-p', '--pic_dir', help='pictures directory', required=True)
	args = parser.parse_args()
	data_dir = args.data_dir
	pic_dir = args.pic_dir

	if not os.path.isdir(data_dir):
		print >> sys.stderr, 'Not a directory ' + data_dir
		sys.exit(1)

	if not os.path.exists(pic_dir):
		os.mkdir(pic_dir)

	files = Set(os.listdir(data_dir))
	while len(files) > 0:
		f = files.pop()
		if f.endswith('_pos.txt'):
			f_pos = f
			f_neg = f_pos[:-len('_pos.txt')] + '_neg.txt'
			files.remove(f_neg)
			annotation = f_neg[:-len('_neg.txt')]
			f_pos = os.path.join(data_dir, f_pos)
			f_neg = os.path.join(data_dir, f_neg)
		elif f.endswith('_neg.txt'):
			f_neg = f
			f_pos = f_neg[:-len('_neg.txt')] + '_pos.txt'
			files.remove(f_pos)
			annotation = f_pos[:-len('_pos.txt')]
			f_neg = os.path.join(data_dir, f_neg)
			f_pos = os.path.join(data_dir, f_pos)
		pos_corr = []
		pos_bg = []
		neg_corr = []
		neg_bg = []
		for line in open(f_pos):
			pos_corr.append(float(line.split()[0]))
			pos_bg.append(float(line.split()[1]))
		for line in open(f_neg):
			neg_corr.append(float(line.split()[0]))
			neg_bg.append(float(line.split()[1]))
		draw_hist([pos_corr, pos_bg, neg_corr, neg_bg], ['tl+', 'bg+', 'tl-', 'bg-'], ['deeppink', 'pink', 'deepskyblue', 'skyblue'], os.path.join(pic_dir, annotation + '.png'), annotation, 35)




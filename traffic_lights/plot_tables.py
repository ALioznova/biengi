#!/usr/bin/env python
# -*- coding: utf-8 -*- 
import os
import sys
import argparse
import numpy
import pylab
from sets import Set

def draw_hist(data, data_labels, colors, pic_path, plot_label, nbins):
	if sum([len(el) for el in data]) == 0:
		return
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

	parser = argparse.ArgumentParser(prog = sys.argv[0], description='plot tables')
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

	all_data = {}
	table_dirs = os.listdir(data_dir)
	for table_dir in table_dirs:
		files = Set(os.listdir(os.path.join(data_dir, table_dir)))
		while len(files) > 0:
			f = files.pop()
			if f == ('Pairs.txt'):
				continue
			if f.endswith('_tl.txt'):
				f_tl = f
				f_bg = f_tl[:-len('_tl.txt')] + '_bg.txt'
				files.remove(f_bg)
				annotation = f_tl[:-len('_tl.txt')]
				f_tl = os.path.join(os.path.join(data_dir, table_dir), f_tl)
				f_bg = os.path.join(os.path.join(data_dir, table_dir), f_bg)
			elif f.endswith('_bg.txt'):
				f_bg = f
				f_tl = f_bg[:-len('_bg.txt')] + '_tl.txt'
				files.remove(f_tl)
				annotation = f_tl[:-len('_tl.txt')]
				f_tl = os.path.join(os.path.join(data_dir, table_dir), f_tl)
				f_bg = os.path.join(os.path.join(data_dir, table_dir), f_bg)
			tl_data = []
			for line in open(f_tl):
				if line.split()[1] != 'True':
					tl_data.append(float(line.split()[1]))
				else:
					tl_data.append(1)
			bg_data = []
			for line in open(f_bg):
				if line.split()[1] != 'True':
					bg_data.append(float(line.split()[1]))
				else:
					bg_data.append(1)
			if not all_data.has_key(annotation):
				all_data[annotation] = {}
			all_data[annotation][table_dir[len('table'):]] = {'tl': tl_data, 'bg': bg_data}
	for annotation in all_data.iterkeys():
		print annotation
		data_pos_tl = []
		data_neg_tl = []
		data_whole_tl = []
		data_pos_bg = []
		data_neg_bg = []
		data_whole_bg = []
		seeds = Set(all_data[annotation].keys())
		while len(seeds) > 0:
			seed = seeds.pop()
			if seed.endswith('_pos'):
				seed_pos = seed
				seed_whole = seed_pos[:-len('_pos')]
				seeds.remove(seed_whole)
				seed_neg = seed_whole + '_neg'
				seeds.remove(seed_neg)
			elif seed.endswith('_neg'):
				seed_neg = seed
				seed_whole = seed_neg[:-len('_neg')]
				seeds.remove(seed_whole)
				seed_pos = seed_whole + '_pos'
				seeds.remove(seed_pos)
			else:
				seed_whole = seed
				seed_pos = seed + '_pos'
				seeds.remove(seed_pos)
				seed_neg = seed + '_neg'
				seeds.remove(seed_neg)
			data_pos_tl.append(all_data[annotation][seed_pos]['tl'])
			data_pos_bg.append(all_data[annotation][seed_pos]['bg'])
			data_neg_tl.append(all_data[annotation][seed_neg]['tl'])
			data_neg_bg.append(all_data[annotation][seed_neg]['bg'])
			data_whole_tl.append(all_data[annotation][seed_whole]['tl'])
			data_whole_bg.append(all_data[annotation][seed_whole]['bg'])
		draw_hist(data_pos_tl + data_pos_bg, ['tl'] * len(data_pos_tl) + ['bg'] * len(data_pos_bg), ['orangered']  * len(data_pos_tl) + ['cyan'] * len(data_pos_bg), os.path.join(pic_dir, annotation + '_pos.png'), annotation + '_pos', 11)
		draw_hist(data_neg_tl + data_neg_bg, ['tl'] * len(data_neg_tl) + ['bg'] * len(data_neg_bg), ['orangered']  * len(data_neg_tl) + ['cyan'] * len(data_neg_bg), os.path.join(pic_dir, annotation + '_neg.png'), annotation + '_neg', 11)
		draw_hist(data_whole_tl + data_whole_bg, ['tl'] * len(data_whole_tl) + ['bg'] * len(data_whole_bg), ['orangered']  * len(data_whole_tl) + ['cyan'] * len(data_whole_bg), os.path.join(pic_dir, annotation + '_whole.png'), annotation + '_whole', 11)





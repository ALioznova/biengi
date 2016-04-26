#!/usr/bin/env python
# -*- coding: utf-8 -*- 
import os
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
from sets import Set

class Data_frame:
	def __init__(self, label, color, data, global_title):
		self.label = label
		self.color = color
		self.data = data
		self.global_title = global_title

def draw_hist(data, suptitle, pic_path, nbins, normed=0):
	plt.suptitle(suptitle, fontsize=16)
	fig, a = plt.subplots(nrows=1,ncols=len(data), figsize=(24,6))
	a = a.ravel()
	for idx,ax in enumerate(a):
		sum_len = 0
		for elem in data[idx].data:
			sum_len += len(elem)
		if sum_len == 0:
			continue
		ax.hist(data[idx].data, bins=nbins, normed=normed, histtype='bar', color=data[idx].color, label=data[idx].label)
		ax.set_title(data[idx].global_title)
	plt.tight_layout()
	fig.savefig(pic_path)
	plt.close(fig)

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
				if f_bg in files:
					files.remove(f_bg)
				annotation = f_tl[:-len('_tl.txt')]
				f_tl = os.path.join(os.path.join(data_dir, table_dir), f_tl)
				f_bg = os.path.join(os.path.join(data_dir, table_dir), f_bg)
			elif f.endswith('_bg.txt'):
				f_bg = f
				f_tl = f_bg[:-len('_bg.txt')] + '_tl.txt'
				if f_tl in files:
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

	new_ann = {}
	new_ann['wgEncodeBroadHistoneGm12878H2azStdPk'] = 'H2az_Gm12878'
	new_ann['wgEncodeBroadHistoneGm12878H3k27acStdPk'] = 'H3k27ac_Gm12878'
	new_ann['wgEncodeBroadHistoneGm12878H3k27me3StdPk'] = 'H3k27me3_Gm12878_BH'
	new_ann['wgEncodeBroadHistoneGm12878H3k27me3StdPkV2'] = 'H3k27me3_Gm12878_BHV2'
	new_ann['wgEncodeBroadHistoneGm12878H3k36me3StdPk'] = 'H3k36me3_Gm12878_BH'
	new_ann['wgEncodeBroadHistoneGm12878H3k4me1StdPk'] = 'H3k4me1_Gm12878'
	new_ann['wgEncodeBroadHistoneGm12878H3k4me2StdPk'] = 'H3k4me2_Gm12878'
	new_ann['wgEncodeBroadHistoneGm12878H3k4me3StdPk'] = 'H3k4me3_Gm12878_BH'
	new_ann['wgEncodeBroadHistoneH1hescH2azStdPk'] = 'H2az_H1hesc'
	new_ann['wgEncodeBroadHistoneH1hescH3k27acStdPk'] = 'H3k27ac_H1hesc'
	new_ann['wgEncodeBroadHistoneH1hescH3k27me3StdPk'] = 'H3k27me3_H1hesc'
	new_ann['wgEncodeBroadHistoneH1hescH3k36me3StdPk'] = 'H3k36me3_H1hesc'
	new_ann['wgEncodeBroadHistoneH1hescH3k4me1StdPk'] = 'H3k4me1_H1hesc'
	new_ann['wgEncodeBroadHistoneH1hescH3k4me2StdPk'] = 'H3k4me2_H1hesc'
	new_ann['wgEncodeBroadHistoneH1hescH3k4me3StdPk'] = 'H3k4me3_H1hesc'
	new_ann['wgEncodeBroadHistoneK562H2azStdPk'] = 'H2az_K562'
	new_ann['wgEncodeBroadHistoneK562H3k27acStdPk'] = 'H3k27ac_K562'
	new_ann['wgEncodeBroadHistoneK562H3k27me3StdPk'] = 'H3k27me3_K562_BH'
	new_ann['wgEncodeBroadHistoneK562H3k36me3StdPk'] = 'H3k36me3_K562_BH'
	new_ann['wgEncodeBroadHistoneK562H3k4me1StdPk'] = 'H3k4me1_K562_BH'
	new_ann['wgEncodeBroadHistoneK562H3k4me2StdPk'] = 'H3k4me2_K562'
	new_ann['wgEncodeBroadHistoneK562H3k4me3StdPk'] = 'H3k4me3_K562_BH'
	new_ann['wgEncodeRegDnaseClustered'] = 'RegDnaseClustered'
	new_ann['wgEncodeSydhHistoneK562H3k27me3bUcdPk'] = 'H3k27me3_K562_SH'
	new_ann['wgEncodeSydhHistoneK562H3k4me1UcdPk'] = 'H3k4me1_K562_SH'
	new_ann['wgEncodeSydhHistoneK562H3k4me3bUcdPk'] = 'H3k4me3_K562_SH'
	new_ann['wgEncodeUwDgfK562Hotspots'] = 'dgf_hotspots_K562'
	new_ann['wgEncodeUwDgfK562Pk'] = 'dgf_hotspots_K562_pk'
	new_ann['wgEncodeUwHistoneGm12878H3k27me3StdHotspotsRep1'] = 'H3k27me3_Gm12878_UHHs1'
	new_ann['wgEncodeUwHistoneGm12878H3k27me3StdHotspotsRep2'] = 'H3k27me3_Gm12878_UHHs2'
	new_ann['wgEncodeUwHistoneGm12878H3k27me3StdPkRep1'] = 'H3k27me3_Gm12878_UHPk1'
	new_ann['wgEncodeUwHistoneGm12878H3k27me3StdPkRep2'] = 'H3k27me3_Gm12878_UHPk2'
	new_ann['wgEncodeUwHistoneGm12878H3k36me3StdHotspotsRep1'] = 'H3k36me3_Gm12878_UHHs1'
	new_ann['wgEncodeUwHistoneGm12878H3k36me3StdHotspotsRep2'] = 'H3k36me3_Gm12878_UHHs2'
	new_ann['wgEncodeUwHistoneGm12878H3k36me3StdPkRep1'] = 'H3k36me3_Gm12878_UHPk1'
	new_ann['wgEncodeUwHistoneGm12878H3k36me3StdPkRep2'] = 'H3k36me3_Gm12878_UHPk2'
	new_ann['wgEncodeUwHistoneGm12878H3k4me3StdHotspotsRep1'] = 'H3k4me3_Gm12878_UHHs1'
	new_ann['wgEncodeUwHistoneGm12878H3k4me3StdHotspotsRep2'] = 'H3k4me3_Gm12878_UHHs2'
	new_ann['wgEncodeUwHistoneGm12878H3k4me3StdPkRep1'] = 'H3k4me3_Gm12878_UHPk1'
	new_ann['wgEncodeUwHistoneGm12878H3k4me3StdPkRep2'] = 'H3k4me3_Gm12878_UHPk2'
	new_ann['wgEncodeUwHistoneK562H3k27me3StdHotspotsRep1'] = 'H3k27me3_K562_UHHs1'
	new_ann['wgEncodeUwHistoneK562H3k27me3StdHotspotsRep2'] = 'H3k27me3_K562_UHHs2'
	new_ann['wgEncodeUwHistoneK562H3k27me3StdPkRep1'] = 'H3k27me3_K562_UHPk1'
	new_ann['wgEncodeUwHistoneK562H3k27me3StdPkRep2'] = 'H3k27me3_K562_UHPk2'
	new_ann['wgEncodeUwHistoneK562H3k36me3StdHotspotsRep1'] = 'H3k36me3_K562_UHHs1'
	new_ann['wgEncodeUwHistoneK562H3k36me3StdHotspotsRep2'] = 'H3k36me3_K562_UHHs2'
	new_ann['wgEncodeUwHistoneK562H3k36me3StdPkRep1'] = 'H3k36me3_K562_UHPk1'
	new_ann['wgEncodeUwHistoneK562H3k36me3StdPkRep2'] = 'H3k36me3_K562_UHPk2'
	new_ann['wgEncodeUwHistoneK562H3k4me3StdHotspotsRep1'] = 'H3k4me3_K562_UHHs1'
	new_ann['wgEncodeUwHistoneK562H3k4me3StdHotspotsRep2'] = 'H3k4me3_K562_UHHs2'
	new_ann['wgEncodeUwHistoneK562H3k4me3StdPkRep1'] = 'H3k4me3_K562_UHPk1'
	new_ann['wgEncodeUwHistoneK562H3k4me3StdPkRep2'] = 'H3k4me3_K562_UHPk2'

	for annotation in all_data.iterkeys():
		if annotation != 'enhancer':
			continue
		if annotation in new_ann.keys():
			ann = new_ann[annotation]
		else:
			ann = annotation
		print ann
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

		df_pos = Data_frame(['tl'] * len(data_pos_tl) + ['bg'] * len(data_pos_bg),  ['orangered']  * len(data_pos_tl) + ['cyan'] * len(data_pos_bg), data_pos_tl + data_pos_bg, '+')
		df_neg = Data_frame(['tl'] * len(data_neg_tl) + ['bg'] * len(data_neg_bg),  ['orangered']  * len(data_neg_tl) + ['cyan'] * len(data_neg_bg), data_neg_tl + data_neg_bg, '-')
		df_whole = Data_frame(['tl'] * len(data_whole_tl) + ['bg'] * len(data_whole_bg),  ['orangered']  * len(data_whole_tl) + ['cyan'] * len(data_whole_bg), data_whole_tl + data_whole_bg, ann)

		if sum([len(set(e)) for e in data_pos_tl + data_pos_bg]) > len([len(set(e)) for e in data_pos_tl + data_pos_bg]):
			if ann == 'enhancer':
				max_enhancer = max([max([max(el) for el in df_pos.data]), max([max(el) for el in df_neg.data]), max([max(el) for el in df_whole.data])])
				enhancer_bins = np.concatenate((np.linspace(0, 300, num=4, endpoint=False), np.linspace(300, max_enhancer, num=7, endpoint=True)))
				draw_hist([df_pos, df_neg, df_whole], ann, os.path.join(pic_dir, ann + '.png'), enhancer_bins, 0)
				draw_hist([df_pos, df_neg, df_whole], ann, os.path.join(pic_dir, ann + '_normed.png'), enhancer_bins, 1)
			else:
				draw_hist([df_pos, df_neg, df_whole], ann, os.path.join(pic_dir, ann + '.png'), 11)
				draw_hist([df_pos, df_neg, df_whole], ann, os.path.join(pic_dir, ann + '_normed.png'), 11, 1)
		else:
			draw_hist([df_pos, df_neg, df_whole], ann, os.path.join(pic_dir, ann + '.png'), 2)



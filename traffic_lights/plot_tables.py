#!/usr/bin/env python
# -*- coding: utf-8 -*- 
import os
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from sets import Set
from scipy import stats

class Data_frame:
	def __init__(self, label, color, data, global_title):
		self.label = label
		self.color = color
		self.data = data
		self.global_title = global_title

def check_diff(tl, bg, total_num = 100015): # total_num = 100015 with clusters, total_num = 37391 without them
	tl_mean = np.mean(tl)
	bg_mean = np.mean(bg)
	is_different = False
	oddsratio, pvalue = stats.fisher_exact([[tl_mean, total_num - tl_mean], [bg_mean, total_num - bg_mean]])
	if pvalue < 0.00005:
		is_different = True
	return is_different

def draw_hist(df, suptitle, pic_path, nbins, normalized=0, check_significance=False):
	plt.rcParams.update({'font.size': 22})
	data = [df['cause_all']['corr_pos'], df['cause_all']['corr_neg'], df['cause_all']['corr_all'],
		df['cause_pos']['corr_pos'], df['cause_pos']['corr_neg'], df['cause_pos']['corr_all'],
		df['cause_neg']['corr_pos'], df['cause_neg']['corr_neg'], df['cause_neg']['corr_all']]
	if check_significance:
		data_worth_plotting = False
		for data_elem in data:
			bg_data = []
			tl_data = []
			cur_data = [len(d) for d in data_elem.data]
			cur_labels = data_elem.label
			for i in xrange(len(cur_labels)):
				if cur_labels[i] == 'tl':
					tl_data.append(cur_data[i])
				elif cur_labels[i] == 'bg':
					bg_data.append(cur_data[i])
				else:
					print 'unknown data label', cur_labels[i]
			is_different = check_diff(tl_data, bg_data)
			data_worth_plotting = (data_worth_plotting or is_different)
			if not is_different:
				new_colors = []
				for clr in data_elem.color:
					if clr == 'cyan':
						new_colors.append('lightcyan')
					elif clr == 'orangered':
						new_colors.append('lightsalmon')
					else:
						new_colors.append(clr)
				data_elem.color = new_colors
		if not data_worth_plotting:
			return
	fig, a = plt.subplots(nrows=3, ncols=3, figsize=(30,36))
	a = a.ravel()
	for idx,ax in enumerate(a):
		sum_len = 0
		for elem in data[idx].data:
			sum_len += len(elem)
		if sum_len == 0:
			continue
		if [sum(el) for el in data[idx].data] != [len(el) for el in data[idx].data]:
			ax.hist(data[idx].data, bins=nbins, normed=normalized, histtype='bar', color=data[idx].color, label=data[idx].label)
			continue
		data_for_plot = [len(el) for el in data[idx].data]
		bars = [np.mean(data_for_plot[:len(data_for_plot)/2]), np.mean(data_for_plot[len(data_for_plot)/2:])]
		errors = [max(data_for_plot[:len(data_for_plot)/2]) - bars[0], max(data_for_plot[len(data_for_plot)/2:]) - bars[1]]
		width = 0.25 # bar width
		ax.set_xlim(0, 1)
		ax.bar(np.array([0.2]), bars[0], width, color=data[idx].color[0]) # first param -- x-coords of the bar
		ax.bar(np.array([0.3]) + width, bars[1], width, color=data[idx].color[-1], yerr=errors[1], error_kw=dict(ecolor='gray', capsize=50, lw=5,capthick=5)) # error_kw -- error bars, caps -- horizontal part
		ax.tick_params(axis='both', which='major', labelsize=40) # y axis text size
		ax.tick_params(axis='both', which='minor', labelsize=30)
		ax.set_xticklabels([])
		ax.set_xticks([])
		ax.set_title(data[idx].global_title, fontsize=50) # title of the subplot
		ttl = ax.title
		ttl.set_position([.5, 1.05]) # h space between subplots
	plt.tight_layout()
	fig.subplots_adjust(top=0.85, hspace = 0.2) # some extra space at teh top of total pic
	plt.suptitle(suptitle.replace('_', ' '), fontsize=70) # title of the whole fig
	TL_patch = mpatches.Patch(color='orangered', label='TL')
	BG_patch = mpatches.Patch(color='cyan', label='BG')
	fig.legend( (TL_patch, BG_patch), ('TL', 'BG'), prop={'size':50} ) # legend
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

	all_data = {}
	cause_dir_all = os.path.join(data_dir, 'all')
	cause_dir_pos = os.path.join(data_dir, 'causality_pos')
	cause_dir_neg = os.path.join(data_dir, 'causality_neg')
	for cur_cause_dir in [cause_dir_all, cause_dir_pos, cause_dir_neg]:
		for table_dir in [os.path.join(cur_cause_dir, el) for el in os.listdir(cur_cause_dir)]:
			files = Set([os.path.join(table_dir, elem) for elem in os.listdir(table_dir)])
			while len(files) > 0:
				f = os.path.basename(files.pop())
				if f == ('Pairs.txt'):
					continue
				if f.endswith('_tl.txt'):
					f_tl = f
					f_bg = f_tl[:-len('_tl.txt')] + '_bg.txt'
					if f_bg in files:
						files.remove(f_bg)
					annotation = f_tl[:-len('_tl.txt')]
					f_tl = os.path.join(os.path.join(cur_cause_dir, table_dir), f_tl)
					f_bg = os.path.join(os.path.join(cur_cause_dir, table_dir), f_bg)
				elif f.endswith('_bg.txt'):
					f_bg = f
					f_tl = f_bg[:-len('_bg.txt')] + '_tl.txt'
					if f_tl in files:
						files.remove(f_tl)
					annotation = f_tl[:-len('_tl.txt')]
					f_tl = os.path.join(os.path.join(cur_cause_dir, table_dir), f_tl)
					f_bg = os.path.join(os.path.join(cur_cause_dir, table_dir), f_bg)
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
				if not all_data[annotation].has_key(os.path.basename(cur_cause_dir)):
					all_data[annotation][os.path.basename(cur_cause_dir)] = {}
				all_data[annotation][os.path.basename(cur_cause_dir)][os.path.basename(table_dir)[len('table'):]] = {'tl': tl_data, 'bg': bg_data}
	print 'Data was read'

	def get_corr_data(cur_cause, annotation, all_data):
		data_pos_tl = []
		data_neg_tl = []
		data_whole_tl = []
		data_pos_bg = []
		data_neg_bg = []
		data_whole_bg = []
		seeds = Set(all_data[annotation][cur_cause].keys())
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
			data_pos_tl.append(all_data[annotation][cur_cause][seed_pos]['tl'])
			data_pos_bg.append(all_data[annotation][cur_cause][seed_pos]['bg'])
			data_neg_tl.append(all_data[annotation][cur_cause][seed_neg]['tl'])
			data_neg_bg.append(all_data[annotation][cur_cause][seed_neg]['bg'])
			data_whole_tl.append(all_data[annotation][cur_cause][seed_whole]['tl'])
			data_whole_bg.append(all_data[annotation][cur_cause][seed_whole]['bg'])
		if (cur_cause).split('_')[-1] == 'pos':
			cause_label = 'cause > 1'
		elif (cur_cause).split('_')[-1] == 'neg':
			cause_label = 'cause < -1'
		if (cur_cause).split('_')[-1] == 'all':
			cause_label = ''
		df_pos = Data_frame(['tl'] * len(data_pos_tl) + ['bg'] * len(data_pos_bg),  ['orangered']  * len(data_pos_tl) + ['cyan'] * len(data_pos_bg), data_pos_tl + data_pos_bg, cause_label + ', corr > 0' if cause_label != '' else 'corr > 0')
		df_neg = Data_frame(['tl'] * len(data_neg_tl) + ['bg'] * len(data_neg_bg),  ['orangered']  * len(data_neg_tl) + ['cyan'] * len(data_neg_bg), data_neg_tl + data_neg_bg, cause_label + ', corr < 0' if cause_label != '' else 'corr < 0')
		df_whole = Data_frame(['tl'] * len(data_whole_tl) + ['bg'] * len(data_whole_bg),  ['orangered']  * len(data_whole_tl) + ['cyan'] * len(data_whole_bg), data_whole_tl + data_whole_bg, cause_label + '')
		data_binary = True
		if sum([len(set(e)) for e in data_pos_tl + data_pos_bg]) > len([len(set(e)) for e in data_pos_tl + data_pos_bg]):
			data_binary = False
		return (df_pos, df_neg, df_whole, data_binary)

	for annotation in all_data.iterkeys():
		if annotation in new_ann.keys():
			ann = new_ann[annotation]
		else:
			ann = annotation
		print ann
		
		df = {'cause_pos': {}, 'cause_neg': {}, 'cause_all': {}}
		(df['cause_all']['corr_pos'], df['cause_all']['corr_neg'], df['cause_all']['corr_all'], data_binary_all) = get_corr_data(os.path.basename(cause_dir_all), annotation, all_data)
		(df['cause_pos']['corr_pos'], df['cause_pos']['corr_neg'], df['cause_pos']['corr_all'], data_binary_pos) = get_corr_data(os.path.basename(cause_dir_pos), annotation, all_data)
		(df['cause_neg']['corr_pos'], df['cause_neg']['corr_neg'], df['cause_neg']['corr_all'], data_binary_neg) = get_corr_data(os.path.basename(cause_dir_neg), annotation, all_data)

		if not (data_binary_all and data_binary_pos and data_binary_neg):
			if ann == 'enhancer':
				max_enhancer = max([max([max([max(el) for el in df[cause_item][corr_item].data]) for corr_item in df[cause_item].keys()]) for cause_item in df.keys()])
				enhancer_bins = np.concatenate((np.linspace(0, 300, num=4, endpoint=False), np.linspace(300, max_enhancer, num=7, endpoint=True)))
				draw_hist(df = df, suptitle = ann, pic_path = os.path.join(pic_dir, ann + '.png'), nbins = enhancer_bins, normalized = 0)
				draw_hist(df = df, suptitle = ann + '_normalized', pic_path = os.path.join(pic_dir, ann + '_normalized.png'), nbins = enhancer_bins, normalized = 1)
				enhancer_bins = np.linspace(0, 300, num=11, endpoint=True)
				draw_hist(df = df, suptitle = ann + ' up to 300', pic_path = os.path.join(pic_dir, ann + '_cut_at_300.png'), nbins = enhancer_bins, normalized = 0)
				draw_hist(df = df, suptitle = ann + '_normalized up to 300', pic_path = os.path.join(pic_dir, ann + '_cut_at_300_normalized.png'), nbins = enhancer_bins, normalized = 1)
			elif 'CL:' in ann:
				draw_hist(df = df, suptitle = ann, pic_path = os.path.join(pic_dir, ann + '.png'), nbins = 1, check_significance = True)
			elif ann in ['enhancerNewborn', 'enhancerNeuroblastoma', 'enhancerColonCarcinoma', 'enhancerAcuteMyeloidLeukemia']:
				draw_hist(df = df, suptitle = ann + '_total', pic_path = os.path.join(pic_dir, ann + '_total' + '.png'), nbins = 1, check_significance = True)
			else:
				draw_hist(df = df, suptitle = ann, pic_path = os.path.join(pic_dir, ann + '.png'), nbins = 11)
				draw_hist(df = df, suptitle = ann + '_normalized', pic_path = os.path.join(pic_dir, ann + '_normalized.png'), nbins = 11, normalized = 1)
				draw_hist(df = df, suptitle = ann + '_total', pic_path = os.path.join(pic_dir, ann + '_total' + '.png'), nbins = 1, check_significance = True)
		else:
			draw_hist(df = df, suptitle = ann, pic_path = os.path.join(pic_dir, ann + '.png'), nbins = 1, check_significance = True)



#!/usr/bin/env python
# -*- coding: utf-8 -*- 
import os
import sys
import argparse
from Bio import SeqIO
from Bio.SeqUtils import GC
from sets import Set
import random
import bisect

class tl_record:
	def __init__(self, gc_content, cpg_content, tl_line, tl_line_num, tl_chr, tl_id, use_clusters_data=False):
		self.tl_id = tl_id
		self.gc_content = gc_content
		self.cpg_content = cpg_content
		self.line_num = tl_line_num # 1-based
		self.chr = tl_chr
		self.gene_info = tl_line.split()[0]
		self.pos = int(tl_line.split()[1])
		self.strand = tl_line.split()[2]
		self.annotation = {}
		for elem in tl_line.split()[3].split(','):
			if len(elem.split('=')) == 1:
				self.annotation[elem] = True
			else:
				self.annotation[elem.split('=')[0]] = float(elem.split('=')[1])
		calculated_data = [elem.split(')')[0].split() for elem in tl_line.split('(')[1:]]
		calculated_data = [[int(elem[0]), float(elem[1]), float(elem[2]), float(elem[3]), float(elem[4])] for elem in calculated_data]
		self.extra = calculated_data
		selected_data = calculated_data[0]
		if use_clusters_data:
			for item in calculated_data:
				if item[3] < selected_data[3]: # FDR for correlation
					selected_data = item
		self.num = selected_data[0]
		self.corr = selected_data[1]
		self.p_corr = selected_data[2]
		self.p_corr_fdr = selected_data[3]
		self.cause = selected_data[4]

def compute_content_for_chr(ref_file_name, tl_file_name, tl_records_scope):
	window = 100
	cur_chr = [seq_record for seq_record in SeqIO.parse(ref_file_name, "fasta")]
	assert len(cur_chr) == 1, 'More than 1 seq in' + ref_file_name
	cur_chr = cur_chr[0]
	print(cur_chr.id)
	tl_id = len(tl_records_scope)
	tl_line_num = 1
	for line in open(tl_file_name):
		pos = (int(line.split()[1]))
		gc_content = (GC(cur_chr.seq[pos-window:pos+window+1]))
		cpg_content = (cur_chr.seq[pos-window:pos+window+1].upper().count('CG'))
		tl_records_scope[tl_id] = tl_record(gc_content, cpg_content, line, tl_line_num, cur_chr.id, tl_id, True)
		tl_id += 1
		tl_line_num += 1

def compute_content(ref_dir, tl_dir):
	tl_records_scope = {}
	tl_dir_list = os.listdir(tl_dir)
	for tl_name in tl_dir_list:
		print 'processing', os.path.basename(tl_name)
		cur_chr_number = os.path.basename(tl_name).split('_')[0]
		if cur_chr_number == '23':
			chr_name = 'chrX.fa'
		elif cur_chr_number == '24':
			chr_name = 'chrY.fa'
		else:
			chr_name = 'chr' + cur_chr_number + '.fa'
		compute_content_for_chr(os.path.join(ref_dir, chr_name), os.path.join(tl_dir, tl_name), tl_records_scope)
	return tl_records_scope

def split_for_background(tl_records_scope):
	p_fdr_threshold = 0.2
	tl_rec = {}
	bg_rec = {}
	for (key, tl) in tl_records_scope.iteritems():
		if tl.p_corr_fdr < p_fdr_threshold:
			tl_rec[key] = tl
		else:
			bg_rec[key] = tl
	return (bg_rec, tl_rec)

def background_for_gc_and_cpg(tl_background):
	background_dict_gc = {}
	background_dict_cpg = {}
	for (tl_id, tl) in tl_background.iteritems():
		if not background_dict_gc.has_key(tl.gc_content):
			background_dict_gc[tl.gc_content] = []
		background_dict_gc[tl.gc_content].append(tl_id)
		if not background_dict_cpg.has_key(tl.cpg_content):
			background_dict_cpg[tl.cpg_content] = []
		background_dict_cpg[tl.cpg_content].append(tl_id)
	return (background_dict_gc, background_dict_cpg)

def find_pair(tl, tl_records_scope):
	(background_dict_gc, background_dict_cpg) = background_for_gc_and_cpg(tl_records_scope)
	gc_max_difference = 0.05
	sorted_gc = sorted(background_dict_gc.keys())
	gc_index = bisect.bisect_left(sorted_gc, tl.gc_content)
	gc_set = Set()
	pos = gc_index
	while pos >= 0 and sorted_gc[pos] >= tl.gc_content * (1 - gc_max_difference):
		for tl_id in background_dict_gc[sorted_gc[pos]]:
			gc_set.add(tl_id)
		pos -= 1
	pos = gc_index
	while pos < len(sorted_gc) and sorted_gc[pos] <= tl.gc_content * (1 + gc_max_difference):
		for tl_id in background_dict_gc[sorted_gc[pos]]:
			gc_set.add(tl_id)
		pos += 1

	cpg_max_difference = 0.05
	sorted_cpg = sorted(background_dict_cpg.keys())
	cpg_index = bisect.bisect_left(sorted_cpg, tl.cpg_content)
	cpg_set = Set()
	pos = cpg_index
	while pos >= 0 and sorted_cpg[pos] >= tl.cpg_content * (1 - cpg_max_difference):
		for tl_id in background_dict_cpg[sorted_cpg[pos]]:
			cpg_set.add(tl_id)
		pos -= 1
	pos = cpg_index
	while pos < len(sorted_cpg) and sorted_cpg[pos] <= tl.cpg_content * (1 + cpg_max_difference):
		for tl_id in background_dict_cpg[sorted_cpg[pos]]:
			cpg_set.add(tl_id)
		pos += 1
	result = gc_set.intersection(cpg_set)
	return result

def build_tl_pairs(main_tl, background_tl, outf):
	outf = open(outf, 'w')
	pairs_list = []
	bg_used_id = Set()
	cur_num = 0
	for (tl_id, tl) in main_tl.iteritems():
		cur_num += 1
		background_records_id = find_pair(tl, background_tl)
		bg_record = None
		while len(background_records_id) > 0:
			possible_id_index = random.randint(0, len(background_records_id)-1)
			possible_id = list(background_records_id)[possible_id_index]
			if not possible_id in bg_used_id:
				bg_record = background_tl[possible_id]
				outf.write(tl.chr + ':' + str(tl.line_num) + '\t' + bg_record.chr + ':' + str(bg_record.line_num) + '\n')
				print 'n=', cur_num, tl.tl_id, bg_record.tl_id
				bg_used_id.add(possible_id)
#				del background_tl[possible_id]
				break
			else:
				background_records_id.remove(possible_id)
		if not bg_record is None:
			pairs_list.append((tl, bg_record))
	outf.close()
	return pairs_list

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print 'Usage:', sys.argv[0], '-r <reference directory> -t <traffic lights directory> -o <output directory>'
		exit()

	random.seed(6)
	
	parser = argparse.ArgumentParser(prog = sys.argv[0], description='GC and CpG content count')
	parser.add_argument('-r', '--ref_dir', help='reference directory', required=True)
	parser.add_argument('-t', '--tl_dir', help='traffic lights directory', required=True)
	parser.add_argument('-o', '--out_dir', help='output directory', required=True)
	args = parser.parse_args()
	ref_dir = args.ref_dir
	tl_dir = args.tl_dir
	out_dir = args.out_dir

	if not os.path.isdir(ref_dir):
		print >> sys.stderr, 'Not a directory ' + ref_dir
		sys.exit(1)

	if not os.path.isdir(tl_dir):
		print >> sys.stderr, 'Not a directory ' + tl_dir
		sys.exit(1)

	if not os.path.exists(out_dir):
		os.mkdir(out_dir)

	tl_records_scope = compute_content(ref_dir, tl_dir)
	(bg, tl) = split_for_background(tl_records_scope)
	outf = os.path.join(out_dir, 'Pairs.txt')
	pairs = build_tl_pairs(tl, bg, outf)
	print 'pairs', len(pairs)


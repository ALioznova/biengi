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
	def __init__(self, gc_content, cpg_content, tl_line, tl_id):
		self.gc_content = gc_content
		self.cpg_content = cpg_content
		self.gene_info = tl_line.split()[0]
		self.pos = int(tl_line.split()[1])
		self.strand = tl_line.split()[2]
		self.annotation = tl_line.split()[3].split(',')
		assert tl_line.split()[4].startswith('(')
		self.num = int(tl_line.split()[4][1:])
		self.corr = float(tl_line.split()[5])
		self.p_corr = float(tl_line.split()[6])
		self.cause = float(tl_line.split()[7])
		assert tl_line.split()[8].endswith(')')
		self.p_cause = float(tl_line.split()[8][:-1])
		self.extra = tl_line.split()[9:]
		self.tl_id = tl_id

def compute_content_for_chr(ref_file_name, tl_file_name, tl_records_scope):
	window = 100
	cur_chr = [seq_record for seq_record in SeqIO.parse(ref_file_name, "fasta")]
	assert len(cur_chr) == 1, 'More than 1 seq in' + ref_file_name
	cur_chr = cur_chr[0]
	print(cur_chr.id)
	tl_id = len(tl_records_scope)
	for line in open(tl_file_name):
		pos = (int(line.split()[1]))
		gc_content = (GC(cur_chr.seq[pos-window:pos+window+1]))
		cpg_content = (cur_chr.seq[pos-window:pos+window+1].upper().count('CG'))
		tl_records_scope[tl_id] = tl_record(gc_content, cpg_content, line, tl_id)
		tl_id += 1

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

def get_annotations(tl_pos_corr, tl_neg_corr):
	annotations = Set()
	for tl_records_scope in (tl_pos_corr, tl_neg_corr):
		for tlr in tl_records_scope.itervalues():
			for a in tlr.annotation:
				annotations.add(a)
	return annotations

def split_for_background(tl_records_scope):
	p_val_threshold = 0.05
	tl_pos_corr = {}
	tl_neg_corr = {}
	tl_background = {}
	for (key, tl) in tl_records_scope.iteritems():
		if tl.p_corr < p_val_threshold:
			if tl.corr > 0:
				tl_pos_corr[key] = tl
			elif tl.corr < 0:
				tl_neg_corr[key] = tl
		else:
			tl_background[key] = tl
	return (tl_background, tl_pos_corr, tl_neg_corr)

def background_for_gc_and_cpg(tl_background):
	background_dict = {}
	for (tl_id, tl) in tl_background.iteritems():
		if not background_dict.has_key(tl.gc_content):
			background_dict[tl.gc_content] = {}
		if not background_dict[tl.gc_content].has_key(tl.cpg_content):
			background_dict[tl.gc_content][tl.cpg_content] = []
		background_dict[tl.gc_content][tl.cpg_content].append(tl_id)
	return background_dict

def get_annotated_records(annotation, tl_records_scope):
	annotated_tl = {}
	for (tl_id, tl) in tl_records_scope.iteritems():
		if annotation in tl.annotation:
			annotated_tl[tl_id] = tl
	return annotated_tl

def find_closest_cg_and_cpg(tl, tl_records_scope):
	background_dict = background_for_gc_and_cpg(tl_records_scope)
	gc_max_difference = 0.05
	sorted_gc = sorted(background_dict.keys())
	gc_index = bisect.bisect_left(sorted_gc, tl.gc_content)
	if sorted_gc[gc_index] != tl.gc_content:
		if gc_index!=0 and gc_index!=len(sorted_gc) and abs(sorted_gc[gc_index-1] - tl.gc_content) < abs(sorted_gc[gc_index] - tl.gc_content) and abs(sorted_gc[gc_index-1] - tl.gc_content) < abs(sorted_gc[gc_index+1] - tl.gc_content):
			gc_index = gc_index-1
		elif gc_index!=0 and gc_index!=len(sorted_gc) and abs(sorted_gc[gc_index+1] - tl.gc_content) < abs(sorted_gc[gc_index] - tl.gc_content) and abs(sorted_gc[gc_index+1] - tl.gc_content) < abs(sorted_gc[gc_index-1] - tl.gc_content):
			gc_index = gc_index+1
	if abs(sorted_gc[gc_index] - tl.gc_content) >= tl.gc_content * gc_max_difference:
		return None
	cpg_max_difference = 0.05
	sorted_cpg = sorted(background_dict[sorted_gc[gc_index]].keys())
	cpg_index = bisect.bisect_left(sorted_cpg, tl.cpg_content)
	if sorted_cpg[cpg_index] != tl.cpg_content:
		if cpg_index!=0 and cpg_index!=len(sorted_cpg) and abs(sorted_cpg[cpg_index-1] - tl.cpg_content) < abs(sorted_cpg[cpg_index] - tl.cpg_content) and abs(sorted_cpg[cpg_index-1] - tl.cpg_content) < abs(sorted_cpg[cpg_index+1] - tl.cpg_content):
			cpg_index = cpg_index-1
		elif cpg_index!=0 and cpg_index!=len(sorted_cpg) and abs(sorted_cpg[cpg_index+1] - tl.cpg_content) < abs(sorted_cpg[cpg_index] - tl.cpg_content) and abs(sorted_cpg[cpg_index+1] - tl.cpg_content) < abs(sorted_cpg[cpg_index-1] - tl.cpg_content):
			cpg_index = cpg_index+1
	if abs(sorted_cpg[cpg_index] - tl.cpg_content) >= tl.cpg_content * cpg_max_difference:
		return None
	closest_tl_id = background_dict[sorted_gc[gc_index]][sorted_cpg[cpg_index]][random.randint(0, len(background_dict[sorted_gc[gc_index]][sorted_cpg[cpg_index]])-1)]
	return tl_records_scope[closest_tl_id]

def build_tl_pairs(main_tl, background_tl):
	main_corr = []
	background_corr = []
	for (tl_id, tl) in main_tl.iteritems():
		background_record = find_closest_cg_and_cpg(tl, background_tl)
		if background_record :
			main_corr.append(tl.corr)
			background_corr.append(background_record.corr)
	return (main_corr, background_corr)

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print 'Usage:', sys.argv[0], '-r <reference directory> -t <traffic lights directory> -o <output directory>'
		exit()

	random.seed(1)
	
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
	(tl_background, tl_pos_corr, tl_neg_corr) = split_for_background(tl_records_scope)
	annotations = get_annotations(tl_pos_corr, tl_neg_corr)
	for annotation in annotations:
		print annotation
		tl_pos_corr_an = get_annotated_records(annotation, tl_pos_corr)
		tl_neg_corr_an = get_annotated_records(annotation, tl_neg_corr)
		tl_background_an = get_annotated_records(annotation, tl_background)
		(pos_corr, background_pos_corr) = build_tl_pairs(tl_pos_corr_an, tl_pos_corr_an)
		out_pos = open(os.path.join(out_dir, annotation + '_pos.txt'), 'w')
		for i in xrange(len(pos_corr)):
			out_pos.write(str(pos_corr[i]) + '\t' + str(background_pos_corr[i]) + '\n')
		out_pos.close()
		(neg_corr, background_neg_corr) = build_tl_pairs(tl_neg_corr_an, tl_neg_corr_an)
		out_neg = open(os.path.join(out_dir, annotation + '_neg.txt'), 'w')
		for i in xrange(len(neg_corr)):
			out_neg.write(str(neg_corr[i]) + '\t' + str(background_neg_corr[i]) + '\n')
		out_neg.close()


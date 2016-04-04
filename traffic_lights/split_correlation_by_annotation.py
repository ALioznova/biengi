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
	def __init__(self, gc_content, cpg_content, tl_line, tl_line_num, tl_chr, tl_id):
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
		assert tl_line.split()[4].startswith('(')
		self.num = int(tl_line.split()[4][1:])
		self.corr = float(tl_line.split()[5])
		self.p_corr = float(tl_line.split()[6])
		self.p_corr_fdr = float(tl_line.split()[7])
		assert tl_line.split()[8].endswith(')')
		self.cause = float(tl_line.split()[8][:-1])
		self.extra = tl_line.split()[9:]
		self.tl_id = tl_id

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
		tl_records_scope[tl_id] = tl_record(gc_content, cpg_content, line, tl_line_num, cur_chr.id, tl_id)
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

def get_annotations(tl_pos_corr, tl_neg_corr):
	annotations = Set()
	for tl_records_scope in (tl_pos_corr, tl_neg_corr):
		for tlr in tl_records_scope.itervalues():
			for a in tlr.annotation.keys():
				annotations.add(a)
	return annotations

def split_for_background(tl_records_scope):
	p_fdr_threshold = 0.2
	tl_pos_corr = {}
	tl_neg_corr = {}
	tl_background = {}
	for (key, tl) in tl_records_scope.iteritems():
		if tl.p_corr_fdr < p_fdr_threshold:
			if tl.corr > 0:
				tl_pos_corr[key] = tl
			elif tl.corr < 0:
				tl_neg_corr[key] = tl
		else:
			tl_background[key] = tl
	return (tl_background, tl_pos_corr, tl_neg_corr)

def background_for_gc_and_cpg_ann(tl_background, annotation):
	background_dict_gc = {}
	background_dict_cpg = {}
	background_dict_ann = {}
	for (tl_id, tl) in tl_background.iteritems():
		if not background_dict_gc.has_key(tl.gc_content):
			background_dict_gc[tl.gc_content] = []
		background_dict_gc[tl.gc_content].append(tl_id)
		if not background_dict_cpg.has_key(tl.cpg_content):
			background_dict_cpg[tl.cpg_content] = []
		background_dict_cpg[tl.cpg_content].append(tl_id)
		if not background_dict_ann.has_key(tl.annotation[annotation]):
			background_dict_ann[tl.annotation[annotation]] = []
		background_dict_ann[tl.annotation[annotation]].append(tl_id)
	return (background_dict_gc, background_dict_cpg, background_dict_ann)

def get_annotated_records(annotation, tl_records_scope):
	annotated_tl = {}
	for (tl_id, tl) in tl_records_scope.iteritems():
		if annotation in tl.annotation.keys():
			annotated_tl[tl_id] = tl
	return annotated_tl

def find_closest_cg_and_cpg_an(tl, tl_records_scope, annotation):
	(background_dict_gc, background_dict_cpg, background_dict_ann) = background_for_gc_and_cpg_ann(tl_records_scope, annotation)
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

	ann_max_difference = 0.05
	sorted_ann = sorted(background_dict_ann.keys())
	ann_set = None
	if len(sorted_ann) > 1:
		ann_index = bisect.bisect_left(sorted_ann, tl.annotation[annotation])
		ann_set = Set()
		pos = ann_index
		while pos >= 0 and sorted_ann[pos] >= tl.annotation[annotation] * (1 - ann_max_difference):
			for tl_id in background_dict_ann[sorted_ann[pos]]:
				ann_set.add(tl_id)
			pos -= 1
		pos = ann_index
		while pos < len(sorted_ann) and sorted_ann[pos] <= tl.annotation[annotation] * (1 + ann_max_difference):
			for tl_id in background_dict_ann[sorted_ann[pos]]:
				ann_set.add(tl_id)
			pos += 1

	result = gc_set.intersection(cpg_set)
	if ann_set:
		result = result.intersection(ann_set)
	return result

def build_tl_pairs(main_tl, background_tl, annotation):
	pairs_list = []
	bg_used_id = Set()
	for (tl_id, tl) in main_tl.iteritems():
		background_records_id = find_closest_cg_and_cpg_an(tl, background_tl, annotation)
		bg_record = None
		while len(background_records_id) > 0:
			possible_id_index = random.randint(0, len(background_records_id)-1)
			possible_id = list(background_records_id)[possible_id_index]
			if not possible_id in bg_used_id:
				bg_record = background_tl[possible_id]
				bg_used_id.add(possible_id)
				break
			else:
				background_records_id.remove(possible_id)
		if not bg_record is None:
			pairs_list.append((tl, bg_record))
	return pairs_list

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
	print 'bg', len(tl_background), '+', len(tl_pos_corr), '-', len(tl_neg_corr)
	annotations = get_annotations(tl_pos_corr, tl_neg_corr)
	for annotation in annotations:
		print annotation
#		if annotation in ['phyloP20way', 'phastCons20way', 'repeats', 'Intron', 'wgEncodeRegDnaseClustered']:
#			continue
		tl_pos_corr_an = get_annotated_records(annotation, tl_pos_corr)
		tl_neg_corr_an = get_annotated_records(annotation, tl_neg_corr)
		tl_background_an = get_annotated_records(annotation, tl_background)
		print 'TL', len(tl_pos_corr_an) + len (tl_neg_corr_an), 'of', len(tl_pos_corr) + len(tl_neg_corr)
		print 'bg', len(tl_background_an), 'of', len(tl_background)
		
		pos_pairs = build_tl_pairs(tl_pos_corr_an, tl_background_an, annotation)
		print 'pos pairs selected'
		out_pos = open(os.path.join(out_dir, annotation + '_pos.txt'), 'w')
		for i in xrange(len(pos_pairs)):
			(tl_record, bg_record) = pos_pairs[i]
			tl_annotation_val = str(tl_record.annotation[annotation])
			tl_corr = str(tl_record.corr)
			bg_corr = str(bg_record.corr)
			tl_chr = tl_record.chr
			tl_line_num = str(tl_record.line_num)
			bg_chr = bg_record.chr
			bg_line_num = str(bg_record.line_num)
			out_pos.write(tl_corr + '\t' + bg_corr + '\t' + tl_annotation_val + '\t' + tl_chr + ':' + tl_line_num + '\t' + bg_chr + ':' + bg_line_num + '\n')
		out_pos.close()
		
		neg_pairs = build_tl_pairs(tl_neg_corr_an, tl_background_an, annotation)
		print 'neg pairs selected'
		out_neg = open(os.path.join(out_dir, annotation + '_neg.txt'), 'w')
		for i in xrange(len(neg_pairs)):
			(tl_record, bg_record) = neg_pairs[i]
			tl_annotation_val = str(tl_record.annotation[annotation])
			tl_corr = str(tl_record.corr)
			bg_corr = str(bg_record.corr)
			tl_chr = tl_record.chr
			tl_line_num = str(tl_record.line_num)
			bg_chr = bg_record.chr
			bg_line_num = str(bg_record.line_num)
			out_neg.write(tl_corr + '\t' + bg_corr + '\t' + tl_annotation_val + '\t' + tl_chr + ':' + tl_line_num + '\t' + bg_chr + ':' + bg_line_num + '\n')
		out_neg.close()
		print

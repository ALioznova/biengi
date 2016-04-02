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
	if len(result)>0:
		return result
	else:
		return None

def build_tl_pairs(main_tl, background_tl, annotation):
	main_corr = []
	background_corr = []
	annotation_val = []
	bg_used_id = Set()
	for (tl_id, tl) in main_tl.iteritems():
		background_records_id = find_closest_cg_and_cpg_an(tl, background_tl, annotation)
		if background_records_id :
			bg_record = None
			while not bg_record and len(background_records_id)>0:
				possible_id_index = random.randint(0, len(background_records_id)-1)
				possible_id = list(background_records_id)[possible_id_index]
				if not possible_id in bg_used_id:
					bg_record = background_tl[possible_id]
					bg_used_id.add(possible_id)
					break
				else:
					background_records_id.remove(possible_id)
			if bg_record:
				main_corr.append(tl.corr)
				background_corr.append(bg_record.corr)
				annotation_val.append(tl.annotation[annotation])
	return (main_corr, background_corr, annotation_val)

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
#		if annotation in ['Intron', 'Encode_TRIM28', 'Encode_MYBL2', 'Encode_SMARCC1', 'Encode_SMARCC2', 'Encode_FAM48A', 'Encode_ESR1', 'Encode_MBD4', 'Encode_GATA1', 'Encode_GATA2', 'Encode_GATA3', 'Encode_BACH1', 'Encode_NR2F2', 'Encode_CHD2', 'Encode_CHD1', 'Encode_TFAP2C', 'Encode_FOSL2', 'Encode_TFAP2A', 'Encode_BATF', 'Encode_eGFP-GATA2', 'Encode_ELK4', 'Encode_FOXP2', 'Encode_GTF3C2', 'Encode_ELK1', 'Encode_SMARCB1', 'Encode_SMC3', 'Encode_IKZF1', 'Encode_FOS', 'Encode_ESRRA', 'Encode_ZNF217', 'Encode_TCF7L2', 'Encode_HDAC2', 'Encode_HDAC1', 'Encode_MXI1', 'Encode_eGFP-HDAC8', 'Encode_HDAC6', 'Encode_ZNF274', 'Encode_JUN', 'Encode_SUZ12', 'Encode_BDP1', 'Encode_KDM5A', 'Encode_GRp20', 'Encode_PHF8', 'Encode_PRDM1', 'Encode_NFE2', 'Encode_ZBTB33', 'Encode_YY1', 'Encode_GABPA', 'Encode_RAD21', 'Exon', 'Encode_BHLHE40', 'Encode_FOXA2', 'Encode_FOXA1', 'Encode_CTCFL', 'Encode_HSF1', 'Encode_HNF4A', 'Encode_BCL11A', 'Encode_EBF1', 'Encode_ZNF263', 'Encode_CEBPD', 'Encode_BCLAF1', 'Encode_CEBPB', 'Encode_PAX5', 'Encode_RCOR1', 'Encode_MTA3', 'Enhancer', 'Encode_RELA', 'Encode_NRF1', 'Encode_SREBP1', 'Encode_MAFF', 'Encode_MAFK', 'CpGIsland', 'Encode_NR2C2', 'CpGBeacon', 'Encode_CBX3', 'Encode_USF1', 'Encode_USF2', 'Encode_THAP1', 'Encode_CTBP2', 'Encode_TAL1', 'Encode_CCNT2', 'Encode_SIX5', 'Encode_MYC', 'Encode_HMGN3', 'Encode_TCF12', 'Encode_POLR2A', 'Encode_KAP1', 'Encode_ZNF143', 'Encode_ZZZ3', 'Encode_ELF1', 'Encode_BRF2', 'Encode_ETS1', 'Encode_EGR1', 'Encode_SP4', 'Encode_SP2', 'Encode_EZH2', 'Encode_SP1', 'Encode_POU5F1', 'Encode_SPI1', 'Encode_eGFP-FOS', 'Encode_MAZ', 'Encode_MAX', 'Encode_NANOG', 'Encode_EP300', 'Encode_TBL1XR1', 'Encode_BCL3', 'Encode_HA-E2F1', 'Encode_FOXM1', 'Encode_TEAD4', 'Encode_RDBP', 'Encode_SIRT6', 'Outside', 'Encode_REST', 'Encode_ZKSCAN1', 'Encode_NR3C1', 'Encode_FOSL1', 'Encode_RXRA', 'Encode_E2F4', 'Encode_E2F6', 'Encode_E2F1', 'Encode_ZEB1', 'Encode_RUNX3', 'Encode_WRNIP1', 'Encode_NFYA', 'Intron', 'Encode_NFYB', 'Encode_ATF3', 'Encode_ATF2', 'Encode_ATF1', 'Encode_SETDB1', 'Encode_CREB1', 'Encode_GTF2F1', 'Encode_TAF1', 'Encode_SIN3AK20', 'Encode_STAT5A', 'Encode_UBTF', 'Encode_ARID3A', 'Encode_SIN3A', 'Encode_SMARCA4', 'Encode_IRF1', 'Encode_IRF3', 'Encode_IRF4', 'Encode_GTF2B', 'Encode_CTCF', 'Encode_TAF7', 'Encode_PML', 'Encode_PPARGC1A', 'Encode_RFX5', 'Encode_eGFP-JUND', 'Encode_BRCA1', 'Encode_SAP30', 'Encode_eGFP-JUNB', 'Encode_RBBP5', 'Promoter', 'Encode_JUND', 'Encode_MEF2A', 'Encode_MEF2C', 'Encode_NFATC1', 'Encode_HNF4G', 'Encode_NFIC', 'Encode_POU2F2', 'Encode_SRF', 'Encode_TBP', 'Encode_KDM5B', 'Encode_ZBTB7A', 'Encode_TCF3', 'Encode_PBX3', 'Encode_STAT1', 'Encode_STAT3', 'Encode_STAT2', 'Encode_RPC155']:
#			continue
		tl_pos_corr_an = get_annotated_records(annotation, tl_pos_corr)
		tl_neg_corr_an = get_annotated_records(annotation, tl_neg_corr)
		tl_background_an = get_annotated_records(annotation, tl_background)
		print 'TL', len(tl_pos_corr_an) + len (tl_neg_corr_an), 'of', len(tl_pos_corr) + len(tl_neg_corr)
		print 'bg', len(tl_background_an), 'of', len(tl_background)
		
		(pos_corr, background_pos_corr, annotation_val_pos) = build_tl_pairs(tl_pos_corr_an, tl_background_an, annotation)
		print 'pos pairs selected'
		out_pos = open(os.path.join(out_dir, annotation + '_pos.txt'), 'w')
		for i in xrange(len(pos_corr)):
			out_pos.write(str(pos_corr[i]) + '\t' + str(background_pos_corr[i]) + 't' + str(annotation_val_pos[i]) + '\n')
		out_pos.close()
		
		(neg_corr, background_neg_corr, annotation_val_neg) = build_tl_pairs(tl_neg_corr_an, tl_background_an, annotation)
		print 'neg pairs selected'
		out_neg = open(os.path.join(out_dir, annotation + '_neg.txt'), 'w')
		for i in xrange(len(neg_corr)):
			out_neg.write(str(neg_corr[i]) + '\t' + str(background_neg_corr[i]) + '\t' + str(annotation_val_neg[i]) + '\n')
		out_neg.close()
		print

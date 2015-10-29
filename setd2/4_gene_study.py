#!/usr/bin/env python
# -*- coding: utf-8 -*- 
import os
import sys
import argparse
import numpy
from sets import Set
from intervaltree import Interval, IntervalTree
from enum import Enum, unique

class Pos_rec:
	def __init__(self, description):
		(chr_name, begin1, end1, begin2, end2, strand) = description.split(':')
		self.desc = description
		self.chr_name = chr_name
		self.beg1 = int(begin1)
		self.end1 = int(end1)
		self.beg2 = int(begin2)
		self.end2 = int(end2)
		if strand == '+':
			self.strand = True
		else:
			self.strand = False

def read_data(in_fn):
	in_f = open(in_fn)
	header_line = in_f.readline()
	tumor_setd2_broken_num = int((header_line.split()[1]).split(':')[1])
	tumor_num = int((header_line.split()[2]).split(':')[1])
	normal_num = int((header_line.split()[3]).split(':')[1])
	pos = []
	tumor_setd2_broken = []
	tumor = []
	normal = []
	for line in in_f:
		pos.append(Pos_rec(line.split()[0]))
		tumor_setd2_broken.append(float(line.split()[1]))
		tumor.append(float(line.split()[2]))
		normal.append(float(line.split()[3]))
	in_f.close()
	return (pos, tumor_setd2_broken_num, tumor_setd2_broken, tumor_num, tumor, normal_num, normal)

class Genes:
	def __init__(self, name, locus, exons):
		self.name = name
		self.locus = locus
		self.exons = exons

class Locus:
	def __init__(self, description):
		(chr_name, coords, strand) = description.split(':')
		self.chr_name = chr_name
		self.beg = int(coords.split('-')[0])
		self.end = int(coords.split('-')[1])
		if strand == '+':
			self.strand = True
		else:
			self.strand = False
		
class Exons:
	def __init__(self, description):
		(chr_name, coords, strand) = description.split(':')
		self.chr_name = chr_name
		if strand == '+':
			self.strand = True
		else:
			self.strand = False
		self.begs = [int(c.split('-')[0]) for c in coords.split(',')]
		self.ends = [int(c.split('-')[1]) for c in coords.split(',')]

def parse_genes_file(genes_file):
	# data source: https://tcga-data.nci.nih.gov/docs/GAF/GAF.hg19.June2011.bundle/outputs/TCGA.hg19.June2011.gaf
	gf = open(genes_file)
	gf.readline()
	genes_data = {}
	exon_dict = {}
	interval_trees = {}
	for line in gf:
		if line.split('\t')[2] == 'gene':
			name = line.split('\t')[1]
			exons = Exons(line.split('\t')[14])
			locus = [Locus(l) for l in line.split('\t')[16].split(';')]
			# filling genes_data: gene_name -> [Genes]
			if len(name.split('|')) > 2: # gene has multiple locuses
				if not genes_data.has_key('|'.join(name.split('|')[:-1])):
					genes_data['|'.join(name.split('|')[:-1])] = []
				genes_data['|'.join(name.split('|')[:-1])].append(Genes(name, locus, exons))
			else:
				genes_data[name] = [Genes(name, locus, exons)] # the only line with gene description
			if not exon_dict.has_key((exons.chr_name, exons.strand)):
				exon_dict[(exons.chr_name, exons.strand)] = {}
			if not interval_trees.has_key((exons.chr_name, exons.strand)):
				interval_trees[(exons.chr_name, exons.strand)] = IntervalTree()
			for i in xrange(len(exons.begs)):
				# filling exon_dict: (chr, strand) -> (exon_start, exon_end) -> Set(gene_name)
				if exons.begs[i] == exons.ends[i]:
					continue
				if not exon_dict[(exons.chr_name, exons.strand)].has_key((exons.begs[i], exons.ends[i])):
					exon_dict[(exons.chr_name, exons.strand)][(exons.begs[i], exons.ends[i])] = Set()
				exon_dict[(exons.chr_name, exons.strand)][(exons.begs[i], exons.ends[i])].add('|'.join(name.split('|')[:2]))
				# filling interval_trees for annotated exons:
				interval_trees[(exons.chr_name, exons.strand)][exons.begs[i] : exons.ends[i]] = '|'.join(name.split('|')[:2])
	gf.close()

	return (genes_data, exon_dict, interval_trees)

def exon_to_gene_names(exon_dict, interval_trees, pos):
	exon_names = {}
	for p in pos:
		if not exon_names.has_key((p.chr_name, p.strand)):
			exon_names[(p.chr_name, p.strand)] = {}
		if exon_dict[(p.chr_name, p.strand)].has_key((p.end1, p.beg2)):
			exon_names[(p.chr_name, p.strand)][(p.end1, p.beg2)] = exon_dict[(p.chr_name, p.strand)][(p.end1, p.beg2)]
		else:
			left_genes = Set([lg[2] for lg in interval_trees[(p.chr_name, p.strand)][p.end1 + 0.1]])
			right_genes = Set([rg[2] for rg in interval_trees[(p.chr_name, p.strand)][p.beg2 - 0.1]])
			if len(left_genes) != 0 or len(right_genes) != 0:
				if len(left_genes.intersection(right_genes)) != 0:
					exon_names[(p.chr_name, p.strand)][(p.end1, p.beg2)] = left_genes.intersection(right_genes)
	return exon_names

def get_dist_exon_gene_beg(exon_to_genes, genes_data, pos_rec):
	def get_overlap(a, b):
		return max(0, min(a[1], b[1]) - max(a[0], b[0]))

	exon_chr = pos_rec.chr_name
	exon_strand = pos_rec.strand
	exon_beg = pos_rec.end1
	exon_end = pos_rec.beg2
	if not exon_to_genes[(exon_chr, exon_strand)].has_key((exon_beg, exon_end)):
		return ([], [], [])
	exon_genes = list(exon_to_genes[(exon_chr, exon_strand)][(exon_beg, exon_end)])
	dist_bp = []
	dist_perc = []
	gene_name = []
	for exon_gene in exon_genes:
		for gene in genes_data[exon_gene]:
			if gene.exons.chr_name != exon_chr:
				continue
			for i in xrange(len(gene.exons.begs)):
				# if there is no intersection with exon_beg, exon_end => it is other exon's locus
				ge_beg = gene.exons.begs[i]
				ge_end = gene.exons.ends[i]
				if get_overlap((ge_beg, ge_end), (exon_beg, exon_end)) != 0:
					for loc in gene.locus:
						if loc.chr_name != gene.exons.chr_name:
							continue
						if get_overlap((loc.beg, loc.end), (ge_beg, ge_end)) != 0:
							if loc.strand == True and exon_strand == True:
								if exon_beg - loc.beg < 0 or exon_beg - loc.beg > loc.end - loc.beg:
									continue
								dist_bp.append(exon_beg - loc.beg)
								dist_perc.append(float(exon_beg - loc.beg)/(loc.end - loc.beg))
								gene_name.append(gene.name)
							elif loc.strand == False and exon_strand == False:
								if loc.end - exon_end < 0 or loc.end - exon_end > loc.end - loc.beg:
									continue
								dist_bp.append(loc.end - exon_end)
								dist_perc.append(float(loc.end - exon_end)/(loc.end - loc.beg))
								gene_name.append(gene.name)
							else:
								continue
	return (dist_bp, dist_perc, gene_name)

def get_gene_dist(exon_to_genes, genes_data, pos):
	gene_beg_dist_bp = []
	gene_beg_dist_perc = []
	gene_names = []
	for i in xrange(len(pos)):
		(dist_bp, dist_perc, gene_name) = get_dist_exon_gene_beg(exon_to_genes, genes_data, pos[i])
		dist_bp = list(Set(dist_bp))
		dist_perc = list(Set(dist_perc))
		gene_name = list(Set(gene_name))
		gene_beg_dist_bp.append(dist_bp)
		gene_beg_dist_perc.append(dist_perc)
		gene_names.append(gene_name)
	return (gene_beg_dist_bp, gene_beg_dist_perc, gene_names)

class OrderedEnum(Enum):
	 def __ge__(self, other):
		 if self.__class__ is other.__class__:
			 return self.value >= other.value
		 return NotImplemented
	 def __gt__(self, other):
		 if self.__class__ is other.__class__:
			 return self.value > other.value
		 return NotImplemented
	 def __le__(self, other):
		 if self.__class__ is other.__class__:
			 return self.value <= other.value
		 return NotImplemented
	 def __lt__(self, other):
		 if self.__class__ is other.__class__:
			 return self.value < other.value
		 return NotImplemented

@unique
class Impact(OrderedEnum):
	high = 3
	moderate = 2
	low = 1
	modifier = 0
	no = -1
	unknown = -2

def get_mutation_impact_dict():
	# information source: http://www.ensembl.org/info/genome/variation/predicted_data.html#consequences
	# useful tool: http://mutationassessor.org/
	# one more tool: http://stothard.afns.ualberta.ca/downloads/NGS-SNP/annotate_SNPs.html
	# and one more tool: ftp://hgdownload.cse.ucsc.edu/apache/htdocs-rr/goldenpath/help/hgVaiHelpText.html
	mutation_impact_dict = {}
	mutation_impact_dict['transcript_ablation'] = Impact.high
	mutation_impact_dict['splice_acceptor_variant'] = Impact.high
	mutation_impact_dict['splice_donor_variant'] = Impact.high
	mutation_impact_dict['stop_gained'] = Impact.high
	mutation_impact_dict['frameshift_variant'] = Impact.high
	mutation_impact_dict['stop_lost'] = Impact.high
	mutation_impact_dict['start_lost'] = Impact.high
	mutation_impact_dict['transcript_amplification'] = Impact.high
	mutation_impact_dict['inframe_insertion'] = Impact.moderate
	mutation_impact_dict['inframe_deletion'] = Impact.moderate
	mutation_impact_dict['missense_variant'] = Impact.moderate
	mutation_impact_dict['missense'] = Impact.moderate # added
	mutation_impact_dict['protein_altering_variant'] = Impact.moderate
	mutation_impact_dict['splice_region_variant'] = Impact.low
	mutation_impact_dict['incomplete_terminal_codon_variant'] = Impact.low
	mutation_impact_dict['stop_retained_variant'] = Impact.low
	mutation_impact_dict['synonymous_variant'] = Impact.low
	mutation_impact_dict['coding_sequence_variant'] = Impact.modifier
	mutation_impact_dict['mature_miRNA_variant'] = Impact.modifier
	mutation_impact_dict['5_prime_UTR_variant'] = Impact.modifier
	mutation_impact_dict['3_prime_UTR_variant'] = Impact.modifier # added
	mutation_impact_dict['non_coding_transcript_exon_variant'] = Impact.modifier
	mutation_impact_dict['intron_variant'] = Impact.modifier
	mutation_impact_dict['NMD_transcript_variant'] = Impact.modifier
	mutation_impact_dict['non_coding_transcript_variant'] = Impact.modifier
	mutation_impact_dict['upstream_gene_variant'] = Impact.modifier
	mutation_impact_dict['downstream_gene_variant'] = Impact.modifier
	mutation_impact_dict['TFBS_ablation'] = Impact.modifier
	mutation_impact_dict['TFBS_amplification'] = Impact.modifier
	mutation_impact_dict['TF_binding_site_variant'] = Impact.modifier
	mutation_impact_dict['regulatory_region_ablation'] = Impact.moderate
	mutation_impact_dict['regulatory_region_amplification'] = Impact.modifier
	mutation_impact_dict['feature_elongation'] = Impact.modifier
	mutation_impact_dict['regulatory_region_variant'] = Impact.modifier
	mutation_impact_dict['feature_truncation'] = Impact.modifier
	mutation_impact_dict['intergenic_variant'] = Impact.modifier
	return mutation_impact_dict

def find_setd2_mutation_impact(s_file_name, mutation_impact_dict):
	def find_mutation_type_column(line):
		# maf common header: https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification
		# oncotator help: https://www.broadinstitute.org/oncotator/help/
		for i in xrange(len(line)): 
			if line[i] == 'i_Ensembl_so_term':
				return i
		return None

	maf_file = open(s_file_name)
	maf_file.readline()
	maf_file.readline()
	maf_file.readline()
	line = maf_file.readline().split('\t')
	mutation_type_index = find_mutation_type_column(line)
	if mutation_type_index == None:
		return Impact.unknown
	impacts = [Impact.no]
	for line in maf_file:
		if line.startswith('SETD2\t'):
			mutation_type = line.split('\t')[mutation_type_index]
			if not mutation_impact_dict.has_key(mutation_type):
				print 'Unknown mutation type:', mutation_type
				continue
			setd2_mutation_impact = mutation_impact_dict[mutation_type]
			impacts.append(setd2_mutation_impact)
	maf_file.close()
	return max(impacts)

def get_setd2_mutation_rate(maf_dir, sample_names):
	setd2_mutation_rate = {}
	mutation_impact_dict = get_mutation_impact_dict()
	for sample in sample_names:
		s_file_name = os.path.join(maf_dir, sample[:15] + '.hg19.oncotator.hugo_entrez_remapped.maf.txt')
		if not os.path.exists(s_file_name):
			continue
		mutation_impact = find_setd2_mutation_impact(s_file_name, mutation_impact_dict)
		setd2_mutation_rate[sample] = mutation_impact
	return setd2_mutation_rate

class Gene_expr:
	def __init__(self, name, tumor, tumor_num, tumor_setd2, tumor_setd2_num, normal, normal_num, val_arr):
		self.name = name
		self.tumor = tumor
		self.tumor_num = tumor_num
		self.tumor_setd2 = tumor_setd2
		self.tumor_setd2_num = tumor_setd2_num
		self.normal = normal
		self.normal_num = normal_num
		self.val_arr = val_arr

class Sample_type(Enum):
	unknown = -1
	norma = 0
	tumor_wild_type = 1
	tumor_unclassified = 2
	tumor_setd2 = 3

def get_sample_classification(sample_names):
	sample_classification = {}
	setd2_mutation_rates = get_setd2_mutation_rate(maf_dir, sample_names)
	# barcodes https://wiki.nci.nih.gov/display/TCGA/TCGA+Barcode
	# code tables: https://tcga-data.nci.nih.gov/datareports/codeTablesReport.htm?codeTable=sample%20type
	sample_type = [int(s.split('-')[3][:2]) for s in sample_names]
	for i in xrange(len(sample_type)):
		if (sample_type[i] >= 1) and (sample_type[i] <= 9): # tumor
			if setd2_mutation_rates.has_key(sample_names[i]):
				if setd2_mutation_rates[sample_names[i]] == Impact.high:
					sample_classification[sample_names[i]] = Sample_type.tumor_setd2
				elif setd2_mutation_rates[sample_names[i]] == Impact.no:
					sample_classification[sample_names[i]] = Sample_type.tumor_wild_type
				else:
					sample_classification[sample_names[i]] = Sample_type.tumor_unclassified
			else:
				sample_classification[sample_names[i]] = Sample_type.tumor_unclassified
		elif (sample_type[i] >= 10) and (sample_type[i] <= 19): # norma
			sample_classification[sample_names[i]] = Sample_type.norma
		else:
			sample_classification[sample_names[i]] = Sample_type.unknown
	return sample_classification

def process_gene_expression_file(gene_expression_fn, maf_dir):
	fin = open(gene_expression_fn)
	sample_names = fin.readline().split('\t')[1:]
	sample_classification = get_sample_classification(sample_names)
	fin.readline()
	gene_expression = {}
	for line in fin:
		gene_name = line.split('\t')[0]
		if gene_name.endswith('_calculated'):
			gene_name = gene_name[:-len('_calculated')]
		data = line.split('\t')[1:]
		tumor_setd2_num = 0
		tumor_num = 0
		normal_num = 0
		tumor_setd2_val = 0
		tumor_val = 0
		normal_val = 0
		val_arr = []
		for i in xrange(len(data)):
			if i % 3 == 2:
				if sample_classification[sample_names[i]] == Sample_type.tumor_setd2:
					tumor_setd2_val += (float(data[i]))
					tumor_setd2_num += 1
					val_arr.append((tumor_setd2_val, sample_names[i]))
				elif sample_classification[sample_names[i]] == Sample_type.tumor_wild_type:
					tumor_val += (float(data[i]))
					tumor_num += 1
					val_arr.append((tumor_val, sample_names[i]))
				elif sample_classification[sample_names[i]] == Sample_type.norma:
					normal_val += (float(data[i]))
					normal_num += 1
					val_arr.append((normal_val, sample_names[i]))
				else:
					continue
		if tumor_setd2_num != 0:
			tumor_setd2_v = (float(tumor_setd2_val)/tumor_setd2_num)
		else:
			tumor_setd2_v = (float('nan'))
		if tumor_num != 0:
			tumor_v = (float(tumor_val)/tumor_num)
		else:
			tumor_v = (float('nan'))
		if normal_num != 0:
			normal_v = (float(normal_val)/normal_num)
		else:
			normal_v = (float('nan'))
		gene_expression[gene_name] = Gene_expr(gene_name, tumor_v, tumor_num, tumor_setd2_v, tumor_setd2_num, normal_v, normal_num, val_arr)
	fin.close()
	return (gene_expression, sample_classification)

def get_gene_expression(gene_expression_fn, exon_to_genes, pos, maf_dir):
	(gene_expression_data, sample_classification) = process_gene_expression_file(gene_expression_fn, maf_dir)
	gene_expr = Set()
	tumor_setd2 = []
	tumor = []
	normal = []
	for p in pos:
		tumor_setd2_cur = []
		tumor_cur = []
		normal_cur = []
		if exon_to_genes[(p.chr_name, p.strand)].has_key((p.end1, p.beg2)):
			for gene in list(exon_to_genes[(p.chr_name, p.strand)][(p.end1, p.beg2)]):
				if not gene_expression_data.has_key(gene):
					continue
				tumor_setd2_cur.append(gene_expression_data[gene].tumor_setd2)
				tumor_cur.append(gene_expression_data[gene].tumor)
				normal_cur.append(gene_expression_data[gene].normal)
				gene_expr.add(gene_expression_data[gene])
		tumor_setd2.append(tumor_setd2_cur)
		tumor.append(tumor_cur)
		normal.append(normal_cur)
	return (gene_expr, sample_classification, tumor_setd2, tumor, normal)

def output_expression_and_dist_file(out_fn, pos, gene_beg_dist_bp, gene_beg_dist_perc, gene_expression_fn, tumor_setd2_gene_expr, tumor_gene_expr, normal_gene_expr, gene_names):
	fout = open(out_fn, 'w')
	fout.write('Pos:\tGene_dist_bp:\tGene_dist_perc:\tExpr_tumorsetd2:\tExpr_tumor_wt:\tExpr_normal:\tGene_name:\n')
	for i in xrange(len(pos)):
		new_line = pos[i].desc + '\t'
		for dst in gene_beg_dist_bp[i]:
			new_line += (str(dst) + ',')
		if len(gene_beg_dist_bp[i]) != 0:
			new_line = new_line[:-1]
		new_line += '\t'
		for dst in gene_beg_dist_perc[i]:
			new_line += (str(dst) + ',')
		if len(gene_beg_dist_perc[i]) != 0:
			new_line = new_line[:-1]
		if gene_expression_fn:
			new_line += '\t'
			for expr in tumor_setd2_gene_expr[i]:
				new_line += (str(expr) + ',')
			if len(tumor_setd2_gene_expr[i]) != 0:
				new_line = new_line[:-1]
			new_line += '\t'
			for expr in tumor_gene_expr[i]:
				new_line += (str(expr) + ',')
			if len(tumor_gene_expr[i]) != 0:
				new_line = new_line[:-1]
			new_line += '\t'
			for expr in normal_gene_expr[i]:
				new_line += (str(expr) + ',')
			if len(normal_gene_expr[i]) != 0:
				new_line = new_line[:-1]
		else:
			new_line += '\t\t\t'
		new_line += '\t'
		for name in gene_names[i]:
			new_line += (name + ',')
		if len(gene_names[i]) != 0:
			new_line = new_line[:-1]
		new_line += '\n'
		fout.write(new_line)
	fout.close()

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print 'Usage:', sys.argv[0], '-d <data directory> -g <gene annotation file>'
		exit()

	parser = argparse.ArgumentParser(prog = sys.argv[0], description='Analyze')
	parser.add_argument('-d', '--data_dir', help='data directory', required=True)
	parser.add_argument('-g', '--genes', help='gene annotation file', required=True)
	args = parser.parse_args()
	data_dir = args.data_dir
	genes_fn = args.genes
	
	if not os.path.isfile(genes_fn):
		print 'No such file', genes_fn
		sys.exit(1)
	
	(genes_data, exon_dict, interval_trees) = parse_genes_file(genes_fn)
	print 'Total genes numer', len(genes_data)

	data_dir_list = [os.path.join(data_dir, d) for d in os.listdir(data_dir)]
	for d in data_dir_list:
		print 'processing', os.path.basename(d)
		in_fn = os.path.join(d, os.path.basename(d) + '_PSI_average.txt')
		(pos, tumor_setd2_broken_num, tumor_setd2_broken, tumor_num, tumor, normal_num, normal) = read_data(in_fn)
		exon_to_genes = exon_to_gene_names(exon_dict, interval_trees, pos)

		(gene_beg_dist_bp, gene_beg_dist_perc, gene_names) = get_gene_dist(exon_to_genes, genes_data, pos)

		gene_expression_fn = None
		(tumor_setd2_gene_expr, tumor_gene_expr, normal_gene_expr) = (None, None, None)
		for elem in os.listdir(d):
			if 'gene_expression' in elem:
				gene_expression_fn = os.path.join(d, elem)
		if gene_expression_fn:
			for elem in os.listdir(d):
				if 'Mutation_Packager_Oncotated_Raw_Calls' in elem:
					maf_dir = os.path.join(d, elem)
			if not os.path.exists(maf_dir):
				print 'no such dir:', maf_dir
			(gene_expr, sample_classification, tumor_setd2_gene_expr, tumor_gene_expr, normal_gene_expr) = get_gene_expression(gene_expression_fn, exon_to_genes, pos, maf_dir)

		output_expression_and_dist_file(os.path.join(d, os.path.basename(d) + '_expression_and_dist.txt'), pos, gene_beg_dist_bp, gene_beg_dist_perc, gene_expression_fn, tumor_setd2_gene_expr, tumor_gene_expr, normal_gene_expr, gene_names)


#!/usr/bin/env python
# -*- coding: utf-8 -*- 
import os
import sys
import argparse
import numpy
from enum import Enum, unique
import random

class Data_frame:
	def __init__(self, description, psi):
		(chr_name, begin1, end1, begin2, end2, strand) = description.split(':')
		self.chr_name = chr_name
		self.beg1 = int(begin1)
		self.end1 = int(end1)
		self.beg2 = int(begin2)
		self.end2 = int(end2)
		self.strand = strand
		self.psi = [float(i) for i in psi]

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

class Sample_type(Enum):
	unknown = -1
	norma = 0
	tumor_wild_type = 1
	tumor_unclassified = 2
	tumor_mutant = 3

class Categorized_data_frame:
	def __init__(self, data, samples_num):
		self.data = data
		self.samples_num = samples_num
		self.current_arr_to_call_mean = []
		
	def count_mean(self):
		if len(self.current_arr_to_call_mean) != 0:
			self.data.append(numpy.nanmean(self.current_arr_to_call_mean))
		else:
			self.data.append(float('nan'))
		self.current_arr_to_call_mean = []

def get_mutation_impact_dict():
	# information source: http://www.ensembl.org/info/genome/variation/predicted_data.html#consequences
	# http://snpeff.sourceforge.net/SnpEff_manual.html
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
	mutation_impact_dict['3_prime_UTR_variant'] = Impact.modifier
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
	mutation_impact_dict['initiator_codon_variant'] = Impact.low
	return mutation_impact_dict

def find_mut_gene_mutation_impact(s_file_name, mutation_impact_dict, mut_gene):
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
	setd2_impacts = [Impact.no]
	for line in maf_file:
		if line.startswith(mut_gene.upper() + '\t'):
			mutation_type = line.split('\t')[mutation_type_index]
			if not mutation_impact_dict.has_key(mutation_type):
				print 'Unknown mutation type:', mutation_type
				continue
			mut_gene_mutation_impact = mutation_impact_dict[mutation_type]
			impacts.append(mut_gene_mutation_impact)
		if line.startswith('SETD2\t'):
			mutation_type = line.split('\t')[mutation_type_index]
			if not mutation_impact_dict.has_key(mutation_type):
				print 'Unknown mutation type:', mutation_type
				continue
			mut_gene_mutation_impact = mutation_impact_dict[mutation_type]
			setd2_impacts.append(mut_gene_mutation_impact)
	mut_gene_impact = max(impacts)
	setd2_impact = max(setd2_impacts)
	if mut_gene.upper() != 'SETD2':
		if setd2_impact != Impact.no:
			mut_gene_impact = Impact.unknown # mutation both in mut_gene and setd2
	return mut_gene_impact

def output_mut_gene_mutation_rates(maf_dir, out_file, mut_gene):
	mutation_impact_dict = get_mutation_impact_dict()
	sample_list = [os.path.join(maf_dir, el) for el in os.listdir(maf_dir)]
	out_f = open(out_file, 'w')
	for s_file_name in sample_list:
		if not '.hg19.oncotator.hugo_entrez_remapped.maf.txt' in (s_file_name):
			print 'Not a maf file', s_file_name
			continue
		sample = os.path.basename(s_file_name).split('.')[0]
		mutation_impact = find_mut_gene_mutation_impact(s_file_name, mutation_impact_dict, mut_gene)
		out_f.write(sample + '\t' + str(mutation_impact) + '\n')
	out_f.close()

def get_mut_gene_mutation_rate(maf_dir, sample_names, mut_gene):
	gene_mutation_rate = {}
	mutation_impact_dict = get_mutation_impact_dict()
	for sample in sample_names:
		s_file_name = os.path.join(maf_dir, sample[:15] + '.hg19.oncotator.hugo_entrez_remapped.maf.txt')
		if not os.path.exists(s_file_name):
			continue
		mutation_impact = find_mut_gene_mutation_impact(s_file_name, mutation_impact_dict, mut_gene)
		gene_mutation_rate[sample] = mutation_impact
	return gene_mutation_rate

def process_psi_file(psi_filename, maf_dir, mut_gene):
	def get_sample_type(sample_type, sample_name, gene_mutation_rates):
		if (sample_type >= 1) and (sample_type <= 9): # tumor
			if gene_mutation_rates.has_key(sample_name):
				if gene_mutation_rates[sample_name] == Impact.high:
					return Sample_type.tumor_mutant
				elif gene_mutation_rates[sample_name] == Impact.no:
					return Sample_type.tumor_wild_type
				else:
					return Sample_type.tumor_unclassified
		elif (sample_type >= 10) and (sample_type <= 19): # norma
			return Sample_type.norma
		else:
			return Sample_type.unknown

	inf = open(psi_filename)
	sample_names = [s.strip() for s in inf.readline().split()[2:]]
	gene_mutation_rates = get_mut_gene_mutation_rate(maf_dir, sample_names, mut_gene)
	# barcodes https://wiki.nci.nih.gov/display/TCGA/TCGA+Barcode
	# code tables: https://tcga-data.nci.nih.gov/datareports/codeTablesReport.htm?codeTable=sample%20type
	sample_type = [int(s.split('-')[3][:2]) for s in sample_names]
	data = []
	pos = []
	for line in inf:
		data.append(Data_frame(line.split()[0], line.split()[1:]))
		pos.append(line.split()[0])
	inf.close()
	types = []
	for i in xrange(len(sample_type)):
		types.append(get_sample_type(sample_type[i], sample_names[i], gene_mutation_rates))
	return (pos, data, types)

def compute(psi_filename, maf_dir, mut_gene, filter_for_equal_size = False):
	(pos, data, sample_type) = process_psi_file(psi_filename, maf_dir, mut_gene)
	if filter_for_equal_size:
		indexes = {Sample_type.norma : [], Sample_type.tumor_wild_type : [], Sample_type.tumor_mutant : []}
		for i in xrange(len(sample_type)):
			if sample_type[i] in indexes.keys():
				indexes[sample_type[i]].append(i)
		equal_size = min(len(indexes[Sample_type.tumor_wild_type]), len(indexes[Sample_type.tumor_mutant]))
		random.shuffle(indexes[Sample_type.tumor_wild_type])
		random.shuffle(indexes[Sample_type.tumor_mutant])
		indexes[Sample_type.tumor_wild_type] = indexes[Sample_type.tumor_wild_type][:equal_size]
		indexes[Sample_type.tumor_mutant] = indexes[Sample_type.tumor_mutant][:equal_size]
		equal_size_sample_type = []
		for i in xrange(len(sample_type)):
			for key in indexes.iterkeys():
				if i in indexes[key]:
					equal_size_sample_type.append(sample_type[i])
		for elem in data:
			equal_size_psi = []
			for i in xrange(len(sample_type)):
				for key in indexes.iterkeys():
					if i in indexes[key]:
						equal_size_psi.append(elem.psi[i])
			elem.psi = equal_size_psi
		sample_type = equal_size_sample_type

	categorized_psi = {Sample_type.norma : Categorized_data_frame([], 0), Sample_type.tumor_wild_type : Categorized_data_frame([], 0), Sample_type.tumor_mutant : Categorized_data_frame([], 0)}
	for cur_type in sample_type:
		if categorized_psi.has_key(cur_type):
			categorized_psi[cur_type].samples_num += 1
	for elem in data:
		for i in xrange(len(sample_type)):
			cur_type = sample_type[i]
			if categorized_psi.has_key(cur_type) and ~numpy.isnan(elem.psi[i]):
				categorized_psi[cur_type].current_arr_to_call_mean.append(elem.psi[i])
		for (cur_type, cur_psi_data) in categorized_psi.iteritems():
			cur_psi_data.count_mean()
	return (pos, categorized_psi)

def output_sample_avarage_arrays(pos, categorized_psi, out_fn):
	out_f = open(out_fn, 'w')
	out_f.write('Pos\tPSI_tumor_mutant:' + str(categorized_psi[Sample_type.tumor_mutant].samples_num) + '\tPSI_tumor_wild_type:' + str(categorized_psi[Sample_type.tumor_wild_type].samples_num) + '\tPSI_norma:' + str(categorized_psi[Sample_type.norma].samples_num) + '\n')
	for i in xrange(len(pos)):
		out_f.write(pos[i] + '\t' + str(categorized_psi[Sample_type.tumor_mutant].data[i]) + '\t' + str(categorized_psi[Sample_type.tumor_wild_type].data[i]) + '\t' + str(categorized_psi[Sample_type.norma].data[i]) + '\n')
	out_f.close()

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print 'Usage:', sys.argv[0], '-d <data directory> -c <computations directory> -m <mutant gene name>'
		exit()

	parser = argparse.ArgumentParser(prog = sys.argv[0], description='Categorize PSI to normal, tumor with mut_gene mutation and tumor wild type for every exon')
	parser.add_argument('-d', '--data_dir', help='data directory', required=True)
	parser.add_argument('-c', '--comp_dir', help='computations directory', required=True)
	parser.add_argument('-m', '--mut_gene', help='mutatn gene name', required=True)
	args = parser.parse_args()
	data_dir = args.data_dir
	comp_dir = args.comp_dir
	mutant_gene = args.mut_gene

	if not os.path.isdir(data_dir):
		print >> sys.stderr, 'Not a directory ' + data_dir
		sys.exit(1)

	random.seed(1)

	data_dir_list = [os.path.join(data_dir, d) for d in os.listdir(data_dir)]
	for d in data_dir_list:
		print 'processing', os.path.basename(d)
		d_comp = os.path.join(comp_dir, os.path.basename(d))
		if not os.path.exists(os.path.join(d_comp, mutant_gene)):
			os.mkdir(os.path.join(d_comp, mutant_gene))

		for elem in os.listdir(d):
			if 'Mutation_Packager_Oncotated_Raw_Calls' in elem:
				maf_dir = os.path.join(d, elem)
		if not os.path.exists(maf_dir):
			print 'no such dir:', maf_dir

		sample_types_fn = os.path.join(os.path.join(d_comp, mutant_gene), os.path.basename(d) + '_' + mutant_gene + '_mutation_impact_for_samples.txt')
		output_mut_gene_mutation_rates(maf_dir, sample_types_fn, mutant_gene)

		psi_filename = os.path.join(d_comp, os.path.basename(d) + '_PSI.txt')
		if not os.path.exists(psi_filename):
			print 'no such file:', psi_filename
			continue
		(pos, categorized_psi) = compute(psi_filename, maf_dir, mutant_gene, True)
		out_fn = os.path.join(os.path.join(d_comp, mutant_gene), os.path.basename(d) + '_PSI_averaged_by_' + mutant_gene + '.txt')
		output_sample_avarage_arrays(pos, categorized_psi, out_fn)


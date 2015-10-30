#!/usr/bin/env python
# -*- coding: utf-8 -*- 
import os
import sys
import argparse
import numpy
from enum import Enum, unique

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

def compute(psi_filename, maf_dir):
	inf = open(psi_filename)
	sample_names = [s.strip() for s in inf.readline().split()[2:]]
	setd2_mutation_rates = get_setd2_mutation_rate(maf_dir, sample_names)
	# barcodes https://wiki.nci.nih.gov/display/TCGA/TCGA+Barcode
	# code tables: https://tcga-data.nci.nih.gov/datareports/codeTablesReport.htm?codeTable=sample%20type
	sample_type = [int(s.split('-')[3][:2]) for s in sample_names]

	data = []
	pos = []
	for line in inf:
		data.append(Data_frame(line.split()[0], line.split()[1:]))
		pos.append(line.split()[0])
	inf.close()

	tumor_setd2_broken_num = 0
	tumor_num = 0
	normal_num = 0

	tumor_setd2_broken_psi_spliced_average = []
	tumor_psi_spliced_average = []
	normal_psi_spliced_average = []
	for elem in data:
		tumor_setd2_broken_psi_cur = []
		tumor_psi_cur = []
		normal_psi_cur = []
		tumor_setd2_broken_num = 0
		tumor_num = 0
		normal_num = 0
		for i in xrange(len(sample_type)):
			if (sample_type[i] >= 1) and (sample_type[i] <= 9): # tumor
				if setd2_mutation_rates.has_key(sample_names[i]):
					if setd2_mutation_rates[sample_names[i]] == Impact.high:
						tumor_setd2_broken_num += 1
						if not numpy.isnan(elem.psi[i]):
							tumor_setd2_broken_psi_cur.append(elem.psi[i])
					elif setd2_mutation_rates[sample_names[i]] == Impact.no:
						tumor_num += 1
						if not numpy.isnan(elem.psi[i]):
							tumor_psi_cur.append(elem.psi[i])
			elif (sample_type[i] >= 10) and (sample_type[i] <= 19): # norma
				normal_num += 1
				if not numpy.isnan(elem.psi[i]):
					normal_psi_cur.append(elem.psi[i])
		if len(tumor_setd2_broken_psi_cur) != 0:
			tumor_setd2_broken_psi_spliced_average.append(numpy.nanmean(tumor_setd2_broken_psi_cur))
		else:
			tumor_setd2_broken_psi_spliced_average.append(float('nan'))
		if len(tumor_psi_cur) != 0:
			tumor_psi_spliced_average.append(numpy.nanmean(tumor_psi_cur))
		else:
			tumor_psi_spliced_average.append(float('nan'))
		if len(normal_psi_cur) != 0:
			normal_psi_spliced_average.append(numpy.nanmean(normal_psi_cur))
		else:
			normal_psi_spliced_average.append(float('nan'))
	return (pos, tumor_setd2_broken_num, tumor_num, normal_num, tumor_setd2_broken_psi_spliced_average, tumor_psi_spliced_average, normal_psi_spliced_average)

def output_sample_avarage_arrays(pos, tumor_setd2_broken_num, tumor_num, normal_num, tumor_setd2_broken_psi_spliced_average, tumor_psi_spliced_average, normal_psi_spliced_average, out_fn):
	out_f = open(out_fn, 'w')
	out_f.write('Pos\tPSI_tumor_setd2:' + str(tumor_setd2_broken_num) + '\tPSI_tumor:' + str(tumor_num) + '\tPSI_norma:' + str(normal_num) + '\n')
	for i in xrange(len(normal_psi_spliced_average)):
		out_f.write(pos[i] + '\t' + str(tumor_setd2_broken_psi_spliced_average[i]) + '\t' + str(tumor_psi_spliced_average[i]) + '\t' + str(normal_psi_spliced_average[i]) + '\n')
	out_f.close()


if __name__ == '__main__':
	if len(sys.argv) == 1:
		print 'Usage:', sys.argv[0], '-d <data directory>'
		exit()

	parser = argparse.ArgumentParser(prog = sys.argv[0], description='Split to normal, tumor with setd2 mutation and tumor without mutation in setd2')
	parser.add_argument('-d', '--data_dir', help='data directory', required=True)
	args = parser.parse_args()
	data_dir = args.data_dir

	if not os.path.isdir(data_dir):
		print >> sys.stderr, 'Not a directory ' + data_dir
		sys.exit(1)

	data_dir_list = [os.path.join(data_dir, d) for d in os.listdir(data_dir)]
	for d in data_dir_list:
		print 'processing', os.path.basename(d)
		for elem in os.listdir(d):
			if 'Mutation_Packager_Oncotated_Raw_Calls' in elem:
				maf_dir = os.path.join(d, elem)
		if not os.path.exists(maf_dir):
			print 'no such dir:', maf_dir

		for elem in os.listdir(d):
			if os.path.isdir(os.path.join(os.path.join(data_dir, d), elem)):
				maf_dir = os.path.abspath(os.path.join(os.path.join(data_dir, d), elem))
		psi_filename = os.path.join(d, os.path.basename(d) + '_PSI.txt')
		if not os.path.exists(psi_filename):
			print 'no such file:', psi_filename
			continue
		
		(pos, tumor_setd2_broken_num, tumor_num, normal_num, tumor_setd2_broken_psi_spliced_average, tumor_psi_spliced_average, normal_psi_spliced_average) = compute(psi_filename, maf_dir)

		out_fn = os.path.join(d, os.path.basename(d) + '_PSI_average.txt')
		
		output_sample_avarage_arrays(pos, tumor_setd2_broken_num, tumor_num, normal_num, tumor_setd2_broken_psi_spliced_average, tumor_psi_spliced_average, normal_psi_spliced_average, out_fn)


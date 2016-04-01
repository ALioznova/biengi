#!/usr/bin/env python
# -*- coding: utf-8 -*- 
import os
import sys
import argparse
import numpy

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print 'Usage:', sys.argv[0], '-t <traffic lights directory> -a <annotation directory> -o <output directory>'
		exit()

	parser = argparse.ArgumentParser(prog = sys.argv[0], description='Collect annotation and add to TL dir')
	parser.add_argument('-t', '--tl_dir', help='traffic lights directory', required=True)
	parser.add_argument('-a', '--ann_dir', help='annotation directory', required=True)
	parser.add_argument('-o', '--out_dir', help='output directory', required=True)
	args = parser.parse_args()
	tl_dir = args.tl_dir
	ann_dir = args.ann_dir
	out_dir = args.out_dir

	if not os.path.isdir(tl_dir):
		print >> sys.stderr, 'Not a directory ' + tl_dir
		sys.exit(1)

	if not os.path.isdir(ann_dir):
		print >> sys.stderr, 'Not a directory ' + ann_dir
		sys.exit(1)

	if not os.path.exists(out_dir):
		os.mkdir(out_dir)

	tl_dir_list = os.listdir(tl_dir)
	ann_dir_list = os.listdir(ann_dir)
	for tl_file in tl_dir_list:
		tlf = open(os.path.join(tl_dir, tl_file))
		tl_data = []
		for line in tlf:
			tl_data.append(line.split())
		tlf.close()
		new_annotation = []
		pos = []
		for i in xrange(len(tl_data)):
			new_annotation.append(tl_data[i][3])
			pos.append(tl_data[i][1])
		for ann_dir_elem in ann_dir_list:
			if ann_dir_elem.startswith(tl_file.split('_')[0] + '_tl'):
				annotation_name = ann_dir_elem.split('_')[2]
				anf = open(os.path.join(ann_dir, ann_dir_elem))
				cur_line_num = 0
				for line in anf:
					cur_pos = line.split()[0]
					assert cur_pos == pos[cur_line_num], "cur_pos " + cur_pos + " pos " + pos[cur_line_num] + " chr " + tl_file.split('_')[0] + " " + annotation_name
					if line.split()[1] == '-':
						pass
					elif line.split()[1] == '+':
						new_annotation[cur_line_num] += ',' + annotation_name
					elif len(line.split()) == 2:
						new_annotation[cur_line_num] += ',' + annotation_name + '=' + line.split()[1]
					else:
						cur_data = line.split()[1:]
						cur_data_f = [float(e) for e in cur_data]
						cur_data_f = cur_data_f[len(cur_data_f)/2-1: len(cur_data_f)/2+1]
						cur_mean = str(numpy.mean(cur_data_f))
						new_annotation[cur_line_num] += ',' + annotation_name + '=' + cur_mean
					cur_line_num += 1
				anf.close()
		for i in xrange(len(tl_data)):
			tl_data[i][3] = new_annotation[i]
		outf = open(os.path.join(out_dir, tl_file), 'w')
		for line in tl_data:
			outf.write(line[0])
			for elem in line[1:]:
				outf.write('\t' + elem)
			outf.write('\n')
		outf.close



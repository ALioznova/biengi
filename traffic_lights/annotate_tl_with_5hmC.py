#!/usr/bin/env python

lift = open('/mnt/ce4404a9-911c-4f95-9384-5af06f90e465/Projects/cpgTrafficLights/data/annotation/methylation_hg38_pos.bed')
lift_dict = {}
for line in lift:
	(chr_name, start, end, pid) = line.split()
	lift_dict[pid] = (chr_name, start, end)
lift.close()

import numpy as np

fout = open('/mnt/ce4404a9-911c-4f95-9384-5af06f90e465/Projects/cpgTrafficLights/data/annotation/5hmC_hg38.txt', 'w')
data = open('/mnt/ce4404a9-911c-4f95-9384-5af06f90e465/Projects/cpgTrafficLights/data/annotation/GSE63179_series_matrix.txt')
not_lifted = 0
nulls = 0
for line in data:
	if line.startswith("!"):
		continue
	if line.startswith("\n"):
		continue
	if line.startswith("\"ID_REF\""):
		continue
#ID_REF GSM1543269 GSM1543270 GSM1543271 GSM1543272 GSM1543273 GSM1543274 GSM1543275 GSM1543276
	(id_ref, bs1, oxbs3, bs2, oxbs4, bs3, bs4, oxbs1, oxbs2) = line.split()
	id_ref = id_ref[1:-1]
	if not lift_dict.has_key(id_ref):
		not_lifted +=1 
		continue
	(chr_name, start, end) = lift_dict[id_ref]
	a5hmC = [float(e) for e in (oxbs1, oxbs2, oxbs3, oxbs4) if e != 'null']
	if len(a5hmC) == 0:
		nulls += 1
		continue
	fout.write(chr_name + '\t' + start + '\t' + end + '\t' + str(np.mean(a5hmC)) + '\t'  + id_ref + '\n')
fout.close()
data.close()
print 'not_lifted', not_lifted
print 'nulls', nulls

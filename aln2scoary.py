#!/usr/bin/env python

# aln2scoary - Finding association of point mutations and phenotypic traits.
# Copyright (C) 2019 Ji Zhang

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import argparse
import re
import os
import pandas as pd
import time
parser = argparse.ArgumentParser()
parser.add_argument("-a", "--alignment", help="the file name of the sequence alignment file")
parser.add_argument("-t", "--traits", help="the file name of the traits file")
args = parser.parse_args()
ALN = args.alignment
A_TSV = ALN + '.tsv1.tmp'
A_TSV2 = ALN + '.tsv'
A_CSV = ALN + '.csv'
acsv_file = open(A_CSV, 'w')
###########################
seq_collecting = 'OFF'
seq = ''
seqs = {}
atsv_file = open(A_TSV, 'w')
with open(ALN, 'r') as file_aln:
	lines = file_aln.readlines()
	for line in lines:
		if re.match('>', line):
			if seq_collecting == 'ON':
				seq = re.sub(r'^-', 'N', seq)
				while re.match('^N+-', seq):
					seq = re.sub(r'N-', 'NN', seq, 1)
				seq_reversed = seq[::-1]
				seq_reversed = re.sub(r'^-', 'N', seq_reversed)
				while re.match('^N+-', seq_reversed):
					seq_reversed = re.sub(r'N-', 'NN', seq_reversed, 1)
				seq = seq_reversed[::-1] # slice syntax [begin:end:step]
				seqs[head] = seq
				seq = ''
			head = re.sub(r'[\n>]', '', line)
		else:
			seq_collecting = 'ON'
			seq_extention = re.sub(r'\s', '', line)
			seq = seq + seq_extention
seqs[head] = seq # don't forget to process the last record

for head, seq in seqs.items():
	bases_list = list(seq)
	bases_str = '\t'.join(bases_list)
#	line = head + ',' + str(len(seq)) + ',' + seq + ','+ bases_str + '\n'
	line = head + '\t'+ bases_str + '\n'
	atsv_file.write(line)
atsv_file.close()

df = pd.read_csv(A_TSV, sep = '\t', header = None)
df = df.T
df.to_csv(A_TSV2, sep = '\t', header=False, index=True)

with open(A_TSV2, 'r') as file_atsv:
	lines = file_atsv.readlines()
	for line in lines:
		line = line.rstrip()
		if re.match('^0\t', line):
			line = re.sub(r'^0\t', '', line)
			isolates = line.split('\t')
			isolates_n = len(isolates)
			title = "Gene,Non-unique Gene name,Annotation,No. isolates,No. sequences,Avg sequences per isolate,Genome Fragment,Order within Fragment,Accessory Fragment,Accessory Order with Fragment,QC,Min group size nuc,Max group size nuc,Avg group size nuc"
			for isolate in isolates:
				title = title + ',' + isolate
			title = title + '\n'
			acsv_file.write(title)
		else:
			poz = line
			poz = re.sub(r'\t.*', '', line)
			bases = line
			bases = re.sub(r'^[0-9]+\t|[Nn]', '', line)
			bases = bases.split('\t')
			bases_filtered = [x for x in bases if x]
			bases_len = len(bases_filtered)
			if bases_len == isolates_n:
				variation = list(set(bases))
				variation_filtered = [x for x in variation if x]
				var_len = len(variation_filtered)
				var = ''.join(variation_filtered)
				if var_len > 1:
					for base in variation_filtered:
						Gene = poz + base
						Annotation = poz + var
						acsv_line = Gene + ',,' + Annotation + ',,,,,,,,,,,'
						counter = -1
						for base2 in bases:
							counter = counter + 1
							if base2 == base:
								isolate = isolates[counter]
								acsv_line = acsv_line + ',' + isolate + '_' + poz + base
							else:
								acsv_line = acsv_line + ','
						acsv_line = acsv_line + '\n'
						acsv_file.write(acsv_line)
			else:
				continue
acsv_file.close()
###########################
TRAITS = args.traits
T_TSV = TRAITS + '.tsv2.tmp'
T_TSV2 = TRAITS + '.tsv3.tmp'
T_CSV = TRAITS + '.csv'

df = pd.read_csv(TRAITS, sep = '\t', header = None)
df = df.T
df.to_csv(T_TSV, sep = '\t', header=False, index=False)

T_tsv2 = open(T_TSV2, 'w')
counter_line = 0
with open(T_TSV, 'r') as file_ttsv:
	lines = file_ttsv.readlines()
	for line in lines:
		line = line.rstrip()
		counter_line = counter_line + 1
		if counter_line == 1:
			isolates = line.split('\t')
			line_out = 'Name' + '\t' + line + '\n'
			T_tsv2.write(line_out)
		else:
			traits = line.split('\t')
			variation = list(set(traits))
			variation_filtered = [x for x in variation if x]
			for trait in variation_filtered:
				line_out = trait
				counter_isolate = -1
				for isolate in isolates:
					counter_isolate = counter_isolate + 1
					trait2 = traits[counter_isolate]
					if trait == trait2:
						line_out = line_out + '\t' + '1'
					else:
						line_out = line_out + '\t' + '0'
				line_out = line_out + '\n'
				T_tsv2.write(line_out)
T_tsv2.close()

df = pd.read_csv(T_TSV2, sep = '\t', header = None)
df = df.T
df.to_csv(T_CSV, sep = ',', header=False, index=False)
###########################
T0 = time.time()
T0 = str(T0)
output_dir = 'output_' + T0
os.system('mkdir ' + output_dir)
os.system('mv ' + A_TSV2 + ' ' + output_dir + '/full_alignment.tsv')
os.system('mv ' + A_CSV + ' ' + output_dir + '/snp_alignment.csv')
os.system('mv ' + T_CSV + ' ' + output_dir + '/traits.csv')
os.system('rm -f *.tmp')
os.chdir(output_dir)


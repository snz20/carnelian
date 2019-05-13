#!/usr/bin/env python

'''
Runs the Carnelian framework for functional binning of metagenomic sequences
using LDPC k-mer hash based features. Based of the paper "Improved functional
binning of metagenomic sequences in clinical studies" by Sumaiya Nazeen and
Bonnie Berger.

This is written by Sumaiya Nazeen <nazeen@mit.edu> and is based off of the
python wrapper of Opal tool originally written by Yun William Yu.

Carnelian uses an adapted version of the implementation of the metagenomic
binning from the source code of K. Vervier, P. Mahe, M. Tournoud, J.-B.
Veyrieras, and J.-P. Vert. Large-scale Machine Learning for Metagenomics
Sequence Classification, Technical report HAL-01151453, May, 2015.
This code is included in the util/ directory, with modifications to enable
using Carnelian's (originally used in Opal framework) Gallagher code based
hashes in util/ldpc.py.

The code from Verview, et al, requires the Genetic Data Analysis Library,
which we have included a copy of under util/ext/ for ease of installation.

This pipeline depends on Python scikit-learn, numpy, pandas, and on Vowpal Wabbit.
Vowpal Wabbit must be properly installed in the system path. It also assumes a working
R installation (version >= 3.3.2) to be available.
'''


# import required packages
from __future__ import print_function
__version__ = "1.0.0"

import os
from os import path, mkdir
from os.path import isdir
import glob
import argparse
import sys
import subprocess
import random
import time
import threading
import pandas as pd
import numpy as np
import csv
import shutil
import operator
from shutil import copyfile
from multiprocessing.dummy import Pool
from sklearn.metrics import precision_score, recall_score
from datetime import datetime
from collections import Counter
from Bio import SeqIO

script_loc = os.path.realpath(__file__)
sys.path.append(os.path.join(os.path.dirname(script_loc),'util'))
import ldpc
import sequtil

# Setting up environment variables
my_env = os.environ.copy()
my_env["PATH"]=(
	os.path.join(os.path.dirname(script_loc),'util') + ":" +
	os.path.join(os.path.dirname(script_loc),'util','ext','gdl-1.1','GDL','bin') + ":" +
	os.path.join(os.path.dirname(script_loc),'util','ext','gdl-1.1','GDL','include') + ":" +
	my_env.get("PATH", ""))
my_env["LD_LIBRARY_PATH"]=(
	os.path.join(os.path.dirname(script_loc),'util','ext','gdl-1.1','GDL','lib') + ":" +
	my_env.get("LD_LIBRARY_PATH", ""))

# Utility classes and functions
class ArgClass:
	def __init__(self, *args, **kwargs):
			self.args = args
			self.kwargs = kwargs


class LabelProb:
        def __init__(self, label, prob):
                self.label = label
                self.prob = prob
        def __getitem__(self,label):
                return self.label
        def __getitem__(self,prob):
                return self.prob
        def __setitem__(self,label,value):
                self.label = value
        def __setitem__(self,prob,value):
                self.prob = value
        def __lt__(self,other):
                return self.prob < other.prob
        def __le__(self,other):
                return self.prob <= other.prob
        def __gt__(self,other):
                return self.prob > other.prob
        def __ge__(self,other):
                return self.prob >= other.prob
        def __eq__(self,other):
                return self.prob == other.prob

def eprint(*args, **kwargs):
	print(*args, file=sys.stderr, **kwargs)

def unique_lines(file):
	'''gets number of unique lines in file'''
	seen = set()
	with open(file) as f:
		for line in f:
			seen.add(line)
	return len(seen)

def safe_makedirs(directory):
	if not os.path.exists(directory):
		os.makedirs(directory)
	else:
		print("Directory already exists!!")
		return 1
	return 0

def get_fasta_and_label(directory):
        '''finds the 'first' fasta file in directory, and returns a tuple with
        it and the matching named label file in the directory if both exist'''
        try:
                fasta = glob.glob(directory + "/*.fasta")[0]
        except IndexError:
                raise RuntimeError("Could not find fasta file in:" + directory)
        labels = os.path.splitext(fasta)[0] + ".label"
        if not os.path.isfile(labels):
                print("Couldn't find matching labels ... creating labels ...")
                name = os.path.splitext(os.path.basename(fasta))[0]
                path = os.path.join(directory,name+'.label')
                cmd = "grep '>' "+ fasta + " > " + path
                print(cmd)
                os.system(cmd)
                labels = glob.glob(directory + "/*.label")[0]
        return [fasta, labels]

def merge_dico(old_dico,new_dico):
        ''' merges old dictionary and new dictionary files while retraining'''
        dold = [line.strip() for line in open(old_dico)]
        dnew = [line.strip() for line in open(new_dico)]
        m = {}
        cur_max = -1
        for l in dold:
                x = l.split('\t')
                key = x[0]
                val = int(x[1])
                m[key] = val
                if val > cur_max:
                        cur_max = val
        for l in dnew:
                x = l.split('\t')
                if x[0] not in m.keys():
                        cur_max += 1
                        m[x[0]] = cur_max
        out_file = open(new_dico,'w')
        for k in m.keys():
                s = k + '\t' + str(m[k]) + '\n'
                out_file.write(s)
        out_file.close()

def parse_extra(parser, namespace):
	namespaces = []
	extra = namespace.extra
	while extra:
		n = parser.parse_args(extra)
		extra = n.extra
		namespaces.append(n)
	return namespaces

def extract_column_two(infile, outfile):
	"""cut -f2 infile > outfile"""
	with open(infile, 'r') as inf:
		with open(outfile, 'w') as outf:
			for line in inf:
				parts = line.split()
				if len(parts) > 1:
					print(parts[1], file=outf)
				else:
					print('',file=outf)

def vw_class_to_label(inputfile, dicofile, outputfile):
	'''Converts vw IDs in a newline delimited list (inputfile) to
	outputfile using the mapping specified in dicofile'''
	dico = {}
	with open(dicofile, "r") as fin:
		for line in fin:
			txid, vwid = line.strip().split()[:2]
			dico[vwid] = txid
	predout = open(outputfile, "w")
	with open(inputfile, "r") as fin:
		for line in fin:
			#print(dico[line.strip()])
			predout.write("%s\n"%(dico[str(int(float(line.strip().split(' ')[0])))]))
	predout.close()

def vw_class_to_label2(inputfile, dicofile, outputfile, cutoff):
        '''Converts vw IDs and probabilities in a newline delimited list (inputfile) to outputfile using the mapping specified in dicofile'''
        dico = {}
        with open(dicofile, "r") as fin:
                for line in fin:
                        txid, vwid = line.strip().split()[:2]
                        dico[vwid] = txid
        predout = open(outputfile, "w")
	probfile = os.path.join(os.path.dirname(outputfile),'label-probabilities.tsv')
        probout = open(probfile,'w')
        i = 1
        with open(inputfile, "r") as fin:
                for line in fin:
                        x = line.strip().split(' ')
                        lp_arr = []
                        for k in x:
                                y = k.split(':')
                                lp_arr.append(LabelProb(y[0],float(y[1])))
                        lp_arr.sort(key=operator.attrgetter('prob'),reverse=True)
                        if lp_arr[0].prob >= cutoff:
                                s = dico[lp_arr[0].label] + '\t' + str(lp_arr[0].prob) + '\n'
                                probout.write(s)
                                predout.write(str(i)+'\t' + dico[lp_arr[0].label]+'\n')
                        else:
                                s = dico[lp_arr[0].label] + '\t' + str(lp_arr[0].prob)+ '\t' + dico[lp_arr[1].label] + '\t' + str(lp_arr[1].prob) + '\t' + dico[lp_arr[2].label] + '\t' + str(lp_arr[2].prob) + '\n'
                                probout.write(s)
                        i += 1
        predout.close()
        probout.close()

def majority_voting(L):
        counts = dict((x,L.count(x)) for x in set(L))
        sorted_counts = sorted(counts.items(), key=operator.itemgetter(1), reverse=True)
        if len(sorted_counts) == 1:
                return sorted_counts[0][0]
	else:
                if sorted_counts[0][1] > sorted_counts[1][1]:
			return sorted_counts[0][0]
		else:
			return 'N.N.N.N'

def resolve_label(out_dir, orf_seq_file, orf_label_file, args):
        orf_labels = [line.strip() for line in open(orf_label_file)]
        read_labels = []

	id_list = []
	with open(orf_seq_file, mode='r') as handle:
		for record in SeqIO.parse(handle,'fasta'):
			id_list.append(record.id.split('_')[0])

	m = {}
	for i in xrange(len(orf_labels)):
		k = id_list[i]
		if k not in m.keys():
			m[k] = [orf_labels[i]]
		else:
			m[k].append(orf_labels[i])

	for k in m.keys():
		if len(m[k]) == 1:
			read_labels.append(m[k][0])
		else:
			read_labels.append(majority_voting(m[k]))

        safe_makedirs(out_dir)
        read_label_path = os.path.join(out_dir,"reads.label")
        of = open(read_label_path,'w')
	for l in read_labels:
		of.write(l+'\n')
        of.close()
        return read_label_path

def resolve_label_wrt_ref(out_dir, orf_seq_file, orf_label_file, read_seq_file, args):
        orf_labels = [line.strip() for line in open(orf_label_file)]
	print(orf_labels)
        ref_id_list = []
        with open(read_seq_file, mode='r') as handle:
                for record in SeqIO.parse(handle,'fasta'):
                        ref_id_list.append(record.id)

	read_labels = ['N.N.N.N' for _ in ref_id_list]

	orf_id_list = []
	with open(orf_seq_file, mode='r') as handle:
		for record in SeqIO.parse(handle,'fasta'):
			orf_id_list.append(record.id.split('_')[0])
	from collections import Counter
	z = Counter(orf_id_list)
	for i in xrange(len(orf_id_list)):
		id = orf_id_list[i]
		ind = ref_id_list.index(id)
		if z[id] == 1:
			read_labels[ind] =  orf_labels[i]
		else:
			oi = [j for j,val in enumerate(orf_id_list) if val == id]
			labels = [orf_labels[j] for j in oi]
			read_labels[ind] = majority_voting(labels)
	#m = {}
        #for i in xrange(len(orf_labels)):
        #        k = orf_id_list[i]
        #        if k not in m.keys():
        #                m[k] = [orf_labels[i]]
        #        else:
        #                m[k].append(orf_labels[i])

        #for k in m.keys():
	#	if k in ref_id_list:
	#		ind = ref_id_list.index(k)
        #        	if len(m[k]) == 1:
        #                	read_labels[ind] = m[k][0]
        #        	else:
        #                	read_labels[ind] = majority_voting(m[k])

        safe_makedirs(out_dir)
        read_label_path = os.path.join(out_dir,"reads_wrt_ref.label")
        of = open(read_label_path,'w')
        for l in read_labels:
                of.write(l+'\n')
        of.close()
        return read_label_path

def get_final_model(directory):
	'''gets a 'final' model from a directory. Note, will match the first
	file ending in _final.model'''
	try:
		model = glob.glob(directory + "/*_final.model")[0]
	except IndexError:
		raise RuntimeError("Could not find final model file in: " + directory)
	return model

def getAbundance(pred_labels, gs_file, mean_read_length, out_dir):
        ''' Calculate raw and effective counts from the predicted
        labels for a given sample 

        pred_labels (string): path to a file containing the functional
                              labels predicted by Carnelian
        gs_file (string): must be the name of the file with gold standard label with
                          full path specification
        mean_read_length (int): mean read length in base pairs
        out_dir (string): path to the output directory where abundance files
                          will be written 
        '''

        lines = [line.strip() for line in open(gs_file)]
        labels = [x.split('\t')[0] for x in lines[1:]]
        lengths = [(float(x.split('\t')[1])-(mean_read_length/3.0)+1)/1000  for x in lines[1:]]
	lengths2 = [float(x.split('\t')[1]) for x in lines[1:]]
        labels = list(set(labels))
        flines = [line.strip() for line in open(pred_labels)]
        counts = Counter(flines)
        raw_counts = [0 for _ in labels]
        eff_counts = [0 for _ in labels]
	rpk_counts = [0 for _ in labels]
        for k in counts.keys():
                if k != "N.N.N.N":
                        r = labels.index(k)
                        raw_counts[r] = counts[k]
                        eff_counts[r] = counts[k]/lengths[r]
			rpk_counts[r] = counts[k]*1000000000/(lengths2[r]*len(flines))
        N = sum(eff_counts)
        #print(N)
        eff_counts2 = [1000000*r/N for r in eff_counts]
        of = open(os.path.join(out_dir,'raw_counts.tsv'),'w')
        of2 = open(os.path.join(out_dir,'effective_counts.tsv'),'w')
	of2.write('EC Label\tRPKM\tTPM\n')
        for i in xrange(len(labels)):
                if raw_counts[i] != 0:
                        of.write(labels[i]+'\t'+str(raw_counts[i])+'\n')
                        of2.write(labels[i]+'\t'+str(rpk_counts[r])+'\t'+str(eff_counts2[i])+'\n')
        of.close()
        of2.close()

def createAbundanceMatrix(predict_dir, aa_dir, gs_file):
	'''creates a samples x predicted_labels abundance matrix for downstream
	analyses given all the label prediction files in predict_dir

	predict_dir (string):  must be a path to the input directory containing label
			       predictions (.label and .vw) files for each sample in its 
			       own sub-directory. If paired end reads, both predictions files
			       should be in the corresponding sample's subdirectory.
	aa_dir (string):  must be a path to an output directory where the abundance
			  matrix will be written. It should not be nested within predict_dir.
	gs_file(string): must be the name of the file with gold standard label with
			 full path specification

	'''
	#samples = os.listdir(predict_dir)
	samples = next(os.walk(predict_dir))[1]
	#print(samples)
	samples.sort()
	lines = [line.strip() for line in open(gs_file)]
	labels = [x.split('\t')[0] for x in lines[1:]]
	labels = list(set(labels))
	#print(labels[0:4])	

	ncol = len(samples)
	nrow = len(labels)
	#print(nrow)
	#print(ncol)

	abmat = np.zeros((nrow,ncol))

	for c in xrange(ncol):
		fpath = os.path.join(predict_dir,samples[c])
		filelist = glob.glob(fpath+'/*.label')
		for f in filelist:
			flines = [line.strip() for line in open(f)]
			counts = Counter(flines)
			for k in counts.keys():
				if k != "N.N.N.N":
					r = labels.index(k)
					abmat[r,c] += counts[k]

	df = pd.DataFrame(abmat, index=labels, columns=samples)
	fname = os.path.join(aa_dir,'raw_counts.tsv')
	
	df.to_csv(fname,sep='\t',index=True,header=True)
	return 0

# Carnelian features
#def translateOne(argument):
#	'''Subroutine for translating one sample on one cpu using transeq'''
#	#print("in translate one")
#	os.system('transeq -frame 6 ' + argument)

def translateSeqs(seq_dir, out_dir, fgsp_loc, args):
	'''
	Find genes in the input reads and translate the coding sequences to ORFs using FragGeneScan 
	using n cpus.

	seq_dir (string):  must be a path to a directory with a nucleotide fasta file
	out_dir (string):  must be a path to an output directory where ORFs will be written
	fgsp_loc(string): must be a path to the directory where FragGeneScan is installed

	Unpacking args:
		ncpus (int):  number of cpus to be used to parallelize the translation
	'''
	ncpus = args.ncpus
	#p=Pool(args.ncpus)
	#my_env["PATH"]=(os.path.dirname(fgsp_loc) + ":" + my_env.get("PATH", ""))
	os.environ["PATH"]=(fgsp_loc + ":" + my_env.get("PATH", ""))
	try:
		fpath = os.path.join(seq_dir,'*fasta')
                fasta_file = [x for x in glob.glob(fpath)][0]
                #name_path = [(name, seq_dir + '/' + name) for name in fasta_filelist]
		first_record = SeqIO.parse(fasta_file, "fasta").next()
		if not sequtil.check_if_nucl(str(first_record.seq)):
			print("Could not find nucleotide fasta file in:" + seq_dir)
			return(1)
	except IndexError:
		raise RuntimeError("Could not find fasta file in:" + seq_dir)
	safe_makedirs(out_dir)
	out_file = os.path.join(out_dir, os.path.basename(fasta_file).rsplit('.',1)[0])
	starttime = datetime.now()
	print('''================================================
Predicting ORFs from  nucleotide fasta file
{:%Y-%m-%d %H:%M:%S}'''.format(starttime))
        sys.stdout.flush()

	#cmd = 'FGS+ -s ' + fasta_file + ' -o ' + out_file + ' -w 0 -r ' + os.path.join(fgsp_loc,'train') + ' -t illumina_1 ' + ' -m 20480'
	print(os.path.join(fgsp_loc,'train/complete'))
	cmd = 'run_FragGeneScan.pl -genome=' + fasta_file + ' -out='+out_file +' -train=complete'+ ' -complete=0' 
	if ncpus > 1:
		#cmd += ' -p ' + str(ncpus)
		cmd += ' -thread=' + str(ncpus)

	os.system(cmd)
	outfile = out_file + '.faa'
	os.system('mv ' + outfile + ' ' + out_file + '.fasta')
	print('''------------------------------------------------
Total wall clock runtime (sec): {}
================================================'''.format(
        (datetime.now() - starttime).total_seconds()))
        sys.stdout.flush()
	return(0)

def frag(test_dir, frag_dir, args):
	'''
	Draw fragments of length l from the fasta file found in the test_dir with
	coverage c. Note that there must be a label file of the same basename with
	matching ids for each of the fasta lines.

	test_dir (string):  must be a path to a directory with a single fasta
						and label file
	frag_dir (string):  must be a path to an output directory

	Unpacking args:
		frag_length (int):  length of fragments to be drawn
		coverage (float):   fraction of times each location is to be covered
							by drawn fragments
	'''
	# Unpack args
	frag_length = args.frag_length
	coverage = args.coverage
	# Finish unpacking args

	fasta, labels = get_fasta_and_label(test_dir)
	safe_makedirs(frag_dir)
	fasta_out = os.path.join(frag_dir, 'test.fragments.fasta')
	gi2label_out = os.path.join(frag_dir, 'test.fragments.gi2label')
	label_out = os.path.join(frag_dir, 'test.fragments.label')
	starttime = datetime.now()
	print('''================================================
Drawing fragments
{:%Y-%m-%d %H:%M:%S}
'''.format(starttime) + '''
frag_length = {frag_length}
coverage = {coverage}
------------------------------------------------
Fasta input:	{fasta}
labels input:   {labels}

Fasta output:   {fasta_out}
gi2label output:{gi2label_out}
labels output:  {label_out}'''.format(
	frag_length=frag_length, coverage=coverage, fasta=fasta,
	labels=labels, fasta_out=fasta_out, gi2label_out=gi2label_out,
	label_out=label_out)
	)
	sys.stdout.flush()
	# set seed (for reproducibility)
	seed = 42
	# draw fragments
	subprocess.check_call(["drawfrag",
		"-i", fasta,
		"-t", labels,
		"-l", str(frag_length),
		"-c", str(coverage),
		"-o", fasta_out,
		"-g", gi2label_out,
		"-s", str(seed)],
		env=my_env)

	# extract labels
	extract_column_two(gi2label_out, label_out)
	print('''------------------------------------------------
Total wall clock runtime (sec): {}
================================================'''.format(
	(datetime.now() - starttime).total_seconds()))
	sys.stdout.flush()
	return 0

def train(ref_dir, model_dir, args):
	'''Draws fragments from the fasta file found in ref_dir. Note that
	there must be a label file of the same basename with matching ids for
	each of the fasta lines.

	ref_dir (string):   must be a path to a directory with a single fasta
						and label file
	model_dir (string): must be a path to an output directory

	Unpacking args:
		frag_length (int):	length of fragments to be drawn
		coverage (float):	fraction of times each location is to be covered
					by drawn fragments
		kmer_length (int):	size of k-mers used
		rweight (int):  	how many positions will be randomly chosen in the
					contiguous k-mer (k-mer length should be multiple
					of row_weight)

		num_hash (int):	 number of hashing functions
		num_batches (int):  number of times to run vowpal_wabbit
		num_passes (int):   number of passes within vowpal_wabbit
		precise (flag): if set trained model will store probabilities for labels
	'''
	# Unpack args
	frag_length = args.frag_length
	coverage = args.coverage
	kmer = args.kmer_length
	row_weight = args.rweight
	hierarchical = args.hweight # only comes into play if > 0
	num_hash = args.num_hash
	num_batches = args.num_batches
	num_passes = args.num_passes
	bits = args.bits
	lambda1 = args.lambda1
	lambda2 = args.lambda2
	# Finish unpacking args
	fasta, labels = get_fasta_and_label(ref_dir)
	starttime = datetime.now()

	if kmer % row_weight != 0:
		raise ValueError("Row weight [{}] must divide into k-mer length [{}].".format(row_weight, kmer))
	if (hierarchical > 0):
		if kmer % hierarchical != 0:
			raise ValueError("Hierarchy middle level [{}] must divide into k-mer length [{}].".format(hierarchical, kmer))
		if hierarchical % row_weight != 0:
			raise ValueError("Row weight[{}] must divide into middle hierarchical structure weight [{}].".format(row_weight, hierarchical))

	print(
	'''================================================
Training using Carnelian + vowpal-wabbit
{:%Y-%m-%d %H:%M:%S}
'''.format(starttime) + '''
frag_length = {frag_length}
coverage:	   {coverage}
k-mer length:   {kmer}'''.format(
	frag_length=frag_length,
	coverage=coverage,
	kmer=kmer))
	if hierarchical > 0:
		print('''hierarchical:   {}'''.format(hierarchical))
	print('''row weight:	 {row_weight}
num hashes:	 {num_hash}
num batches:	{num_batches}
num passes:	 {num_passes}
------------------------------------------------
Fasta input:	{fasta}
labels input:   {labels}
------------------------------------------------'''.format(
	row_weight=row_weight,
	num_hash=num_hash,
	num_batches=num_batches,
	num_passes=num_passes,
	fasta=fasta,
	labels=labels)
	)
	sys.stdout.flush()
	num_labels = unique_lines(labels)
	print("Number labels:  {}".format(num_labels))
	sys.stdout.flush()

	safe_makedirs(model_dir)

	# define output "dictionary" : label <--> vw classes
	dico = os.path.join(model_dir, "vw-dico.txt")

	# define model prefix
	model_prefix = os.path.join(model_dir, "vw-model")

	# generate LDPC spaced pattern
	pattern_file = os.path.join(model_dir, "patterns.txt")
	ldpc.ldpc_write(k=kmer, t=row_weight, _m=num_hash, d=pattern_file)

	seed = 42
	for i in range(num_batches):
		seed = seed + 1
		batch_prefix = os.path.join(model_dir, "train.batch-{}".format(i))
		fasta_batch = batch_prefix + ".fasta"
		gi2label_batch = batch_prefix + ".gi2label"
		label_batch = batch_prefix + ".label"

		# draw fragments
		subprocess.check_call(["drawfrag",
			"-i", fasta,
			"-t", labels,
			"-l", str(frag_length),
			"-c", str(coverage),
			"-o", fasta_batch,
			"-g", gi2label_batch,
			"-s", str(seed)],
			env=my_env)
		# extract labels
		extract_column_two(gi2label_batch, label_batch)
		#cherry = [line.strip() for line in open(label_batch)]
		#print(len(cherry))
		print("calling fasta2skm for batch {}".format(i))
		# learn model
		fasta2skm_param_list = ["fasta2skm",
			"-i", fasta_batch,
			"-t", label_batch,
			"-k", str(kmer),
			"-d", dico,
			"-p", pattern_file]
		print("Getting training set ...")
		sys.stdout.flush()
		training_list = subprocess.check_output(
				fasta2skm_param_list, env=my_env).splitlines()
		print("Shuffling training set ...")
		sys.stdout.flush()
		random.shuffle(training_list)
		curr_model = model_prefix + "_batch-{}.model".format(i)
		prev_model = model_prefix + "_batch-{}.model".format(i-1) # May not exist if first run
		vw_param_base = ["vw",
			"--random_seed", str(seed),
			"-f", curr_model,
			"--cache_file", batch_prefix + ".cache",
			"--passes", str(num_passes),
			"--save_resume"]

		if args.precise: 
                        vw_param_base += ["--loss_function=logistic", "--probabilities"]

		vw_param_firstrun = [
			"--oaa", str(num_labels),
			"--bit_precision", str(bits),
			"--l1", str(lambda1),
			"--l2", str(lambda2)]
 		if i > 0:
			vw_param_list = vw_param_base + ["-i", prev_model]
		else:
			vw_param_list = vw_param_base + vw_param_firstrun
		print(vw_param_list)
		sys.stdout.flush()
		vwps = subprocess.Popen(vw_param_list, env=my_env,
				stdin=subprocess.PIPE, stdout=subprocess.PIPE,
				stderr=subprocess.STDOUT)

		gsp = vwps.communicate(input='\n'.join(training_list))
		print(gsp)
		while vwps.poll() is None:
			l = vwps.stdout.readline()
			sys.stdout.write(l)
			sys.stdout.flush()
#		thread.join() # This shouldn't be necessary, but just being safe.

		if i > 0:
			os.remove(prev_model)
		if i == num_batches - 1:
			os.rename(curr_model, model_prefix + "_final.model")
		os.remove(batch_prefix + ".cache")
		os.remove(fasta_batch)
		os.remove(label_batch)
		os.remove(gi2label_batch)
		print("Finished batch {}".format(i))
	print('''------------------------------------------------
Total wall clock runtime (sec): {}
================================================'''.format(
	(datetime.now() - starttime).total_seconds()))
	sys.stdout.flush()
	return 0

def retrain(old_model_dir, new_model_dir, new_examples_dir, args):
        '''Draws fragments from the fasta file found in ref_dir. Note that
        there must be a label file of the same basename with matching ids for
        each of the fasta lines.

        old_model_dir (string):   must be a path to a directory with old vowpal wabbit model
        new_model_dir (string):   must be a path to a directory that will contain the new model
        new_examples_dir (string): must be a path to a directory containing the new training samples and labels 

        Unpacking args:
                frag_length (int):  length of fragments to be drawn
                coverage (float):   fraction of times each location is to be covered
                            by drawn fragments
                kmer_length (int):         size of k-mers used
                row_weight (int):   how many positions will be randomly chosen in the
                            contiguous k-mer (k-mer length should be multiple
                            of row_weight)

                num_hash (int):     number of hashing functions
                num_batches (int):  number of times to run vowpal_wabbit
                num_passes (int):   number of passes within vowpal_wabbit
		precise (flag):	if set trained model will store probabilities for labels
        '''
        frag_length = args.frag_length
        coverage = args.coverage
        kmer = args.kmer_length
        num_batches = args.num_batches
        num_passes = args.num_passes

        fasta, labels = get_fasta_and_label(new_examples_dir)
        starttime = datetime.now()

        print('''================================================
Retraining using Carnelian + vowpal-wabbit
{:%Y-%m-%d %H:%M:%S}
'''.format(starttime) + 
'''num batches:    {num_batches}
num passes:     {num_passes}
------------------------------------------------
Fasta input:    {fasta}
labels input:   {labels}
------------------------------------------------'''.format(
        num_batches=num_batches,
        num_passes=num_passes,
        fasta=fasta,
        labels=labels)
        )
        sys.stdout.flush()
        num_labels = unique_lines(labels)
        print("Number labels:  {}".format(num_labels))
        sys.stdout.flush()

        safe_makedirs(new_model_dir)

        old_dico = os.path.join(old_model_dir,"vw-dico.txt")
        dico = os.path.join(new_model_dir, "vw-dico.txt")

        # define model prefix
        prev_model = os.path.join(old_model_dir,"vw-model_final.model")
        model_prefix = os.path.join(new_model_dir, "vw-model")

        # copy previously used LDPC spaced pattern
        old_pattern_file = os.path.join(old_model_dir, "patterns.txt") 
        pattern_file = os.path.join(new_model_dir, "patterns.txt")
        copyfile(old_pattern_file, pattern_file)

        seed = 42
        for i in range(num_batches):
                seed = seed + 1
                batch_prefix = os.path.join(new_model_dir, "train.batch-{}".format(i))
                fasta_batch = batch_prefix + ".fasta"
                gi2label_batch = batch_prefix + ".gi2label"
                label_batch = batch_prefix + ".label"

                # draw fragments
                subprocess.check_call(["drawfrag",
                "-i", fasta,
                "-t", labels,
                "-l", str(frag_length),
                "-c", str(coverage),
                "-o", fasta_batch,
                "-g", gi2label_batch,
                "-s", str(seed)],
                env=my_env)
                # extract labels
                extract_column_two(gi2label_batch, label_batch)

                # learn model
                fasta2skm_param_list = ["fasta2skm",
                "-i", fasta_batch,
                "-t", label_batch,
                "-k", str(kmer),
                "-d", dico,
                "-p", pattern_file]
                print("Getting new training examples ...")
                sys.stdout.flush()
                training_list = subprocess.check_output(
                        fasta2skm_param_list, env=my_env).splitlines()
                #print(training_list)
                print("Shuffling training set ...")
                sys.stdout.flush()
                random.shuffle(training_list)
                curr_model = model_prefix + "_batch-{}.model".format(i)
                if i > 0:
                        prev_model = model_prefix + "_batch-{}.model".format(i-1) # May not exist if first run
                vw_param_base = ["vw",
                "--random_seed", str(seed),
                "-f", curr_model,
                "--cache_file", batch_prefix + ".cache",
                "--passes", str(num_passes),
                "--save_resume"]
		if args.precise:
			vw_param_base += ["--loss_function=logistic", "--probabilities"]

                vw_param_list = vw_param_base + ["-i", prev_model]
                print(vw_param_list)
                sys.stdout.flush()
                vwps = subprocess.Popen(vw_param_list, env=my_env,
                        stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                        stderr=subprocess.STDOUT)

                gsp = vwps.communicate(input='\n'.join(training_list))
                print(gsp)
                while vwps.poll() is None:
                        l = vwps.stdout.readline()
                        sys.stdout.write(l)
                        sys.stdout.flush()
                #thread.join() # This shouldn't be necessary, but just being safe.

                if i > 0:
                        os.remove(prev_model)
                if i == num_batches - 1:
                        os.rename(curr_model, model_prefix + "_final.model")
                os.remove(batch_prefix + ".cache")
                os.remove(fasta_batch)
                os.remove(label_batch)
                os.remove(gi2label_batch)

        merge_dico(old_dico,dico)
        print('''------------------------------------------------
Total wall clock runtime (sec): {}
================================================'''.format(
        (datetime.now() - starttime).total_seconds()))
        sys.stdout.flush()
        return 0


def predict(model_dir, test_dir, predict_dir, args):
        '''Predicts functional labels for each record in  fasta file found in test_dir
        using the vowpal-wabbit model. Note that there must be a label file of the 
        same basename with matching ids for each of the fasta lines.

        test_dir (string):   must be a path to a directory with a single fasta
                                                and label file
        model_dir (string): must be a path to a directory with a vw model file
        predict_dir (string):output directory of predictions

        Unpacking args:
                kmer_length (int): size of k-mers used
                ncpus (int):    number of cpus to be used
                precise (flag): if set, two prediction files will be generated -- one with read 
                                ids and labels that have probabilities above cutoff, other with 
                                label probabilities for strong predictions and three top possible
                                labels with probabilities for weak predictions
                cutoff (float): user-defined probability cut off for precise mode

        Returns a tuple with (reffile, predicted_labels_file) for easy input
        into evaluate_predictions.
        '''
        # Unpack args
        kmer = args.kmer_length
        ncpus = args.ncpus
	cutoff = 0.0
        if args.precise:
                cutoff = args.cutoff
        try:
                fasta = [os.path.basename(x) for x in glob.glob(os.path.join(test_dir,'*.fasta'))][0]
        except IndexError:
                raise RuntimeError("Could not find fasta file in:" + test_dir)

        safe_makedirs(predict_dir)
        starttime = datetime.now()
        print(
        '''================================================
Predicting using Carnelian + vowpal-wabbit
{:%Y-%m-%d %H:%M:%S}
'''.format(starttime) + '''
k-mer length:   {kmer}
------------------------------------------------
Fasta input:    {fasta}
------------------------------------------------'''.format(
        kmer=kmer,
        fasta=fasta))

        sys.stdout.flush()
        if ncpus > 1:
                p=Pool(ncpus)
                tmp_path = os.path.join(predict_dir, 'tmp')
                safe_makedirs(tmp_path)
                infastapath = os.path.join(test_dir, fasta)
                splits = sequtil.split_fasta2(infastapath, tmp_path, ncpus, my_env)
                #splits = sequtil.split_fasta(infastapath, tmp_path, ncpus)
		print(splits)
                test_dirlist = [os.path.dirname(x) for x in glob.glob(os.path.join(tmp_path,'*/*.fasta'))]
                pred_dirlist = [os.path.join(f,'predict') for f in test_dirlist]
                arglist=[model_dir+'|'+test_dirlist[i]+'|'+pred_dirlist[i]+'|'+str(kmer) for i in range(ncpus)]
                if args.precise:
                        arglist = [x+'|'+str(args.precise)+'|'+str(cutoff) for x in arglist]
                p.map(predict_unpack, arglist)
                filelist = [glob.glob(os.path.join(d,'*.label'))[0] for d in pred_dirlist]
                #print(filelist)
                prefix = sequtil.merge_files2(filelist,predict_dir,fasta,'label', my_env)
		vw_prefix = sequtil.merge_files2(filelist,predict_dir,fasta,'vw', my_env)
		if args.precise:
			filelist = [glob.glob(os.path.join(d,'*.tsv'))[0] for d in pred_dirlist]
			prob_file = sequtil.merge_files2(filelist,predict_dir,fasta,'tsv',my_env)
                shutil.rmtree(tmp_path, ignore_errors=True)
        else:
                prefix = predictOne(model_dir, test_dir, predict_dir, kmer, args.precise, cutoff)
		pred_file = os.path.join(os.path.dirname(prefix),fasta.split('.')[0]+'.label')
		os.system("mv "+prefix+' '+pred_file)
		orig_file = prefix.rsplit('.',1)[0]+'.vw'
		vw_file = os.path.join(os.path.dirname(prefix),fasta.split('.')[0]+'.vw')
		os.system("mv "+orig_file+' '+vw_file)
		if args.precise:
			prob_file = os.path.join(os.path.dirname(prefix),fasta.split('.')[0]+'.tsv')
			orig_file = os.path.join(os.path.dirname(prefix),'label-probabilities.tsv')
			os.system("mv "+orig_file+' '+prob_file)
		prefix = pred_file

        print('''------------------------------------------------
Predicted labels:   {pl}
Total wall clock runtime (sec): {s}
================================================'''.format(
        pl=prefix,
        s=(datetime.now() - starttime).total_seconds()))
        sys.stdout.flush()
        return (prefix.strip())

def predict_unpack(args):
        s = args.split('|')
        model_dir = s[0]
        test_dir = s[1]
        predict_dir = s[2]
        k = int(s[3])
	precise = False
        cutoff=0.0
        if len(s) > 5: 
		precise = bool(s[4])
                cutoff = float(s[5])
        predictOne(model_dir, test_dir, predict_dir, k, precise, cutoff)

def predictOne(model_dir, test_dir, predict_dir, kmer, precise, cutoff=0.0):
        # Don't need to get labels until eval
        #fasta, labels = get_fasta_and_label(test_dir)
        try:
                fasta = glob.glob(test_dir + "/*.fasta")[0]
        except:
                raise RuntimeError("No fasta file found in: " + test_dir)
        model = get_final_model(model_dir)
        dico = os.path.join(model_dir, "vw-dico.txt")
        pattern_file = os.path.join(model_dir, "patterns.txt")
        safe_makedirs(predict_dir)
        prefix = os.path.join(predict_dir, "test.fragments-db")
                # get vw predictions
        fasta2skm_param_list = ["fasta2skm",
                "-i", fasta,
                "-k", str(kmer),
                "-p", pattern_file]
        vw_param_list = ["vw", "-t",
                "-i", model,
                "-p", prefix + ".preds.vw"]
        if precise:
                vw_param_list.append("--probabilities")
        ps = subprocess.Popen(fasta2skm_param_list, env=my_env,
                        stdout=subprocess.PIPE)
        vwps = subprocess.Popen(vw_param_list, env=my_env,
                        stdin=ps.stdout, stdout=subprocess.PIPE,
                        stderr=subprocess.STDOUT)
        while vwps.poll() is None:
                l = vwps.stdout.readline()
                sys.stdout.write(l)
                sys.stdout.flush()
        if cutoff > 0:
                vw_class_to_label2(prefix + '.preds.vw', dico, prefix + '.preds.label', cutoff)
        else:
                vw_class_to_label(prefix + '.preds.vw', dico, prefix + '.preds.label')
        return (prefix + '.preds.label')


def evaluatePreds(ref_file, pred_file):
	'''Evaluates how good a predicted list is compared to a reference gold standard'''
	with open(pred_file, "r") as fin:
		pred = fin.read().splitlines()
	with open(ref_file, "r") as fin:
		ref = fin.read().splitlines()

	correct = [False for _ in ref]
	for i in xrange(len(ref)):
		if ref[i] == pred[i]:
			correct[i] = True

	perf = pd.DataFrame({"pred":pred, "ref":ref, "correct":correct})
	tmp = perf.groupby("ref")
        recall_i = tmp["correct"].agg(np.mean)
        micro_recall = np.mean(correct)
        macro_recall = np.mean(recall_i)
        median_recall = np.median(recall_i)
        tmp = perf.groupby("pred")
        precision_i = tmp["correct"].agg(np.mean)
        micro_precision = micro_recall
        macro_precision = np.mean(precision_i)
        median_precision = np.median(precision_i)
        micro_f1 = 2*micro_precision*micro_recall/(micro_precision+micro_recall)
        f1score_i = 2*recall_i*(precision_i+0.00000001)/(recall_i+precision_i+0.00000001)
        macro_f1 = np.mean(f1score_i)
  
        print("micro recall = {:.2f}".format(micro_recall*100))
        print("macro recall = {:.2f}".format(macro_recall*100))
        print("median recall = {:.2f}".format(median_recall*100))

        print("micro precision = {:.2f}".format(micro_precision*100))
        print("macro precision = {:.2f}".format(macro_precision*100))
        print("median precision = {:.2f}".format(median_precision*100))

        print("micro F1 = {:.2f}".format(micro_f1*100))
        print("macro F1 = {:.2f}".format(macro_f1*100))
	sys.stdout.flush()
	return 0

def evaluatePredsPrecise(reffile, predfile):
        '''Evaluates how good a predicted list is compared to a reference gold standard'''
        with open(predfile, "r") as fin:
                pred_l = fin.read().splitlines()
        with open(reffile, "r") as fin:
                ref = fin.read().splitlines()

        pred = ['NA' for _ in ref]
        for l in pred_l:
                x = l.strip().split('\t')
                id = int(x[0]) - 1
                pred[id] = x[1]

        correct = [False for _ in ref]
        for i in xrange(len(ref)):
                if ref[i] == pred[i]:
                        correct[i] = True

        perf = pd.DataFrame({"pred":pred, "ref":ref, "correct":correct})
        tmp = perf.groupby("ref")
        recall_i = tmp["correct"].agg(np.mean)
        micro_recall = np.mean(correct)
        macro_recall = np.mean(recall_i)
        median_recall = np.median(recall_i)
        tmp = perf.groupby("pred")
        precision_i = tmp["correct"].agg(np.mean)
        micro_precision = micro_recall
        macro_precision = np.mean(precision_i)
        median_precision = np.median(precision_i)
        micro_f1 = 2*micro_precision*micro_recall/(micro_precision+micro_recall)
        f1score_i = 2*recall_i*precision_i/(recall_i+precision_i)
        macro_f1 = np.mean(f1score_i)
  
        print("micro recall = {:.2f}".format(micro_recall*100))
        print("macro recall = {:.2f}".format(macro_recall*100))
        print("median recall = {:.2f}".format(median_recall*100))

        print("micro precision = {:.2f}".format(micro_precision*100))
        print("macro precision = {:.2f}".format(macro_precision*100))
        print("median precision = {:.2f}".format(median_precision*100))

	print("micro F1 = {:.2f}".format(micro_f1*100))
        print("macro F1 = {:.2f}".format(macro_f1*100))
        sys.stdout.flush()
	return 0


def calcAbundance(in_dir, out_dir, mapping_file, gs_file):
	'''
	Performs abundance estimation on the predicttions in the directory in_dir
	and outputs the raw counts and effective counts matrices in directory out_dir. 
	.mapping file is a tab separated file with three columns: <sample_id,group,fraglen>
	
	in_dir (string):  must be a path to an input directory containing the
			  predicted labels for each sample in its own subdirectory
	out_dir (string):  must be a path to an output directory
	mapping_file (string): must be a path to a tab separated file with three 
			       columns: <sample_id,group,fraglen>. The first line 
			       of the file should contain headers.
	gs_file (string):  must be a path to a tab-separated input file containing the gold
			   standard protein labels with average protein lengths under each label
	'''
	safe_makedirs(out_dir)
	starttime = datetime.now()
        print('''================================================
Performing abundance estimation with Carnelian and R
Creating raw counts matrix
{:%Y-%m-%d %H:%M:%S}'''.format(starttime))
	sys.stdout.flush()

	df = createAbundanceMatrix(in_dir, out_dir, gs_file)
	counts_file = os.path.join(out_dir,'raw_counts.tsv')

	print(" "+counts_file)

	command =  counts_file + ' ' + gs_file + ' ' + mapping_file + ' ' + out_dir
	rscript_path = os.path.dirname(script_loc)+'/scripts/abundance_estimation.R'
	print(rscript_path)
	os.system('Rscript '+rscript_path+ ' ' + command)
	print('''------------------------------------------------
Total wall clock runtime (sec): {}
================================================'''.format(
    	(datetime.now() - starttime).total_seconds()))
    	sys.stdout.flush()

	return 0

# Carnelian main function
def main(argv):
	parser = argparse.ArgumentParser(description='Perform functional binning of metagenomic sequences using gene annotations')
	# Shared arguments
	#parser.add_argument('--version', action='version', version='%(prog)s {version}'.format(version=__version__))
	ncpus_arg = ArgClass('-n', dest='ncpus', default=1, help='Number of parallel processors to be used', type=int)
	kmer_arg = ArgClass('-k', dest='kmer_length', default=8, help='length of k-mers used', type=int)
	num_hash_arg = ArgClass('--num_hash', dest='num_hash', default=1, help='number of k-mer hashing functions to get features', type=int)
	frag_length_arg = ArgClass('-l', dest='frag_length', default=16, help='length of fragments to be drawn from fasta', type=int)
	coverage_arg = ArgClass('-c', dest='coverage', default=1.0, help='number/fraction of times a position in a fragment should be covered by the k-mer', type=float)
	hierarchical_arg = ArgClass('--hweight', help='intermediate organization of positions chosen in the k-mer in row_weight; should be a multiple of row_weight and a divisor of k-mer length if set', type=int, default=-1)
	row_weight_arg = ArgClass('--rweight', help='the number of positions that will be randomly chosen in the contiguous k-mer; k-mer length should be a multiple of row_weight', type=int, default=4)
	num_batches_arg = ArgClass('--num-batches', help='Number of times to generate a random batch of training data for VW', type=int, default=1)
	num_passes_arg = ArgClass('--num-passes', help='Number of VW passes in each training batch', type=int, default=1)
	bits_arg = ArgClass('--bits', help='Number of bits used in VW model', type=int, default=31)
	lambda1_arg = ArgClass('--lambda1', help='VW model lambda1 training parameter', type=float, default=0.)
	lambda2_arg = ArgClass('--lambda2', help='VW model lambda2 training parameter', type=float, default=0.)
	type_arg = ArgClass('--type', dest='type', default='nucleotide', help='Sequence type. "prot" if the reference fasta files are for amino acid sequences. "nucl" if nucleotide sequences. In the latter case the nucleotide sequences will be translated. (default:  %(default)s)')
	cutoff_arg = ArgClass('--cutoff', help='Probability cutoff for VW predictions', type=float, default=0.)

	# Subparsers
	subparsers = parser.add_subparsers(help='sub-commands', dest='mode')

	# Translation Args
	parser_translate = subparsers.add_parser('translate', help='Find genes in the input reads.fasta file and generate all possible ORFs', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser_translate.add_argument('seq_dir', help='Input directory for fasta file to be translated')
	parser_translate.add_argument('out_dir', help='Output directory for fasta translated fasta file')
	parser_translate.add_argument('fgsp_loc', help='Path where FragGeneScan is installed')
	parser_translate.add_argument(*ncpus_arg.args, **ncpus_arg.kwargs)

	# Fragmentation Args
	parser_frag = subparsers.add_parser('frag', help='Fragment a fasta file into substrings for training/testing', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser_frag.add_argument('test_dir', help='Input directory for fasta and label files to be fragmented')
	parser_frag.add_argument('frag_dir', help='Output directory for fasta fragments and corresponding labels')
	parser_frag.add_argument(*frag_length_arg.args, **frag_length_arg.kwargs)
	parser_frag.add_argument(*coverage_arg.args, **coverage_arg.kwargs)
	parser_frag.add_argument(*type_arg.args, **type_arg.kwargs)

	# Training Args
	parser_train = subparsers.add_parser('train', help='Train a Vowpal Wabbit model using carnelian hash-based features', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser_train.add_argument('train_dir', help='Input directory for training data')
	parser_train.add_argument('model_dir', help='Output directory for VW model')
	parser_train.add_argument('--precise', help='If set, will train model to output probabilities', action='store_true')
	parser_train.add_argument(*frag_length_arg.args, **frag_length_arg.kwargs)
	parser_train.add_argument(*coverage_arg.args, **coverage_arg.kwargs)
	#parser_train.add_argument(*type_arg.args, **type_arg.kwargs)
	parser_train.add_argument(*kmer_arg.args, **kmer_arg.kwargs)
	parser_train.add_argument(*num_batches_arg.args, **num_batches_arg.kwargs)
	parser_train.add_argument(*num_passes_arg.args, **num_passes_arg.kwargs)
	parser_train.add_argument(*num_hash_arg.args, **num_hash_arg.kwargs)
	parser_train.add_argument(*row_weight_arg.args, **row_weight_arg.kwargs)
	parser_train.add_argument(*hierarchical_arg.args, **hierarchical_arg.kwargs)
	parser_train.add_argument(*bits_arg.args, **bits_arg.kwargs)
	parser_train.add_argument(*lambda1_arg.args, **lambda1_arg.kwargs)
	parser_train.add_argument(*lambda2_arg.args, **lambda2_arg.kwargs)

	parser_retrain = subparsers.add_parser("retrain", help="Incrementally retrain an existing Vowpal Wabbit model using new examples", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser_retrain.add_argument("old_model_dir", help="Input directory for old VW model")
	parser_retrain.add_argument("new_model_dir", help="Output directory for new VW model")
	parser_retrain.add_argument("new_examples_dir", help="Input directory for fasta and label files containing new examples")
	parser_retrain.add_argument('--precise', help='If set, will train model to output probabilities', action='store_true')
	parser_retrain.add_argument(*frag_length_arg.args, **frag_length_arg.kwargs)
	parser_retrain.add_argument(*coverage_arg.args, **coverage_arg.kwargs)
	parser_retrain.add_argument(*kmer_arg.args, **kmer_arg.kwargs)
	parser_retrain.add_argument(*num_batches_arg.args, **num_batches_arg.kwargs)
	parser_retrain.add_argument(*num_passes_arg.args, **num_passes_arg.kwargs)

	# Prediction Args
	parser_predict = subparsers.add_parser('predict', help='Predict functional labels with or without probabilities given a carnelian model', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser_predict.add_argument('model_dir', help='Input directory of VW model obtained from training')
	parser_predict.add_argument('test_dir', help='Input directory containing already fragmented test data')
	parser_predict.add_argument('predict_dir', help='Output directory for predictions')
	parser_predict.add_argument('--precise', help='If set, label probabilities will be generated', action='store_true')
	parser_predict.add_argument(*kmer_arg.args, **kmer_arg.kwargs)
	parser_predict.add_argument(*cutoff_arg.args, **cutoff_arg.kwargs)
	parser_predict.add_argument(*ncpus_arg.args, **ncpus_arg.kwargs)

	# Performance Evaluation Args
	parser_eval = subparsers.add_parser('eval', help='Evaluate quality of predictions given gold standard labels', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser_eval.add_argument('reference_file', help='Gold standard labels')
	parser_eval.add_argument('predicted_labels', help='Predicted labels')
	parser_eval.add_argument('--precise', help='If set, expect predicted labels file to have two tab-separated columns <readID, label>', action='store_true')

	# Abundance Analysis Args
	parser_aa = subparsers.add_parser('abundance', help='Generate abundance estimates of functional terms', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser_aa.add_argument('in_dir', help='Input directory containing labeled reads for each sample in its own sub-directory')
	parser_aa.add_argument('out_dir', help='Output directory with abundance results')
	parser_aa.add_argument('mapping_file', help='Mapping file containing sample ids, case-control labels, and nominal read lengths')
	parser_aa.add_argument('gs_file', help='File containing gold standard functional labels and average gene lengths')

	# Simulate Args
	parser_simulate = subparsers.add_parser('simulate', help='''Run a full pipeline of frag, train, predict, and eval to
determine how good a model is under particular parameter ranges''', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	parser_simulate.add_argument('test_dir', help='Input directory for test data')
	parser_simulate.add_argument('train_dir', help='Input directory for train data')
	parser_simulate.add_argument('out_dir', help='Output directory for all steps')
	parser_simulate.add_argument('--do-not-translate', help='If set, will use test_dir fasta files as is without translating', action='store_true')
	parser_simulate.add_argument('--precise', help='If set, label probabilities will be generated', action='store_true')
	parser_simulate.add_argument('fgsp_loc', help='Path where FragGeneScanPlus is installed')
        parser_simulate.add_argument(*ncpus_arg.args, **ncpus_arg.kwargs)
	parser_simulate.add_argument(*frag_length_arg.args, **frag_length_arg.kwargs)
	parser_simulate.add_argument(*coverage_arg.args, **coverage_arg.kwargs)
	parser_simulate.add_argument(*kmer_arg.args, **kmer_arg.kwargs)
	parser_simulate.add_argument(*num_batches_arg.args, **num_batches_arg.kwargs)
	parser_simulate.add_argument(*num_passes_arg.args, **num_passes_arg.kwargs)
	parser_simulate.add_argument(*num_hash_arg.args, **num_hash_arg.kwargs)
	parser_simulate.add_argument(*row_weight_arg.args, **row_weight_arg.kwargs)
	parser_simulate.add_argument(*hierarchical_arg.args, **hierarchical_arg.kwargs)
	parser_simulate.add_argument(*bits_arg.args, **bits_arg.kwargs)
	parser_simulate.add_argument(*lambda1_arg.args, **lambda1_arg.kwargs)
	parser_simulate.add_argument(*lambda2_arg.args, **lambda2_arg.kwargs)
	parser_simulate.add_argument(*cutoff_arg.args, **cutoff_arg.kwargs)

	#  Annotate Args
	parser_annotate = subparsers.add_parser('annotate', help='''End-to-end annotation of a sample staring from translation, 
prediction, and label resolution''', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser_annotate.add_argument('sample_dir', help='Input directory containing one read file')
        parser_annotate.add_argument('model_dir', help='Input directory for trained model')
        parser_annotate.add_argument('out_dir', help='Output directory for all steps')
        parser_annotate.add_argument('--precise', help='If set, label probabilities will be generated', action='store_true')
        parser_annotate.add_argument('fgsp_loc', help='Path where FragGeneScan is installed')
        parser_annotate.add_argument('--do-not-translate', help='If set, will use test_dir fasta files as is without translating', action='store_true')
        parser_annotate.add_argument(*ncpus_arg.args, **ncpus_arg.kwargs)
        parser_annotate.add_argument(*kmer_arg.args, **kmer_arg.kwargs)
        parser_annotate.add_argument(*cutoff_arg.args, **cutoff_arg.kwargs)

	args=parser.parse_args(argv)
	print(args)
	sys.stdout.flush()

	mode = args.mode

	if (mode == 'simulate'):
		st_time = datetime.now()
		print('Simulation for Performance')
		print("{:%Y-%m-%d %H:%M:%S}".format(st_time))

		output_dir = args.out_dir
		trans_dir = os.path.join(output_dir, '1orf')
		model_dir = os.path.join(output_dir, '2model')
		predict_dir = os.path.join(output_dir, '3predict')

		if args.do_not_translate:
			train(args.train_dir, model_dir, args)
			pf = predict(model_dir, args.test_dir, predict_dir, args)
			_, rf = get_fasta_and_label(args.test_dir)
		else:
			translateSeqs(args.test_dir, trans_dir, args.fgsp_loc, args)
			train(args.train_dir, model_dir, args)
			orf_label = predict(model_dir, trans_dir, predict_dir, args)
			orf_seq = [x for x in glob.glob(os.path.join(trans_dir,"*.fasta"))][0]
			read_seq, rf = get_fasta_and_label(args.test_dir)
			pf = resolve_label_wrt_ref(predict_dir, orf_seq, orf_label, read_seq, args)

		print("Evaluation reference file: " + rf)
		sys.stdout.flush()
		if args.precise:
			v = evaluatePredsPrecise(rf, pf)
		else:
			v = evaluatePreds(rf, pf)
		print("Total full sim wall clock runtime (sec): {}".format((datetime.now() - st_time).total_seconds()))

	elif mode == 'annotate':
		st_time = datetime.now()
                print('Annotating one sample')
                print("{:%Y-%m-%d %H:%M:%S}".format(st_time))

		output_dir = args.out_dir
                trans_dir = os.path.join(output_dir, '1orf')
                predict_dir = os.path.join(output_dir, '2predict')
		model_dir = args.model_dir

		if args.do_not_translate:
			pf = predict(model_dir, args.sample_dir, predict_dir, args)
			getAbundance(pf, gs_file, read_len, predict_dir)
		else:
			translateSeqs(args.sample_dir, trans_dir, args.fgsp_loc, args)
			orf_label = predict(model_dir, trans_dir, predict_dir, args)
                        orf_seq = [x for x in glob.glob(os.path.join(trans_dir,"*.fasta"))][0]

		print("Total full annot wall clock runtime (sec): {}".format((datetime.now() - st_time).total_seconds()))


	elif mode == 'translate':
		translateSeqs(args.seq_dir, args.out_dir, args.fgsp_loc, args)

	elif mode == 'frag':
		frag(args.test_dir, args.frag_dir, args)

	elif mode == 'train':
		train(args.train_dir, args.model_dir, args)

	elif mode == "retrain":
        	retrain(args.old_model_dir, args.new_model_dir, args.new_examples_dir, args)

	elif mode == 'predict':
		predict(args.model_dir, args.test_dir, args.predict_dir, args)

	elif mode == 'eval':
		if args.precise:
			v = evaluatePredsPrecise(args.reference_file, args.predicted_labels)
		else:
			v = evaluatePreds(args.reference_file, args.predicted_labels)

	elif mode == 'abundance':
		calcAbundance(args.in_dir, args.out_dir, args.mapping_file, args.gs_file)

if __name__ == "__main__":
	main(sys.argv[1:])

# import required packages
from __future__ import print_function
__version__ = "0.6.0"

import unittest
import importlib
import os
import filecmp
import tempfile
from os import path, mkdir
from os.path import isdir
import glob
import argparse
import sys
import numpy as np
import shutil
from distutils.spawn import find_executable


import config

script_loc = os.path.realpath('../carnelian.py')
util_path = os.path.join(os.path.dirname(script_loc),'util')
sys.path.append(os.path.dirname(script_loc))
fgsp_loc = os.path.join(util_path,'ext/FragGeneScan')

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

from carnelian import unique_lines
from carnelian import safe_makedirs
from carnelian import merge_dico
from carnelian import extract_column_two
from carnelian import vw_class_to_label
from carnelian import vw_class_to_label2
from carnelian import evaluatePreds
from carnelian import evaluatePredsPrecise
from carnelian import createAbundanceMatrix
from carnelian import majority_voting
from carnelian import resolve_label
from carnelian import resolve_label_wrt_ref
from carnelian import frag
from carnelian import translateSeqs
from carnelian import train
from carnelian import retrain
from parcarn import predict
from carnelian import calcAbundance

import ldpc
import sequtil

class TestCarnelian_Advanced(unittest.TestCase):
	'''
	Test the functions in carnelian.py
	'''

	def initialize(self):
		logging.getLogger('carnelian_advanced_tests').addHandler(logging.NullHandler())

	def test2_carnelian_one_sample(self):
		'''
		Takes a nucleotide fasta file, translates the reads, predicts the ORF labels, 
		and resolves the ORF  labels to get the read labels.
		'''

		input_dir = config.ttest_dir
		model_dir = config.model_dir
		tmp_trans_dir = os.path.join(config.ttest_dir, 'translated')
		tmp_pred_dir = os.path.join(config.ttest_dir, 'pred')
		tmp_resolve_dir = os.path.join(config.ttest_dir, 'resolved')

		a = translateSeqs(input_dir, tmp_trans_dir, fgsp_loc, args)
		self.assertEqual(a,0)
		pep_file = [x for x in glob.glob(os.path.join(tmp_trans_dir,"*.fasta"))][0]
		self.assertTrue(os.path.isfile(pep_file))
                self.assertTrue(os.stat(pep_file).st_size > 0)

		orf_label_file = predict(model_dir, tmp_trans_dir, tmp_pred_dir, args)
		self.assertTrue(os.path.isfile(orf_label_file))
		self.assertTrue(os.stat(orf_label_file).st_size > 0)

		resolve_label(tmp_resolve_dir, pep_file, orf_label_file, args)
		read_label_file = os.path.join(tmp_resolve_dir, "reads.label")
		self.assertTrue(os.path.isfile(read_label_file))
		self.assertTrue(os.stat(orf_label_file).st_size > 0)

		shutil.rmtree(tmp_trans_dir)
		shutil.rmtree(tmp_pred_dir)
		shutil.rmtree(tmp_resolve_dir)

	def test1_carnelian_end_to_end(self):
                '''
                First trains a model. Then for each of the 6 samples, takes  a nucleotide fasta file, 
		translates the reads, predicts the ORF labels, and resolves the ORF  labels to get the read labels.
		Finally generates an effective counts matrix for the samples.
                '''
		ft_dir = config.ftest_dir
		train_dir = config.data_dir
		print(train_dir)
		out_dir = os.path.join(ft_dir, 'tmp')
		safe_makedirs(out_dir)
		model_dir = os.path.join(train_dir,'model')
		samples_dir = config.samples_dir
		list_of_samples = [os.path.join(samples_dir, x) for x in next(os.walk(samples_dir))[1]]

		# train
		expected_filelist = ['patterns.txt','vw-dico.txt','vw-model_final.model']
        	v = train(ft_dir, model_dir, args)
        	self.assertEqual(v,0)
        	out = os.popen('ls '+model_dir).read()
        	model_filelist = out.split()
		print(model_filelist)
		# check if the expected files are generated and the files are non-zero
        	self.assertTrue(set(model_filelist) == set(expected_filelist))
		for f in model_filelist:
			self.assertTrue(os.stat(os.path.join(model_dir,f)).st_size > 0)

		# prediction for each samples
		trans_dir = os.path.join(out_dir, 'translated')
		pred_dir = os.path.join(out_dir, 'predicted')
		resolve_dir = os.path.join(out_dir, 'resolved')

		for sample in list_of_samples:
			input_dir = sample
			print(input_dir)
                	tmp_trans_dir = os.path.join(trans_dir,os.path.basename(sample))
                	tmp_pred_dir = os.path.join(pred_dir,os.path.basename(sample))
                	tmp_resolve_dir = os.path.join(resolve_dir,os.path.basename(sample))

                	a = translateSeqs(input_dir, tmp_trans_dir, fgsp_loc, args)
                	self.assertEqual(a,0)
			print(tmp_trans_dir)
                	pep_file = [x for x in glob.glob(os.path.join(tmp_trans_dir,"*.fasta"))][0]
                	self.assertTrue(os.path.isfile(pep_file))
                	self.assertTrue(os.stat(pep_file).st_size > 0)

                	orf_label_file = predict(model_dir, tmp_trans_dir, tmp_pred_dir, args)
                	self.assertTrue(os.path.isfile(orf_label_file))
                	self.assertTrue(os.stat(orf_label_file).st_size > 0)

                	resolve_label(tmp_resolve_dir, pep_file, orf_label_file, args)
                	read_label_file = os.path.join(tmp_resolve_dir, "reads.label")
                	self.assertTrue(os.path.isfile(read_label_file))
                	self.assertTrue(os.stat(orf_label_file).st_size > 0)

		aa_dir = os.path.join(out_dir,'abundance')
		calcAbundance(resolve_dir, aa_dir, config.map_file, config.dummy_gs_file)
		expected_filelist = ['raw_counts.tsv','effective_counts.tsv']
		out = os.popen('ls '+aa_dir).read()
		aa_filelist = out.split()
		self.assertTrue(set(aa_filelist) == set(expected_filelist))
                for f in aa_filelist:
                        self.assertTrue(os.stat(os.path.join(aa_dir,f)).st_size > 0)
		shutil.rmtree(out_dir)

if __name__ == '__main__':
    	parser = argparse.ArgumentParser()
    	parser.add_argument('-l',dest='frag_length', default=16, type=int)
    	parser.add_argument('-c',dest='coverage', default=0.1, type=float)
    	parser.add_argument('-n',dest='ncpus', default=5, type=int)
	parser.add_argument('--precise', action='store_true')
	parser.add_argument('-k', dest='kmer_length', default=8, type=int)
	parser.add_argument('--num-batches', default=1, type=int)
	parser.add_argument('--num-passes', default=1, type=int)
	parser.add_argument('--num_hash', dest='num_hash', default=1, type=int)
	parser.add_argument('--rweight', default=4, type=int)
	parser.add_argument('--hweight', default=-1, type=int)
	parser.add_argument('--bits', default=31, type=int)
	parser.add_argument('--lambda1', default=0., type=float)
	parser.add_argument('--lambda2', default=0., type=float)
    	parser.add_argument('--cutoff', default=0., type=float)
    	parser.add_argument('unittest_args', nargs='*')
    	args = parser.parse_args()
    	#fgsp_loc = os.path.abspath(args.fgsp_loc)
    	#if fgsp_loc.endswith('FGS+'):
        #	fgsp_loc = os.path.dirname(fgsp_loc)
    	sys.argv[1:] = args.unittest_args
    	unittest.main()

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
from carnelian import predict


import ldpc
import sequtil

class TestCarnelian_Basic(unittest.TestCase):
    '''
    Test the functions in carnelian.py
    '''
    
    def initialize(self):
        logging.getLogger('carnelian_basic_tests').addHandler(logging.NullHandler())

    def test01_unique_lines(self):
        '''
        Test the unique_lines function on a label file
        '''
        num_unique_labels = unique_lines(config.small_label_file)
        self.assertEqual(num_unique_labels, config.small_label_unique_labels)

    def test02_safe_makedirs_existing_dir(self):
        val = safe_makedirs(config.data_dir)
        self.assertEqual(val,1)   
    
    def test03_safe_makedirs_new_dir(self):
        val = safe_makedirs(os.path.join(config.data_dir,'tmp'))
        self.assertEqual(val,0)
        print(os.path.join(config.data_dir,'tmp'))
        shutil.rmtree(os.path.join(config.data_dir,'tmp'))
        
    def test04_merge_dico_overlapping(self):
	of = open(config.old_dico_file1,'w')
	of.write('2.7.11.1\t1\n2.7.13.3\t2\n2.7.7.6\t3\n3.6.4.12\t4')
	of.close()
	of = open(config.new_dico_file1,'w')
        of.write('3.6.4.12\t1\n4.1.1.11\t2')
        of.close()
        merge_dico(config.old_dico_file1, config.new_dico_file1)
        lines = [line.strip() for line in open(config.new_dico_file1)].sort()
        ref_lines = [line.strip() for line in open(config.merged_dico_file)].sort()
        self.assertTrue(lines == ref_lines)
	os.remove(config.old_dico_file1)
	os.remove(config.new_dico_file1)
        
    def test05_merge_dico_nonoverlapping(self):
	of = open(config.old_dico_file2,'w')
        of.write('2.7.11.1\t1\n2.7.13.3\t2\n2.7.7.6\t3')
        of.close()
        of = open(config.new_dico_file2,'w')
        of.write('3.6.4.12\t1\n4.1.1.11\t2')
        of.close()
        merge_dico(config.old_dico_file2, config.new_dico_file2)
        lines = [line.strip() for line in open(config.new_dico_file2)].sort()
        ref_lines = [line.strip() for line in open(config.merged_dico_file)].sort()
        self.assertTrue(lines == ref_lines)
	os.remove(config.old_dico_file2)
        os.remove(config.new_dico_file2)
    
    def test06_extract_column_two(self):
        safe_makedirs(os.path.join(config.data_dir,'tmp'))
        outpath = os.path.join(config.data_dir,'tmp')
        outfile = os.path.join(outpath,'tmp.txt')
        extract_column_two(config.two_col_file, outfile)
        self.assertTrue(filecmp.cmp(outfile, config.one_col_file))
        shutil.rmtree(outpath)
        
    def test07_vw_class_to_label(self):
        safe_makedirs(os.path.join(config.data_dir,'tmp'))
        outpath = os.path.join(config.data_dir,'tmp')
        outfile = os.path.join(outpath,'tmp.label')
        vw_class_to_label(config.small_vw_file, config.merged_dico_file, outfile)
        self.assertTrue(filecmp.cmp(outfile, config.small_label_file))
        shutil.rmtree(outpath)
        
    def test08_vw_class_to_label2(self):
        safe_makedirs(os.path.join(config.data_dir,'tmp'))
        outpath = os.path.join(config.data_dir,'tmp')
        outfile = os.path.join(outpath,'tmp.label')
        vw_class_to_label2(config.small_vw_prob_file, config.merged_dico_file, outfile, config.cutoff)
        self.assertTrue(filecmp.cmp(outfile, config.small_label_prob_file))
        self.assertTrue(filecmp.cmp(os.path.join(outpath,'label-probabilities.tsv'), config.small_label_prob_tsv_file))
        shutil.rmtree(outpath)
        
    def test09_vw_exists(self):
        self.assertTrue(find_executable('vw') is not None)
    
    def test10_fgsp_exists(self):
        self.assertTrue(find_executable(os.path.join(fgsp_loc,'FragGeneScan')) is not None)

    def test11_split_fasta(self):
        safe_makedirs(os.path.join(config.data_dir,'tmp'))
        outpath = os.path.join(config.data_dir,'tmp')
        ncpus = 5
        sequtil.split_fasta2(config.small_fasta_file,outpath,ncpus,my_env)
        fasta_filelist = [os.path.join(outpath,os.path.basename(x)) for x in glob.glob(os.path.join(outpath,'*.fasta'))]
        ref_filelist = [os.path.join(config.splits_dir,os.path.basename(x)) for x in glob.glob(os.path.join(config.splits_dir,'*.fasta'))]
	ind = [int(os.path.basename(f).split('.')[0]) for f in fasta_filelist]
        sorted_ind = np.argsort(ind)
        flist = [fasta_filelist[i] for i in sorted_ind]
	ind = [int(os.path.basename(f).split('.')[0]) for f in ref_filelist]
        sorted_ind = np.argsort(ind)
        rlist = [ref_filelist[i] for i in sorted_ind]
        for i in xrange(len(flist)):
            print('comparing... '+flist[i]+'\t'+rlist[i])
            self.assertTrue(filecmp.cmp(flist[i],rlist[i]))
        shutil.rmtree(outpath)
        
    def test12_merge_fasta(self):
        safe_makedirs(os.path.join(config.data_dir,'tmp'))
        outpath = os.path.join(config.data_dir,'tmp')
        filelist = [x for x in glob.glob(os.path.join(config.splits_dir,'*.fasta'))]
        ind = [int(os.path.basename(f).split('.')[0]) for f in filelist]
        sorted_ind = np.argsort(ind)
        flist = [filelist[i] for i in sorted_ind]
        sequtil.merge_files2(flist,outpath,'merged.peptides.fasta','fasta',my_env)
        outfile = os.path.join(outpath,'merged.peptides.fasta')
        self.assertTrue(filecmp.cmp(outfile,config.small_fasta_file))
        shutil.rmtree(outpath)

    def test13_translateSeqs_prot(self):
        self.assertEqual(translateSeqs(config.ftest_dir, config.pep_dir, fgsp_loc, args),1)

    def test14_translateSeqs_nucl(self):
	v = translateSeqs(config.ttest_dir, config.pep_dir, fgsp_loc, args)
        self.assertEqual(v,0)
        pep_fasta = [x for x in glob.glob(os.path.join(config.pep_dir,'*.fasta'))][0]
	self.assertTrue(os.path.isfile(pep_fasta))
	self.assertTrue(os.stat(pep_fasta).st_size > 0)
	shutil.rmtree(config.pep_dir)

    def test15_drawfrag_exists(self):
        self.assertTrue(find_executable(os.path.join(util_path,'drawfrag')) is not None)
        
    def test16_frag(self):
	v=frag(config.ftest_dir, config.frag_dir, args)
        self.assertEqual(v,0)
	frag_fasta = [x for x in glob.glob(os.path.join(config.frag_dir,'*.fasta'))][0]
        with open(frag_fasta) as f:
		first_line = f.readline()
        	self.assertTrue('>' in first_line)
        shutil.rmtree(config.frag_dir)
        
    def test17_fasta2skm_exists(self):
        self.assertTrue(find_executable(os.path.join(util_path,'fasta2skm')) is not None)
    
    def test18_train(self):
	args.precise=False
	if os.path.exists(config.model_dir):
		shutil.rmtree(config.model_dir)
        expected_filelist = ['patterns.txt','vw-dico.txt','vw-model_final.model']
        v = train(config.ftest_dir, config.model_dir, args)
        self.assertEqual(v,0)
        out = os.popen('ls '+config.model_dir).read()
        model_filelist = out.split()
        self.assertTrue(set(model_filelist) == set(expected_filelist))
        with open(os.path.join(config.model_dir,'vw-dico.txt')) as f:
            first_line = f.readline()
            x = first_line.split('\t')
            self.assertEqual(len(x),2)

    def test19_train_precise(self):
	args.precise=True
	if os.path.exists(config.model_precise_dir):
                shutil.rmtree(config.model_precise_dir)
        expected_filelist = ['patterns.txt','vw-dico.txt','vw-model_final.model']
        v = train(config.ftest_dir, config.model_precise_dir, args)
        self.assertEqual(v,0)
        out = os.popen('ls '+config.model_precise_dir).read()
        model_filelist = out.split()
        self.assertTrue(set(model_filelist) == set(expected_filelist))
        with open(os.path.join(config.model_precise_dir,'vw-dico.txt')) as f:
            first_line = f.readline()
            x = first_line.split('\t')
            self.assertEqual(len(x),2)

    def test20_predict(self):
	args.precise=False
	print(args.ncpus)
        label_file = predict(config.model_dir, config.ftest_dir, config.predict_dir, args)
	print(label_file)
        self.assertTrue(not(label_file is None))

    def test21_predict_precise(self):
	args.precise = True
	print(args.ncpus)
	args.cutoff = 0.80
        label_file = predict(config.model_precise_dir, config.ftest_dir, config.predict_precise_dir, args)
	self.assertTrue(not(label_file is None))
	with open(label_file.strip()) as f:
            first_line = f.readline()
            x = first_line.split('\t')
            self.assertEqual(len(x),2)

    def test22_evaluate_preds(self):
        val = evaluatePreds(config.small_label_file,config.small_label_pred_file)
        self.assertEqual(val,0)

    def test22_evaluate_preds_precise(self):
        val = evaluatePredsPrecise(config.small_label_file,config.small_label_prob_pred_file)
        self.assertEqual(val,0)

    def test23_retrain(self):
        args.precise=False
        if os.path.exists(config.new_model_dir):
                shutil.rmtree(config.new_model_dir)
        expected_filelist = ['patterns.txt','vw-dico.txt','vw-model_final.model']
        v = retrain(config.model_dir, config.new_model_dir, config.new_examples_dir, args)
        self.assertEqual(v,0)
        out = os.popen('ls '+config.new_model_dir).read()
        model_filelist = out.split()
        self.assertTrue(set(model_filelist) == set(expected_filelist))
        with open(os.path.join(config.new_model_dir,'vw-dico.txt')) as f:
            first_line = f.readline()
            x = first_line.split('\t')
            self.assertEqual(len(x),2)

    def test24_retrain_precise(self):
        args.precise=True
        if os.path.exists(config.new_model_precise_dir):
                shutil.rmtree(config.new_model_precise_dir)
        expected_filelist = ['patterns.txt','vw-dico.txt','vw-model_final.model']
        v = retrain(config.model_precise_dir, config.new_model_precise_dir, config.new_examples_dir, args)
        self.assertEqual(v,0)
        out = os.popen('ls '+config.new_model_precise_dir).read()
        model_filelist = out.split()
        self.assertTrue(set(model_filelist) == set(expected_filelist))
        with open(os.path.join(config.new_model_precise_dir,'vw-dico.txt')) as f:
            first_line = f.readline()
            x = first_line.split('\t')
            self.assertEqual(len(x),2)

    def test25_majority_voting_winner(self):
        L1 = ['a','a','a','a','a','a']
        L2 = ['a','a','a','a','a','c']
        L3 = ['a','a','a','a','b','b']
        L4 = ['a','a','a','b','b','d']
        L5 = ['a','a','b','c','d','e']
        self.assertEqual(majority_voting(L1),'a')
        self.assertEqual(majority_voting(L2),'a')
        self.assertEqual(majority_voting(L3),'a')
        self.assertEqual(majority_voting(L4),'a')
        self.assertEqual(majority_voting(L5),'a')

    def test26_majority_voting_nowinner(self):
        L1 = ['a','a','a','b','b','b']
        L2 = ['a','a','b','b','c','c']
        L3 = ['a','b','c','d','e','f']
        L4 = ['a','a','b','b','c','d']
        self.assertEqual(majority_voting(L1),'N.N.N.N')
        self.assertEqual(majority_voting(L2),'N.N.N.N')
        self.assertEqual(majority_voting(L3),'N.N.N.N')
        self.assertEqual(majority_voting(L4),'N.N.N.N')

    def test27_resolve_label(self):
	args.precise = False
        outpath = os.path.join(config.data_dir,'tmp')
        readpath = resolve_label(outpath,config.orf_seq,config.orf_label, args)
	print(readpath)
	ref_lines = [line.strip() for line in open(config.reads_label)]
	out_lines = [line.strip() for line in open(readpath)]
	self.assertEqual(ref_lines.sort(),out_lines.sort())
        shutil.rmtree(outpath)

    def test28_resolve_label_precise(self):
	args.precise = True
        outpath = os.path.join(config.data_dir,'tmp')
        readpath = resolve_label(outpath,config.orf_seq,config.orf_precise_label, args)
        ref_lines = [line.strip() for line in open(config.reads_precise_label)]
        out_lines = [line.strip() for line in open(readpath)]
	self.assertEqual(ref_lines.sort(),out_lines.sort())
        shutil.rmtree(outpath)

#    def test_resolve_label_wrt_ref(self):
#	 args.precise = False
#        outpath = os.path.join(config.data_dir,'tmp')
#        readpath = resolve_label_wrt_ref(outpath,config.orf_seq,config.orf_label, config.read_seq, args)
#        self.assertTrue(filecmp.cmp(config.reads_ref_label,readpath))
#        shutil.rmtree(outpath)

#    def test_resolve_label_wrt_ref_precise(self):
#	args.precise = True
#        outpath = os.path.join(config.data_dir,'tmp')
#        readpath = resolve_label_wrt_ref(outpath,config.orf_seq,config.orf_precise_label, config.read_seq, args)
#        self.assertTrue(filecmp.cmp(config.reads_ref_label_precise,readpath))
#        shutil.rmtree(outpath)

    def test29_R_installation(self):
        self.assertTrue(find_executable('Rscript') is not None)

    def test30_create_abundance_matrix(self):
        safe_makedirs(os.path.join(config.data_dir,'tmp')) 
        outpath = os.path.join(config.data_dir,'tmp')
        createAbundanceMatrix(config.dummy_aa_in_dir, outpath, config.dummy_gs_file)
        outfile = os.path.join(outpath,'raw_counts.tsv')
        reffile = os.path.join(config.dummy_aa_out_dir,'raw_counts.tsv')
        self.assertTrue(filecmp.cmp(outfile,reffile))
        shutil.rmtree(outpath)

    def test31_abundance_estimation_rscript(self):
        safe_makedirs(os.path.join(config.data_dir,'tmp')) 
        outpath = os.path.join(config.data_dir,'tmp')
        rscript_loc = '../scripts/abundance_estimation.R'
        counts_file = os.path.join(config.dummy_aa_out_dir,'raw_counts.tsv')
    	command =  counts_file + ' ' + config.dummy_gs_file + ' ' + config.dummy_mapping_file + ' ' + outpath
    	os.system('Rscript '+rscript_loc+ ' ' + command)
        outfile = os.path.join(outpath,'effective_counts.tsv')
        reffile = os.path.join(config.dummy_aa_out_dir,'effective_counts.tsv')
        self.assertTrue(filecmp.cmp(outfile,reffile))
        shutil.rmtree(outpath)


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
    	sys.argv[1:] = args.unittest_args
    	unittest.main()

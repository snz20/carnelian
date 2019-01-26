# -*- coding: utf-8 -*-
'''
splits a large fasta file into n smaller ones
'''
import math
import os
import argparse
from Bio import SeqIO
import re
import random

def batch_iterator(iterator, batch_size):
    """Returns lists of length batch_size.

    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.

    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
    entry = True  # Make sure we loop once
    while entry:
        batch = []
        while len(batch) < batch_size:
            try:
                entry = iterator.next()
            except StopIteration:
                entry = None
            if entry is None:
                # End of file
                break
            batch.append(entry)
        if batch:
            yield batch

def new_names(oldname, n, suffix):
    p = oldname.rfind(suffix)
    if p != -1:
        pattern = oldname[:p] + "%s." + oldname[p:]
    else:
        pattern = oldname + ".%s"
    if n==1:
        return [oldname]
    else:
        width = len(str(n))
        names = [pattern % str(i).rjust(width, '0') for i in range(n)]
        return names


def safe_makedirs(directory):
        if not os.path.exists(directory):
                os.makedirs(directory)
        else:
                print("Directory already exists!!")
                return 1
        return 0

def check_if_nucl(strg, search=re.compile(r'[^ATCGN\*-]').search):
	return not bool(search(strg))

def split_fasta(filename, dest_dir, num_splits):
	record_iter = SeqIO.parse(open(filename),"fasta")
	num_recs = sum(1 for _ in record_iter)
	batch_size = int(math.ceil(num_recs*1.0/num_splits))
	record_iter = SeqIO.parse(open(filename),"fasta")
	num_counts = []
	width = len(str(num_splits))
	names = new_names(os.path.basename(filename),num_splits,"fasta")
	for i, batch in enumerate(batch_iterator(record_iter, batch_size)):
		dest_path = os.path.join(os.path.abspath(dest_dir), 'split'+str(i).rjust(width, '0'))
        	safe_makedirs(dest_path)
    		filename = os.path.join(dest_path,names[i])
    		with open(filename, "w") as handle:
        		count = SeqIO.write(batch, handle, "fasta")
    			print("Wrote %i records to %s" % (count, filename))
			num_counts.append(count)
	return num_counts


def split_labels(filename, dest_dir, splits):
	lines = [line.strip() for line in open(filename)]
	num_splits = len(splits)
	#dest_path = os.path.abspath(dest_dir)
	old_st = 0
	old_chunk = 0
	width = len(str(num_splits))
	names = new_names(os.path.basename(filename),num_splits,"label")
	for i in xrange(num_splits):
		dest_path = os.path.join(os.path.abspath(dest_dir), 'split'+str(i).rjust(width, '0'))
		outfile = os.path.join(dest_path,names[i])
		new_st = old_st + old_chunk
		with open(outfile, 'w') as handle:
			handle.writelines('\n'.join(lines[new_st:new_st+splits[i]])+'\n')
			handle.close()
		old_st = new_st
		old_chunk = splits[i]

def merge_fasta(filelist, out_dir,prefix):
	print("merged fasta: " + os.path.join(out_dir,prefix)+'.peptides.fasta')
	out_path = os.path.join(out_dir,prefix)
	with open(out_path+'.peptides.fasta', 'w') as w_file:
		for filen in filelist:
			with open(filen, 'rU') as o_file:
            			seq_records = SeqIO.parse(o_file, 'fasta')
            			SeqIO.write(seq_records, w_file, 'fasta')

def merge_files(filelist, out_dir, infile_name, ext):
	print(infile_name)
	prefix = infile_name.rsplit('.',1)[0];
	out_path = os.path.join(out_dir,prefix)
	with open(out_path+ext, 'w') as w_file:
		for fl in filelist:
			with open(fl,'rU') as o_file:
				for line in o_file:
					w_file.write(line)
	return out_path+ext

def read_to_orf_label(label_file, out_dir, prefix):
	print(prefix)
	out_path = os.path.join(out_dir,prefix+'.label')
	print(out_path)
	with open(out_path,'w') as w_file:
		with open(label_file, 'rU') as in_file:
			for line in in_file:
				for i in range(6):
					w_file.write(line)

def estimate_avg_frag_len(fasta_file):
	print("Estimating read length for fasta: "+ fasta_file)
	recs = list(SeqIO.parse(fasta_file,'fasta'))
	num_sample = int(len(recs)*0.01)
	samp_recs = random.sample(recs,num_sample)
	sum_len = 0
	for seq_record in samp_recs:
		sum_len += len(seq_record.seq)
	avg_len = sum_len*1.0/num_sample
	return avg_len

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('-f', nargs=1)
	parser.add_argument('-o', nargs=1)
	parser.add_argument('-n', nargs=1)
	args = parser.parse_args()
	#print args.k, args.t, args.m
	filename = os.path.abspath(args.f[0])
	dest_dir = args.o[0]
	num_splits = float(args.n[0])
	print filename
	num_counts = split_fasta(filename, dest_dir, num_splits)
	print(num_counts)

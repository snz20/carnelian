# -*- coding: utf-8 -*-
'''
splits a large fasta file into n smaller ones
'''
import math
import os
from Bio import SeqIO

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




def split_fasta(filename, dest_dir, num_splits):
	record_iter = SeqIO.parse(open(filename),"fasta")
	num_recs = sum(1 for _ in record_iter)
	batch_size = int(math.ceil(num_recs*1.0/num_splits))
	record_iter = SeqIO.parse(open(filename),"fasta")
	for i, batch in enumerate(batch_iterator(record_iter, batch_size)):
    		filename = os.path.join(dest_dir,"%i.fasta" % (i + 1))
    		with open(filename, "w") as handle:
        		count = SeqIO.write(batch, handle, "fasta")
    		print("Wrote %i records to %s" % (count, filename))


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
	prefix = infile_name.split('.')[0];
	out_path = os.path.join(out_dir,prefix)
	with open(out_path+ext, 'w') as w_file:
		for fl in filelist:
			with open(fl,'rU') as o_file:
				for line in o_file:
					w_file.write(line)
 

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('-f', nargs=1)
	parser.add_argument('-o', nargs=1)
	parser.add_argument('-n', nargs=1)
	args = parser.parse_args()
	#print args.k, args.t, args.m
	filename = args.f
	dest_dir = args.o
	num_splits = args.n
	split_fasta(filename, dest_dir, num_splits)

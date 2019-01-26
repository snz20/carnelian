#!/usr/bin/python
#merge_pairs.py
# based on Jeremy Leipzig's mergePairs.py
#merge overlapping paired end reads. In the case of no overlap, link them by a sequence of Ns which has the length of inner gap between two reads.

from optparse import OptionParser
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.Seq import Seq
from Bio import pairwise2
import itertools
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import SeqIO
import numpy as np
from StringIO import StringIO


usage = "usage: %prog [options] FILE1 FILE2"
parser = OptionParser(usage=usage)

parser.add_option("-i", "--insert", action="store", type="int", dest="insert",
                  help="insert length of fragment from end to end (|s1---s2|)",metavar="insert length (default:%default)",default=125)
parser.add_option("-m", "--minimum", action="store", type="int", dest="minimum",
                  help="absolute minimum acceptable bp overlap for pairs to be merged (default: 5)",default=5,metavar="minimum")
parser.add_option("-g", "--merged", action="store", type="string", dest="merged",
                  help="output file for successfully merged pairs (default merged.fq)",metavar="mergedFileName.fq",default="merged.fq")
parser.add_option("-d", "--identity", action="store", type="int", dest="identity",
                  help="minimum percent sequence identity within overlap (default: 90)",default=90,metavar="sequence identity")

(options, args) = parser.parse_args()
if len(args) != 2:
	parser.error("incorrect number of arguments")
else:
	fwd=args[0]
	rev=args[1]

insert = int(options.insert)
merged_handle = open(options.merged, "w")

def merge_pairs(seq1, id1, q1, seq2, id2, q2):
		global merged
		exact_pos = seq1.find(seq2[0:options.minimum])
		if exact_pos >= 0:
				seq2_region = seq2[0:len(seq2)-(len(seq1)-exact_pos)]
				#this matrix is necessary otherwise N-N alignments are considered matches
				jerm={('A', 'A'): 1, ('A', 'C'): -1, ('A', 'T'): -1, ('A', 'G'): -1,
				      ('G', 'G'): 1, ('G', 'C'): -1, ('G', 'T'): -1,
				      ('C', 'C'): 1, ('C', 'T'): -1, 
				      ('T', 'T'): 1,
				      ('N','N'): -1,('N','A'): -1,('N','C'): -1,('N','G'): -1,('N','T'): -1}

                                #+1/-1 scoring is somehow necessary, a 1/0 scoring tends to produce awful end-alignments
				alignments = pairwise2.align.globalds(seq1,seq2,jerm,-1,-1,penalize_end_gaps=False,one_alignment_only=True)
				if len(alignments) < 1:
					printUnmerged(id1,seq1,q1,id2,seq2,q2)
				for seq_a, seq_b, score, start, end in alignments:
					overlap=len(seq1)+len(seq2)-(end-start)
					endgaps=(end-start)-overlap
					if score>=overlap-overlap*2*(1-(options.identity/100.0)):
						apos=0
						bpos=0
						seq=''
						qual=''
						for (a1,a2) in itertools.izip(seq_a,seq_b):
							if (a1=='-' or a1 == 'N') and (a2=='-' or a2 == 'N'):
									seq=''
									qual=''
									printUnmerged(id1,seq1,q1,id2,seq2,q2, insert)
									break
							else:
								if (a1=='-'):
									seq+=a2
									qual+=q2[bpos]
									bpos+=1
								elif (a2=='-'):
									seq+=a1
									qual+=q1[apos]
									apos+=1
								elif (a1=='N'):
									seq+=a2
									qual+=q1[bpos]
									apos+=1
									bpos+=1
								elif (a2=='N'):
									seq+=a1
									qual+=q1[apos]
									apos+=1
									bpos+=1
								elif a1!='-' and a2!='-' and a1!='N' and a2!='N':
									if q1[apos]>q2[bpos]:
										seq+=a1
										qual+=q1[apos]
									else:
										seq+=a2
										qual+=q2[bpos]
									apos+=1
									bpos+=1
						if(len(seq)>0):
							printMerged(id1,seq,qual)
							assert(apos==len(seq1))
							assert(bpos==len(seq2))
							merged+=1
					else:
						printUnmerged(id1,seq1,q1,id2,seq2,q2, insert)
		else:
			printUnmerged(id1,seq1,q1,id2,seq2,q2, insert)
			

def printUnmerged(id1,seq1,q1,id2,seq2,q2,insert):
		id = id1
		mid = ''.join(['N' for _ in xrange(int(insert))])
		seq = seq1+mid+seq2
		fastq_string = ">%s\n%s\n" % (id1, seq)
		merged_handle.write(fastq_string)

def printPreserved(id1,seq1,q1,id2,seq2,q2):
		fastq_string = "@%s\n%s\n+\n%s\n" % (id1, seq1, q1)
		record = SeqIO.read(StringIO(fastq_string), formats[options.intypenum])
		preserve_handle.write(str(record.format(formats[options.outtypenum])))
		fastq_string = "@%s\n%s\n+\n%s\n" % (id2, seq2, q2)
		record = SeqIO.read(StringIO(fastq_string), formats[options.intypenum])
		preserve_handle.write(str(record.format(formats[options.outtypenum])))
		
def printMerged(id1,seq,qual):
		fastq_string = ">%s\n%s\n" % (id1, seq)
		merged_handle.write(fastq_string)


merged=0
count=0
formats={1:'fastq-sanger',2:'fastq-illumina',3:'fastq-solexa',4:'fasta'}
f_iter = FastqGeneralIterator(open(fwd,"rU"))
r_iter = FastqGeneralIterator(open(rev,"rU"))
for (f_id, f_seqstr, f_q), (r_id, r_seqstr, r_q) in itertools.izip(f_iter,r_iter):
	f_seq = Seq(f_seqstr)
	r_seq = Seq(r_seqstr)
	count += 2
	print f_seq
	print r_seq.reverse_complement()
	merge_pairs(str(f_seq).upper(),f_id,f_q,str(r_seq.reverse_complement()).upper(),r_id,r_q[::-1])


print "%i records processed. %s merged." % (count, merged)


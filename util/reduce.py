# -*- coding: utf-8 -*-
'''
utility functions for reducing a fasta file using simpler amino acid alphabet
author: Sumaiya Nazeen
usage:
python reduce.py infile outfile alphabet
alph can be any of the following set = {Murphy15, Murphy10, Murphy8, Murphy4, PC5, HPModel}
'''
from Bio.Seq import Seq
from Bio import Alphabet
from Bio.Alphabet import Reduced
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys
import random

def reduceSeq(infile, outfile, alph):
        rec = []
        with open(infile,'rU') as input_handle:
                for record in SeqIO.parse(input_handle,"fasta"):
                        rec.append(record)
        new_p = []
        for r in rec:
                n_p = Seq('',Alphabet.ProteinAlphabet())
                if alph == 'Murphy10':
                        n_p = Seq('',Alphabet.Reduced.Murphy10())
                elif alph == 'Murphy15':
                        n_p = Seq('',Alphabet.Reduced.Murphy15())
                elif alph == 'Murphy8':
                        n_p = Seq('',Alphabet.Reduced.Murphy8())
                elif alph == 'Murphy4':
                        n_p = Seq('',Alphabet.Reduced.Murphy4())
                elif alph == 'PC5':
                        n_p = Seq('',Alphabet.Reduced.PC5())
                elif alph == 'HPModel':
                        n_p = Seq('',Alphabet.Reduced.HPModel())
                for aa in r:
                        if aa != '*' and aa != '-' and aa != 'U':
                                if aa == 'X':
                                        aa = random.sample(set('ACDEFGHIKLMNPQRSTVWY'), 1)[0]
                                n_p += Alphabet.Reduced.murphy_10_tab[aa]
                        else:
                                n_p += aa
                x = SeqRecord(n_p)
                x.id = r.id
                x.description = r.description
                new_p.append(x)
        SeqIO.write(new_p, outfile, "fasta")

def main():
        if len(sys.argv) < 3:
                print "Not enough arguments!"
                exit(1)
        infile  = sys.argv[1]
        outfile = sys.argv[2]
        alph = sys.argv[3]
        reduceSeq(infile, outfile,alph)

if __name__ == "__main__":
        main()

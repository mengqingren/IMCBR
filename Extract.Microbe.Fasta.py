#!/usr/bin/python
import sys
import argparse
from Bio import SeqIO
import os
from collections import defaultdict
from collections import Counter

parser = argparse.ArgumentParser(description='This manuscript for extract microbe.fa reference for PathSeq')
parser.add_argument("--Kraken2DB", type=str, default="./")
parser.add_argument("--TaxidFile", type=str, default="Project.Kraken2.Species.ID.txt")
parser.add_argument("--Output", type=str, default="microbe.fa")
parser.add_argument("--OutputRef", type=str, default="Ref_Microbe.txt")
parser.add_argument('-q', '--quiet', action='store_true', help='Print quiet')
parser.add_argument('-v', '--verbose', action='store_true', help='Print verbose')
args = parser.parse_args()

Taxid_Species = defaultdict()
with open(args.TaxidFile,"r") as f:
	for line in f:
		line=line.strip().split("\t")
		Taxid_Species[line[1]] = line[0].replace(" ","_")+"|"+line[1]
#print(Taxid_Species.keys())
os.popen("touch "+args.Output)
os.popen("touch "+args.OutputRef)
data = open(args.OutputRef,"w")
with open(args.Output,"w") as file:
	for seq_record in SeqIO.parse(args.Kraken2DB+"/library/bacteria/library.fna", "fasta"):
            taxid=(seq_record.id.strip().split("|"))[1]
            ReferrenceID=seq_record.id.strip().split("|")[2]
            if taxid in Taxid_Species.keys():
                file.write(">"+ReferrenceID+" "+Taxid_Species[taxid]+"\t"+taxid+"\n")
                data.write(ReferrenceID+"\t"+Taxid_Species[taxid]+"\t"+taxid+"\n")
                file.write(str(seq_record.seq)+"\n")

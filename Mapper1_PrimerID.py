#!/usr/bin/env python
from Bio import SeqIO
import sys
from Bio.Seq import Seq
from distance import hamming

def readref(reffile):
  refdict = {}
  infile = open(reffile)
  for record in SeqIO.parse(infile,'fasta'):
    refdict[str(record.id)] = str(record.seq)
  infile.close()
  return refdict

def mapping(seq,primers):
  offsets = []; direction = 'NA'
  #read1: barcode_F:barcode:barcode_R / read2: primerID_F: primerID : primerID_R : rc(barcode_R) : rc(barcode) : rc(barcode_F)
  primer = 'Barcode_F'
  priseq = primers[primer]
  for offset_pri in range(0,3):
    for offset_seq in range(35,55):
      pri_window = priseq[:len(priseq)-offset_pri]
      seq_window = seq[offset_seq:offset_seq+len(pri_window)]
      if len(pri_window) != len(seq_window) : continue
      if hamming(pri_window,seq_window) < 3:
        direction = 'F'
        offsets.append(offset_seq+len(pri_window))
        break
    if direction == 'F':
      break
  if direction == 'F':
    primer = 'Barcode_R'
    priseq = primers[primer]
    for offset_pri in range(0,3):
      for offset_seq in range(offsets[0]+10,offsets[0]+30):
        pri_window = priseq[:len(priseq)-offset_pri]
        seq_window = seq[offset_seq:offset_seq+len(pri_window)]
        if len(pri_window) != len(seq_window) : continue
        if hamming(pri_window,seq_window) < 3:
          offsets.append(offset_seq)
          return direction, offsets
  if len(offsets) == 1:
    return direction, offsets
  if len(offsets) == 0:
    primer = 'Primer_F'
    priseq = primers[primer]
    for offset_pri in range(0,3):
      for offset_seq in range(7,17):
        pri_window = priseq[:len(priseq)-offset_pri]
        seq_window = seq[offset_seq:offset_seq+len(pri_window)]
        if len(pri_window) != len(seq_window) : continue
        if hamming(pri_window,seq_window) < 3:
          direction = 'R'
          offsets.append(offset_seq+len(pri_window))
          break
      if direction == 'R':
        break
  if direction == 'R':
    primer = 'Primer_R'
    priseq = primers[primer]
    for offset_pri in range(0,3):
      for offset_seq in range(offsets[0]+7,offsets[0]+17):
        pri_window = priseq[:len(priseq)-offset_pri]
        seq_window = seq[offset_seq:offset_seq+len(pri_window)]
        if len(pri_window) != len(seq_window) : continue
        if hamming(pri_window,seq_window) < 3:
          offsets.append(offset_seq)
          break
      if len(offsets) == 2:
        break
    if len(offsets) == 1:
      return direction, offsets
    primer = 'Barcode_R'
    priseq = str(Seq(primers[primer]).reverse_complement())
    for offset_pri in range(0,3):
      for offset_seq in range(offsets[1]+30,offsets[1]+50):
        pri_window = priseq[:len(priseq)-offset_pri]
        seq_window = seq[offset_seq:offset_seq+len(pri_window)]
        if len(pri_window) != len(seq_window) : continue
        if hamming(pri_window,seq_window) < 3:
          offsets.append(offset_seq+len(pri_window))
          break
      if len(offsets) == 3:
        break
    if len(offsets) == 2:
      return direction, offsets
    primer = 'Barcode_F'
    priseq = str(Seq(primers[primer]).reverse_complement())
    for offset_pri in range(0,3):
      for offset_seq in range(offsets[2]+16,offsets[2]+26):
        pri_window = priseq[:len(priseq)-offset_pri]
        seq_window = seq[offset_seq:offset_seq+len(pri_window)]
        if len(pri_window) != len(seq_window) : continue
        if hamming(pri_window,seq_window) < 3:
          offsets.append(offset_seq)
          break
      if len(offsets) == 4:
        break
    if len(offsets) == 3:
      return direction, offsets
  return direction, offsets

def consens2(seq1,seq2,qual1,qual2):
  if len(seq1) != len(seq2) or len(qual1) != len(qual2): return 'NA'
  seq = ''
  for i in range(len(seq1)):
    if seq1[i] == seq2[i]:
      seq += seq1[i]
    elif qual1[i] > qual2[i] and qual2[i] > 15:
      seq += seq1[i]
    elif qual2[i] > qual1[i] and qual1[i] > 15:
      seq += seq2[i]
    else:
      seq = 'NA'
      break
  return seq

def main():
  workpath = '/u/scratch/t/tianhao/HiSeq080119/'
  infile1 = open(workpath+'split/'+sys.argv[1])
  infile2 = open(workpath+'split/'+sys.argv[1].replace('_R1_','_R2_'))
  outfile = open(workpath+'barcode/'+sys.argv[1].replace('.fastq','.txt'),'w')  #sample_barcode : virus_barcode : primerID
  primers = readref('../Fasta/primer.fasta')
  inhandle1 = SeqIO.parse(infile1,'fastq')
  inhandle2 = SeqIO.parse(infile2,'fastq')
  counts = [0,0,0,0,0]; readcount = 0; goodcount = 0
  #error counts: 0. unpaired tag; 1.incomplete mapping; 2.same strand mapping; 3.unmapped; 4.unpaired barcode
  archfile1 = {}
  archfile2 = {}
  tagdict = {}
  tagfile = open('../Fasta/tag_Jocelyn_080119.txt')
  for line in tagfile:
    line = line.rstrip().rsplit('\t')
    tagdict[line[0]] = line[1]
    archfile1[line[0]] = open(workpath+'archive/'+sys.argv[1].replace('.fastq','')+'_'+line[1]+'.fastq','w')
    archfile2[line[0]] = open(workpath+'archive/'+sys.argv[1].replace('.fastq','').replace('_R1_','_R2_')+'_'+line[1]+'.fastq','w')
  tagfile.close()
  for record1 in inhandle1:
    readcount += 1
    if readcount % 100000 == 0: print('readline:'+str(readcount)+' Good count:'+str(goodcount)+' Error count:',counts)
    record2 = inhandle2.next()
    seq1 = str(record1.seq)
    seq2 = str(record2.seq)
    qual1 = record1.letter_annotations["phred_quality"]
    qual2 = record2.letter_annotations["phred_quality"]
    tag = consens2(seq1[:6],seq2[:6],qual1[:6],qual2[:6])
    if tag == 'NA' or tag not in tagdict: counts[0] += 1; continue
    archfile1[tag].write(record1[7:].format('fastq'))
    archfile2[tag].write(record2[7:].format('fastq'))
    direction1, offsets1 = mapping(seq1,primers)
    direction2, offsets2 = mapping(seq2,primers)
    #print(direction1,offsets1)
    #print(direction2,offsets2)
    if direction1 == 'NA' or direction2 == 'NA': counts[3] += 1; continue
    if direction1 == 'F' and len(offsets1) != 2: counts[1] += 1; continue
    if direction1 == 'R' and len(offsets1) != 4: counts[1] += 1; continue
    if direction2 == 'F' and len(offsets2) != 2: counts[1] += 1; continue
    if direction2 == 'R' and len(offsets2) != 4: counts[1] += 1; continue
    if direction1 == direction2: counts[2] += 1; continue
    if direction1 == 'R': 
      UMI  = seq1[offsets1[0]:offsets1[1]]
      bcrc = seq1[offsets1[2]:offsets1[3]]
      bc_F = seq2[offsets2[0]:offsets2[1]]
      bcrc_qual = qual1[offsets1[2]:offsets1[3]]
      qual_F    = qual2[offsets2[0]:offsets2[1]]
    if direction2 == 'R':
      UMI  = seq2[offsets2[0]:offsets2[1]]
      bcrc = seq2[offsets2[2]:offsets2[3]]
      bc_F = seq1[offsets1[0]:offsets1[1]]
      bcrc_qual = qual2[offsets2[2]:offsets2[3]]
      qual_F    = qual1[offsets1[0]:offsets1[1]]
    bc_R = str(Seq(bcrc).reverse_complement())
    qual_R = bcrc_qual[::-1]
    bc = consens2(bc_F,bc_R,qual_F,qual_R)
    if bc == 'NA': counts[4] += 1; continue
    goodcount += 1
    outfile.write(tag+'\t'+bc+'\t'+UMI+'\n')
  for tag in tagdict:
    archfile1[tag].close()
    archfile2[tag].close()
  infile1.close()
  infile2.close()
  outfile.close()
     
if __name__ == '__main__':
  main()

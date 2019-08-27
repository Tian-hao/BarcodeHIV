#!/usr/bin/env python
import math
import glob
import itertools

def hamming(str1, str2):
  return sum(itertools.imap(str.__ne__, str1, str2))

def main():
  bcfiles = sorted(glob.glob('analysis/HC5_cumfreq90/*.txt')) 
  bcdict = {}
  for bcfile in bcfiles:
    if '_R1' in bcfile or '_R2' in bcfile: continue
    inhandle = open(bcfile)
    sample = bcfile.rsplit('/')[-1].rsplit('.txt')[0]
    sample = '_'.join(sample.rsplit('_')[2:7])
    sample = sample.replace('#','')
    bcdict[sample] = {}
    sumcount = 0
    for line in inhandle:
      line = line.rstrip().rsplit('\t')
      barcode = line[0]
      readcount = int(line[1])
      bcdict[sample][barcode] = readcount
      sumcount += readcount
    for barcode in bcdict[sample]:
      bcdict[sample][barcode] /= float(sumcount)
    inhandle.close()

  outfile_bc = open('analysis/contamination/HC5_cumfreq90/PID_overlap_bc.txt','w')
  outfile_fr = open('analysis/contamination/HC5_cumfreq90/PID_overlap_freq.txt','w')
  outfile_bc.write('samples')
  outfile_fr.write('samples')
  samplelist = []
  for sample in bcdict:
    samplelist.append(sample)
  print(samplelist)
  samplelist = sorted(samplelist, key=lambda x: x.rsplit('_')[0])
  samplelist = sorted(samplelist, key=lambda x: x.rsplit('_')[2])
  samplelist = sorted(samplelist, key=lambda x: x.rsplit('_')[1])
  for sample in samplelist:
    outfile_bc.write('\t'+sample)
    outfile_fr.write('\t'+sample)
  outfile_bc.write('\n')
  outfile_fr.write('\n')
  for i in range(len(samplelist)):
    sample1 = samplelist[i]
    outfile_bc.write(sample1)
    outfile_fr.write(sample1)
    for j in range(len(samplelist)):
      sample2 = samplelist[j]
      ovdict = overlap(bcdict[sample1],bcdict[sample2],0)
      sumfreq = 0
      for bc in ovdict:
        sumfreq += ovdict[bc]
      outfile_bc.write('\t'+str(len(ovdict.keys())))
      outfile_fr.write('\t'+str(sumfreq))
    outfile_bc.write('\n')
    outfile_fr.write('\n')
  outfile_bc.close()
  outfile_fr.close()

def overlap(bcdict1,bcdict2,distance):
  ovdict = {}
  for bc1 in bcdict1:
    for bc2 in bcdict2:
      if hamming(bc1,bc2) <= distance:
        ovdict[bc1] = math.sqrt(bcdict1[bc1]*bcdict2[bc2])
  return ovdict
  

if __name__ == '__main__':
  main()

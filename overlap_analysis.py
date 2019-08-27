#!/usr/bin/env python
import glob
from math import sqrt

def main():
  plasmasamples = sorted(glob.glob('analysis/HC5_Clustered_freq0001/*plasma*.txt'))
  spleensamples = sorted(glob.glob('analysis/HC5_Clustered_freq0001/*spleen*.txt')) 
  compare_overlap(plasmasamples,spleensamples,'HC5_Clustered_freq0001','PID')
  
  plasmasamples = sorted(glob.glob('analysis/HC5_Clustered_cumfreq90/*plasma*.txt'))
  spleensamples = sorted(glob.glob('analysis/HC5_Clustered_cumfreq90/*spleen*.txt')) 
  compare_overlap(plasmasamples,spleensamples,'HC5_Clustered_cumfreq90','PID')

  plasmasamples = sorted(glob.glob('analysis/HC5_Clustered/*plasma*.txt'))
  spleensamples = sorted(glob.glob('analysis/HC5_Clustered/*spleen*.txt')) 
  compare_overlap(plasmasamples,spleensamples,'HC5_Clustered','PID')

  plasmasamples = sorted(glob.glob('analysis/HC5_freq0001/*plasma*.txt'))
  spleensamples = sorted(glob.glob('analysis/HC5_freq0001/*spleen*.txt')) 
  compare_overlap(plasmasamples,spleensamples,'HC5_freq0001','PID')
  
  plasmasamples = sorted(glob.glob('analysis/HC5_cumfreq90/*plasma*.txt'))
  spleensamples = sorted(glob.glob('analysis/HC5_cumfreq90/*spleen*.txt')) 
  compare_overlap(plasmasamples,spleensamples,'HC5_cumfreq90','PID')

  plasmasamples = sorted(glob.glob('analysis/HC5/LRA/*plasma*.txt')) + sorted(glob.glob('analysis/HC5/CTL/*plasma*.txt'))
  spleensamples = sorted(glob.glob('analysis/HC5/LRA/*spleen*.txt')) + sorted(glob.glob('analysis/HC5/CTL/*spleen*.txt'))
  compare_overlap(plasmasamples,spleensamples,'HC5','PID')

def compare_overlap(plasmasamples,spleensamples,annotation,annot2):
  #plasmasamples = sorted(glob.glob('analysis/HC5_Clustered_freq0001/*plasma*.txt'))
  #spleensamples = sorted(glob.glob('analysis/HC5_Clustered_freq0001/*spleen*.txt')) 
  samplelist = []
  spdict = {}
  pldict = {}
  for plasmasample in plasmasamples:
    if '_R1' in plasmasample or '_R2' in plasmasample: continue
    sampleNo = plasmasample.rsplit('plasma_')[-1].rsplit('.txt')[0]
    _paired = 0
    for ss in spleensamples:
      if sampleNo in ss: 
        spleensample = ss
        samplelist.append(sampleNo)
        _paired = 1
    if _paired == 0: continue
    sphandle = open(spleensample)
    plhandle = open(plasmasample)
    spdict[sampleNo] = {}
    pldict[sampleNo] = {}
    for line in sphandle:
      line = line.rstrip().rsplit('\t')
      spdict[sampleNo][line[0]] = float(line[1])
    sphandle.close()
    for line in plhandle:
      line = line.rstrip().rsplit('\t')
      pldict[sampleNo][line[0]] = float(line[1])
    sphandle.close()
    plhandle.close()
  #print barcodes sequences for each sample
  outfiles = {}
  totalcount = {}
  ovdict = {}
  simdict = {}
  for sampleNo in samplelist:
    outfiles[sampleNo] = open('analysis/overlap/'+annotation+'/plasma_spleen_'+annot2+'_barcodes_'+sampleNo+'.txt','w')
    outfiles[sampleNo].write('barcode\tplasma_PIDcount\tspleen_PIDcount\ttotal_PIDcount\n')
    totalcount[sampleNo] = {}
    ovdict[sampleNo] = {}
    for bc in spdict[sampleNo]:
      if bc in pldict[sampleNo]:
        totalcount[sampleNo][bc] = spdict[sampleNo][bc]+pldict[sampleNo][bc]
        ovdict[sampleNo][bc] = spdict[sampleNo][bc]+pldict[sampleNo][bc]
      else:
        totalcount[sampleNo][bc] = spdict[sampleNo][bc]
    for bc in pldict[sampleNo]:
      if bc not in spdict[sampleNo]:
        totalcount[sampleNo][bc] = pldict[sampleNo][bc]
    for bc,count in sorted(totalcount[sampleNo].items(), key=lambda kv: -kv[1]):
      if bc in spdict[sampleNo]:
        spcount = spdict[sampleNo][bc]
      else:
        spcount = 0
      if bc in pldict[sampleNo]:
        plcount = pldict[sampleNo][bc] 
      else:
        plcount = 0
      outfiles[sampleNo].write(bc+'\t'+str(int(plcount))+'\t'+str(int(spcount))+'\t'+str(int(count))+'\n')
    outfiles[sampleNo].close()
  #read treatment info
  reffile = open('analysis/PrimerID_summary.txt')
  header = reffile.readline()
  groupdict = {}
  for line in reffile:
    line = line.rstrip().rsplit('\t')
    sample = line[0].lstrip('spleen_').lstrip('plasma_')
    group = line[2]
    groupdict[sample] = group
  #print summary table
  outfile = open('analysis/overlap/'+annotation+'/plasma_spleen_'+annot2+'_barcodes_summary.tsv','w')
  outfile.write('sample\tgroup\tplasma\tspleen\ttotal\toverlap\tsimilarity\n')
  for sampleNo in samplelist:
    outfile.write(sampleNo+'\t'+groupdict[sampleNo]+'\t'+str(len(pldict[sampleNo]))+'\t'+str(len(spdict[sampleNo]))+'\t'+str(len(totalcount[sampleNo]))+'\t'+str(len(ovdict[sampleNo]))+'\t')
    sumsp = sum(spdict[sampleNo].values())
    sumpl = sum(pldict[sampleNo].values())
    similarity = 0
    freq01any = 0
    freq01both = 0
    for bc in totalcount[sampleNo]:
      if bc in spdict[sampleNo]:
        freqsp = spdict[sampleNo][bc]/sumsp
      else:
        freqsp = 0
      if bc in pldict[sampleNo]:
        freqpl = pldict[sampleNo][bc]/sumpl
      else:
        freqpl = 0
      if freqpl > 0.001 and freqsp > 0.001: freq01both += 1
      if freqpl > 0.001 or freqsp > 0.001: freq01any += 1
      similarity += sqrt(freqsp*freqpl)
    outfile.write(str(similarity)+'\n')
  outfile.close() 
  
  unionfile = open

if __name__ == '__main__':
  main()

#!/usr/bin/env python
from scipy.stats import mannwhitneyu
from scipy.stats import ttest_ind

def main():
  infile = open('analysis/PrimerID_summary.txt')
  outfile = open('analysis/PrimerID_stats.txt','w')
  header = infile.readline().rstrip().rsplit('\t')[3:]
  outfile.write('test\t'+'\t'.join(header)+'\n')
  countdict = {}
  for line in infile:
    line = line.rstrip().rsplit('\t')
    sample = line[0]
    tissue = line[1]
    group  = line[2]
    countdict[sample] = {}
    countdict[sample]['tissue'] = tissue
    countdict[sample]['group'] = group
    for i,count in enumerate(line[3:]):
      countdict[sample][header[i]] = float(count)
  infile.close()
  ttdict = readtt(countdict,header)


  #mannwhitney test plasma+spleen LRA v.s. CTL
  lradict = {}; ctldict = {}
  for sample in countdict:
    samplet = 'total_'+'_'.join(sample.rsplit('_')[1:])
    if countdict[sample]['group'] == 'DMSO':
      ctldict[sample] = ttdict[samplet]
    else:
      lradict[sample] = ttdict[samplet]
  print_test(lradict,ctldict,'mannwhitneyu_total',outfile,header)
  #mannwhitney test spleen LRA v.s. CTL
  lradict = {}; ctldict = {}
  for sample in countdict:
    if countdict[sample]['tissue'] == 'plasma': continue
    if countdict[sample]['group'] == 'DMSO':
      ctldict[sample] = countdict[sample]
    else:
      lradict[sample] = countdict[sample]
  print_test(lradict,ctldict,'mannwhitneyu_spleen',outfile,header)
  #mannwhitney test plasma LRA v.s. CTL
  lradict = {}; ctldict = {}
  for sample in countdict:
    if countdict[sample]['tissue'] == 'spleen': continue
    if countdict[sample]['group'] == 'DMSO':
      ctldict[sample] = countdict[sample]
    else:
      lradict[sample] = countdict[sample]
  print_test(lradict,ctldict,'mannwhitneyu_plasma',outfile,header)
  #mannwhitney test union (spleen and plasma) LRA v.s. CTL
  

  #t test plasma+spleen LRA v.s. CTL
  lradict = {}; ctldict = {}
  for sample in countdict:
    samplet = 'total_'+'_'.join(sample.rsplit('_')[1:])
    if countdict[sample]['group'] == 'DMSO':
      ctldict[sample] = ttdict[samplet]
    else:
      lradict[sample] = ttdict[samplet]
  print_test(lradict,ctldict,'t_test_total',outfile,header)
  #t test spleen LRA v.s. CTL
  lradict = {}; ctldict = {}
  for sample in countdict:
    if countdict[sample]['tissue'] == 'plasma': continue
    if countdict[sample]['group'] == 'DMSO':
      ctldict[sample] = countdict[sample]
    else:
      lradict[sample] = countdict[sample]
  print_test(lradict,ctldict,'t_test_spleen',outfile,header)
  #t test plasma LRA v.s. CTL
  lradict = {}; ctldict = {}
  for sample in countdict:
    if countdict[sample]['tissue'] == 'spleen': continue
    if countdict[sample]['group'] == 'DMSO':
      ctldict[sample] = countdict[sample]
    else:
      lradict[sample] = countdict[sample]
  print_test(lradict,ctldict,'t_test_plasma',outfile,header)

def readtt(countdict,header):
  ttdict = {}
  hc5file = open('analysis/overlap/HC5/plasma_spleen_PID_barcodes_summary.tsv')
  m01file = open('analysis/overlap/HC5_freq0001/plasma_spleen_PID_barcodes_summary.tsv')
  m09file = open('analysis/overlap/HC5_cumfreq90/plasma_spleen_PID_barcodes_summary.tsv')
  clfile  = open('analysis/overlap/HC5_Clustered/plasma_spleen_PID_barcodes_summary.tsv')
  c01file = open('analysis/overlap/HC5_Clustered_freq0001/plasma_spleen_PID_barcodes_summary.tsv')
  c09file = open('analysis/overlap/HC5_Clustered_cumfreq90/plasma_spleen_PID_barcodes_summary.tsv')
  hc5dict = readlist(hc5file)
  m01dict = readlist(m01file)
  m09dict = readlist(m09file)
  cldict  = readlist(clfile)
  c01dict = readlist(c01file)
  c09dict = readlist(c09file)
  for sample in countdict:
    sample  = '_'.join(sample.rsplit('_')[1:])
    samplet = 'total_'+sample
    ttdict[samplet] = {}
    if sample not in hc5dict:
      for i in header[1:]:
        ttdict[samplet][i] = 0 
      continue
    ttdict[samplet][header[1]] = hc5dict[sample]
    ttdict[samplet][header[2]] = m01dict[sample]
    ttdict[samplet][header[3]] = m09dict[sample]
    ttdict[samplet][header[4]] = cldict[sample]
    ttdict[samplet][header[5]] = c01dict[sample]
    ttdict[samplet][header[6]] = c09dict[sample]
  return ttdict
  
def readlist(infile):
  header = infile.readline()
  tmpdict = {}
  for line in infile:
    line = line.rstrip().rsplit('\t')
    tmpdict[line[0]] = float(line[4])
  infile.close()
  return tmpdict

def print_test(lradict,ctldict,testname,outfile,header):
  outfile.write(testname)
  for index in header:
    lralist = []; ctllist = []
    if 'total' in testname and 'PID_count' == index:
      outfile.write('\tNA')
      continue
    for sample in lradict:
      lralist.append(lradict[sample][index])
    for sample in ctldict:
      ctllist.append(ctldict[sample][index])
    if 'mann' in testname:
      u,p = mannwhitneyu(lralist,ctllist,use_continuity=False, alternative='two-sided')
    if 't_test' in testname:
      t,p = ttest_ind(lralist,ctllist)
    outfile.write('\t'+str(p))
  outfile.write('\n')

if __name__ == '__main__':
  main()

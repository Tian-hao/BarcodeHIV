#!/usr/bin/env python
import glob
import os
import sys

def main():
  #filtered to 5 UMI 
  infiles = sorted(glob.glob('analysis/HC5/LRA/PID_count_*.txt'))+sorted(glob.glob('analysis/HC5/CTL/PID_count_*.txt'))
  samplelist = []
  sumfile = open('analysis/PrimerID_summary_d5.txt','w')
  sumfile.write('sample\ttissue\tgroup\tPID_count\tBarcode_count\tBarcode_0001_count\tBarcode_90_count\tCluster_count\tCluster_0001_count\tCluster_90_count\n')
  PIDcount = {}
  HC5count = {}
  bcdict   = {}
  M01count = {}
  M09count = {}
  Clucount = {}
  CluM01   = {}
  CluM09   = {}
  groupdict= {}
  for infile in infiles:
    if '_R1' in infile or '_R2' in infile: continue
    sample = infile.rsplit('.txt')[0].rsplit('_')
    sample = '_'.join(sample[2:])
    samplelist.append(sample)
    if 'LRA' in infile:
      groupdict[sample] = 'SUW133'
    else:
      groupdict[sample] = 'DMSO'
    inhandle = open(infile)
    PIDcount[sample] = 0
    HC5count[sample] = 0
    Clucount[sample] = 0
    CluM01[sample]   = 0
    CluM09[sample]   = 0
    bcdict[sample]   = {}
    for line in inhandle:
      line = line.rstrip().rsplit('\t')
      if int(line[1]) >= 5: 
        HC5count[sample] += 1
        PIDcount[sample] += int(line[1])
        bcdict[sample][line[0]] = int(line[1])

  for sample in bcdict:
    m01file = open('analysis/HC5_freq0001/PID_count_'+sample+'.txt','w')
    m09file = open('analysis/HC5_cumfreq90/PID_count_'+sample+'.txt','w')
    newdict = bcdict[sample]
    sort_newdict = sorted(newdict.items(), key=lambda kv: -kv[1])
    M01count[sample] = 0; M09count[sample] = 0
    for bc,count in sort_newdict:
      if newdict[bc]/PIDcount[sample] > 0.001:
        M01count[sample] += 1
        m01file.write(bc+'\t'+str(int(newdict[bc]))+'\n')
    cumfreq = 0
    for bc,count in sort_newdict:
      freq = float(count)/PIDcount[sample]
      cumfreq += freq
      M09count[sample] += 1
      m09file.write(bc+'\t'+str(int(newdict[bc]))+'\n')
      if cumfreq > 0.9:
        break
    m01file.close()
    m09file.close()
    tmpin = open('analysis/tmpin','w')
    bcindex = 0
    tmpindexdict = {}
    for bc,count in sort_newdict:
      bcindex += 1
      tmpindexdict[bcindex] = bc
      tmpin.write('>bc'+str(bcindex)+'\n'+bc+'\n')
    tmpin.close()
    os.system('~/Tools/cdhit/cdhit-master/cd-hit -i analysis/tmpin -o analysis/tmpout -c 0.75')
    tmpfile = open('analysis/tmpout.clstr')
    Cluster = {}
    for line in tmpfile:
      if line[0] == '>': continue
      if line[0] == '0':
        key = line.rsplit('>bc')[-1].rsplit('...')[0]
        Cluster[key] = [key]
      else:
        kid = line.rsplit('>bc')[-1].rsplit('...')[0]
        Cluster[key].append(kid)
    tmpfile.close()
    clusterdict = {}
    for bc in Cluster:
      clusterdict[bc] = 0
      for kid in Cluster[bc]:
        clusterdict[bc] += sort_newdict[int(kid)-1][1]
    sort_cluster = sorted(clusterdict.items(), key=lambda kv: -kv[1])
    clusterfile = open('analysis/HC5_Clustered_d5/PID_clustered_count_'+sample+'.txt','w')
    clu01file = open('analysis/HC5_Clustered_freq0001_d5/PID_clustered_count_'+sample+'.txt','w')
    clu09file = open('analysis/HC5_Clustered_cumfreq90_d5/PID_clustered_count_'+sample+'.txt','w')
    for bc,count in sort_cluster:
      Clucount[sample] += 1
      if float(count)/PIDcount[sample] > 0.001:
        CluM01[sample] += 1
        clu01file.write(tmpindexdict[int(bc)]+'\t'+str(count)+'\n')
      clusterfile.write(tmpindexdict[int(bc)]+'\t'+str(count)+'\n')
    cumfreq = 0
    for bc,count in sort_cluster:
      freq = float(count)/PIDcount[sample]
      cumfreq += freq
      CluM09[sample] += 1
      clu09file.write(tmpindexdict[int(bc)]+'\t'+str(count)+'\n')
      if cumfreq > 0.9: 
        break
    clusterfile.close()
    clu01file.close()
    clu09file.close()
  for sample in samplelist:
    tissue = sample.rsplit('_')[0]
    group  = groupdict[sample]
    sumfile.write(sample+'\t'+tissue+'\t'+group+'\t'+str(PIDcount[sample])+'\t'+str(HC5count[sample])+'\t'+str(M01count[sample])+'\t'+str(M09count[sample])+'\t'+str(Clucount[sample])+'\t'+str(CluM01[sample])+'\t'+str(CluM09[sample])+'\n')
  sumfile.write('plasma_A16_S10\tplasma\tSUW133\t0\t0\t0\t0\t0\t0\t0\n')
  sumfile.write('plasma_A16_S15\tplasma\tSUW133\t0\t0\t0\t0\t0\t0\t0\n')
  sumfile.write('plasma_A16_S21\tplasma\tSUW133\t0\t0\t0\t0\t0\t0\t0\n')
  sumfile.write('plasma_A16_S18\tplasma\tSUW133\t0\t0\t0\t0\t0\t0\t0\n')
  sumfile.write('spleen_A16_S10\tspleen\tSUW133\t0\t0\t0\t0\t0\t0\t0\n')
  sumfile.write('spleen_A16_S15\tspleen\tSUW133\t0\t0\t0\t0\t0\t0\t0\n')
  sumfile.write('spleen_A16_S21\tspleen\tSUW133\t0\t0\t0\t0\t0\t0\t0\n')
  sumfile.write('spleen_A16_S18\tspleen\tSUW133\t0\t0\t0\t0\t0\t0\t0\n')
  sumfile.close()

if __name__ == '__main__':
  main()

#!/usr/bin/env python
import glob
from matplotlib import pyplot as plt

def main():
  d1file = open('analysis/PrimerID_summary_d1.txt')
  d2file = open('analysis/PrimerID_summary_d2.txt')
  d3file = open('analysis/PrimerID_summary.txt')
  d4file = open('analysis/PrimerID_summary_d4.txt')
  d5file = open('analysis/PrimerID_summary_d5.txt')
  countdict = {}
  header = d1file.readline()
  for line in d1file:
    line = line.rstrip().rsplit('\t')
    sample = line[0]
    countdict[sample] = []
    countdict[sample].append(int(line[4]))
    countdict[sample].append(int(line[7]))
  d1file.close()
  coundict = readsum(d2file,countdict)
  coundict = readsum(d3file,countdict)
  coundict = readsum(d4file,countdict)
  coundict = readsum(d5file,countdict)
 
  plt.figure()
  for sample in countdict:
    plt.plot(range(6),countdict[sample],color='grey',alpha=0.5)
  plt.xlabel('hamming distance of clustering')
  plt.ylabel('number of clusters')
  plt.savefig('figures/cluster_distance.png')
  plt.close('all')

def readsum(sumfile,countdict):
  header = sumfile.readline()
  for line in sumfile:
    line = line.rstrip().rsplit('\t')
    sample = line[0]
    countdict[sample].append(int(line[7]))
  sumfile.close()
  return countdict
    
  

if __name__ == '__main__':
  main()

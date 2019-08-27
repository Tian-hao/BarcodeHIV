#!/usr/bin/env python
import glob

def most_common(lst):
    return max(set(lst), key=lst.count)

def main():
  workpath = '/u/scratch/t/tianhao/HiSeq080119/'
  infiles = sorted(glob.glob(workpath+'barcode/*.txt'))
  tagfile = open('../Fasta/tag_Jocelyn_080119.txt')
  tagdict = {}
  for line in tagfile:
    line = line.rstrip().rsplit('\t')
    tagdict[line[0]] = line[1]
  tagfile.close()
  outfiles = {}
  for tag in tagdict:
    library = tagdict[tag].replace('#','_S')
    outfiles[tag] = open('../result/PrimerID_Jocelyn_batch2_'+library+'.txt','w')
  sumfile = open('../result/PrimerID_Jocelyn_batch2_depth_summary.txt','w')
  sumfile.write('sample\tdepth\n')
  bcdict = {}
  #bcdict[tag] = {UMI:[bc1,bc1,bc1,bc2,...]}
  for tag in tagdict:
    bcdict[tag] = {}
  for infile in infiles:
    print('reading '+infile)
    inhandle = open(infile)
    for line in inhandle:
      line = line.rstrip().rsplit('\t')
      tag = line[0]
      if tag not in tagdict: continue
      bc = line[1]
      UMI = line[2]
      if UMI not in bcdict[tag]: bcdict[tag][UMI] = []
      bcdict[tag][UMI].append(bc)
    inhandle.close()
  
  for tag in bcdict:
    depth = 0
    countdict = {}
    print('writing sample '+tagdict[tag])
    for UMI in bcdict[tag]:
      depth += len(bcdict[tag][UMI])
      bc = most_common(bcdict[tag][UMI])
      if bc not in countdict: countdict[bc] = 0
      countdict[bc] += 1
    for bc in countdict:
      outfiles[tag].write(bc+'\t'+str(countdict[bc])+'\n')
    outfiles[tag].close()
    sumfile.write(tagdict[tag]+'\t'+str(depth)+'\n')
  sumfile.close()
   
     
if __name__ == '__main__':
  main()

#!/usr/bin/env python

# Modified 12/29/18 RPH to:
# 1. Iterate through a sam file (obtained via a samtools view command, as piped stdin), 
#    storing matching and non-matching IDs into dict with best alignment score.  
#    Only the 10X BAM file contains the =corrected= barcode.
#    The reads matching the barcode list are assumed to be from the chosen genome 
#    and these will be selected over alternate genome.
# 3. Open and interate through paired fastq files, writing to two sets of 
#    fastq files based on stored, partitioned IDs.
#
# Original code based on snippets from JMG 12/2018

import sys,os
from Bio import SeqIO
import itertools
import datetime

def checkOptions():
  import getopt
  global sampleID
  global barcodeFile
  global sp1
  global sp2
  global sampleName
  global vString
  global oFname
  
  try:
    opts, args = getopt.getopt(sys.argv[1:], "hvg:s:a:o:", ["help", "version", "genome=","sample=","altgenome=","output="])
  except getopt.GetopError, err:
    print str(err)
    usage()
    sys.exit(2)
  sp1 = 'hg19' #default
  sp2 = 'mm10' #default
  sampleName = None #default
  vString = 'Version 0.2 - RPH 03-13-2019'
  oFname = None #default
  for o, a in opts:
    if o in ("-v", "--version"):
      print vString
      usage()
      sys.exit()
    elif o in ("-h", "--help"):
      usage()
      printhelp()
      sys.exit()
    elif o in ("-g", "--genome"):
      sp1 = a
    elif o in ("-s","--sample"):
      sampleName = a.upper()
    elif o in ("-a", "--altgenome"):
      sp2 = a
    elif o in ("-o", "--output"):
    	oFname = a
    else:
      assert False, "unhandled option"

def printhelp():
  print '''\
  The goal is to separate single-cell RNAseq fastq files into two 
  species specific subsets.  We iterate through the 
  mixed-genome bam file to select read IDs matching the extracted barcodes
  along with alignment scores, choosing the best score by genome.  
  Finally, open paired-end fastq files and output barcode-matching records 
  to one genome, non-matching to an alternate genome. Requires that 
  IDs match in paired files.
  '''

def usage():
  print '''\
  Usage: Usage: samtools view <BAM> | splitByScore.py 
  Options include:
  --version, -v    : Prints the version number
  --sample, -s     : Sample ID for insertion into file name (default None)
  --genome, -g     : Selected genome identifier (default hg19)
  --altgenome, -a  : Non-matching genome identifier (default mm10)
  --help, -h       : Display more information
  '''


def getTag(tag, lis):
  '''
  Get optional tag from a SAM record.
  '''
  for t in lis:
    spl = t.split(':')
    if spl[0] == tag:
      return spl[-1]
  return None

def main():
  '''Main.'''
  checkOptions()
  if sampleName is None:
    print "Sample name is required."
    usage()
    sys.exit(1)	
  if oFname:  #if -o flag gives output file, otherwise stdout
  	sys.stdout = open(oFname,'w')
  #print settings to log this run
  print "splitByScore run started " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
  print vString
  print "Chosen species is " + sp1 + ", alternate species is " + sp2
  
  # initialize counts
  ReadCount = 0 # total reads input
  GoodCount = 0 # total matching reads input
  sp1Count = 0 # reads assigned to sp1 
  sp2Count = 0 # reads assigned to sp2
  
  sp1IDs = dict() #dict to store current-best IDs and scores for sp1
  sp2IDs = dict() #dict to store current-best IDs and scores for sp2
  
  
  # parse SAM (stdin), extract best alignment scores, store in dict
  f = sys.stdin
  for line in f:
    if line[0] == '@':
      continue #skip header
    ReadCount += 1 #found a read line
    spl = line.rstrip().split('\t') #split into fields
    flag = int(spl[1])
    if flag & 0x4: 
      continue #skip this record because read is unmapped

    #read is mapped
    GoodCount += 1
    # retrieve alignment score
    aScore = getTag('AS', spl[11:])
    if not aScore:
      print('missing AS: ' + spl[0] + '\n')
      continue
    aScore = int(aScore)
    
    readGenome = spl[2][:len(sp1)]
    if readGenome == sp1:  #this read is from sp1
      if spl[0] in sp1IDs:  #seen this read in sp1 before
        if aScore > sp1IDs[spl[0]]:
          # this read has better AS
          sp1IDs[spl[0]] = aScore  #replace with better AS
          continue
        else:
          # previous read had equal or better AS -- do nothing
          continue
      elif spl[0] in sp2IDs: #was is seen in sp2 before?
        if aScore > sp2IDs[spl[0]]:  #this new sp1 read is better than the old sp2 read
          del sp2IDs[spl[0]]  #remove from sp2 dict
          sp2Count += -1  #decrement sp2 count
          sp1IDs[spl[0]] = aScore  #add to sp1 dict
          sp1Count += 1  #increment sp1 count
          continue
        else: # score in sp2 is better than in sp1, do nothing
          continue
      else:
        #not seen this before - add it to sp1 dict
        sp1IDs[spl[0]] = aScore
        sp1Count += 1
    elif readGenome == sp2:  #this read is from sp2
      if spl[0] in sp2IDs:  #seen this read in sp2 before
        if aScore > sp2IDs[spl[0]]:
          #this read has better score
          sp2IDs[spl[0]] = aScore  # replace with better AS
          continue
        else:
          # previous read had equal or better AS -- do nothing
          continue
      elif spl[0] in sp1IDs:  # was this seen in sp1 before?
        if aScore > sp1IDs[spl[0]]: #this new sp2 read is better than the old sp1 read
          del sp1IDs[spl[0]] #remove from sp1 dict
          sp1Count += -1
          sp2IDs[spl[0]] = aScore #add to sp2 dict
          sp2Count += 1
          continue
        else: # score in sp1 is better than in sp2 -- do nothing
          continue
      else: #not seen this before, add it to sp2 dict
        sp2IDs[spl[0]] = aScore
        sp2Count += 1


  #Write status of ID matching to stdout or redirected to oFname
  print("Total " + str(ReadCount) + " reads, " + str(GoodCount) + " are aligned")
  print(str(sp1Count) + " reads matching " + sp1)
  print(str(sp2Count) + " reads matching " + sp2)

  f.close()
  
  # The read IDs are now stored in the dicts.  Open paired-end fastq files.
  fastqs = [f for f in os.listdir(".") if f.endswith('.fastq')]
  if (len(fastqs) != 2):
    print("Fastq files not found in current directory, exiting...")
    sys.exit(1)
  R1 = ''.join([f for f in fastqs if "_R1_" in f])
  R2 = ''.join([f for f in fastqs if "_R2_" in f])
  #open both ends
  try:
    fastq_iter1 = SeqIO.parse(open('./' + R1),"fastq")
  except Exception as e:
    print("Unable to open " + R1 + " fastq file, exiting...")
    sys.exit(1)
  try:
    fastq_iter2 = SeqIO.parse(open('./' + R2),"fastq")
  except Exception as e:
    print("Unable to open " + R2 + " fastq file with error " + e + ",exitting...")
    sys.exit(1)
  
  # Create new directories for two sets of output files (sp1 and sp2)
  if (not os.path.exists('../../' + sp1)):
    os.mkdir('../../' + sp1)
  if (not os.path.exists('../../' + sp2)):
    os.mkdir('../../' + sp2)
  
  if sampleName is not None: #if specify sample ID, insert it into file name (expected by cellranger)
    R1=R1.replace("L001",sampleName+"_L001")
    R2=R2.replace("L001",sampleName+"_L001")
    # Extract subdir name from sample ID
    subdirName = "SI-GA-A" + sampleName[-1:]
  
  # Create new sample-specific sudirectories if not exist already
  if (not os.path.exists('../../' + sp1 + '/' + subdirName)):
    os.mkdir('../../' + sp1 + '/' + subdirName)
  if (not os.path.exists('../../' + sp2 + '/' + subdirName)):
    os.mkdir('../../' + sp2 + '/' + subdirName)
  
  
  # Create new output fastq files
  sp1OutFileR1 = open('../../' + sp1 + '/' + subdirName + '/' + R1,'at')
  sp1OutFileR2 = open('../../' + sp1 + '/' + subdirName + '/' + R2,'at')
  sp2OutFileR1 = open('../../' + sp2 + '/' + subdirName + '/' + R1,'at')
  sp2OutFileR2 = open('../../' + sp2 + '/' + subdirName + '/' + R2,'at')
  
  #iterate through paired reads from fastq files
  #and partition based on whether id is in sp1ID keys
  ctr=0
  sp1Ctr=0
  sp2Ctr=0
  for rec1, rec2 in itertools.izip(fastq_iter1, fastq_iter2):
    if rec1.id != rec2.id: #found mismatching records, stop and report
      print("Mismatch found: R1 ID is " + rec1.id + " R2 ID is " + rec2.id)
      print("Reporting incomplete results and stopping...")
      break
    ctr += 1
    if (rec1.id in sp1IDs): # record from chosen genotype
      sp1Ctr += 1
      SeqIO.write(rec1,sp1OutFileR1,'fastq')
      SeqIO.write(rec2,sp1OutFileR2,'fastq')
      #remove this id from dict to prevent duplicates
      del sp1IDs[rec1.id]
    elif (rec1.id in sp2IDs):  # record from alternate genotype
      sp2Ctr += 1
      SeqIO.write(rec1,sp2OutFileR1,'fastq')
      SeqIO.write(rec2,sp2OutFileR2,'fastq')
    else: # record not in either dict
      continue
    
  
  #clean up
  fastq_iter1.close()
  fastq_iter2.close()
  sp1OutFileR1.close()
  sp1OutFileR2.close()
  sp2OutFileR1.close()
  sp2OutFileR2.close()
  
  #report counts
  print("Found " + str(ctr) + " records in fastq input")
  print("Output " + str(sp1Ctr) + " records to " + sp1 + "/" + subdirName + " folder")
  print("Output " + str(sp2Ctr) + " records to " + sp2 + "/" + subdirName + " folder")
  print "splitByScore run ended " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M")

if __name__ == '__main__':
  main()

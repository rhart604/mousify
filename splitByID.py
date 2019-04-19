#!/usr/bin/env python

# Modified 12/29/18 RPH to:
# 1. Load a file of barcodes based on clustering with mixed genomes. Extract barcode
#    specific for the current sample (using argv for sample ID).
# 2. Iterate through a sam file (obtained via a samtools view command, as piped stdin), 
#    storing matching and non-matching IDs into dicts with alignment scores.  
#    Only the 10X BAM file contains the =corrected= barcode.
#    The reads matching the barcode list are assumed to be from the chosen genome 
#    and these will be selected over alternate genome.
# 3. Open and interate through paired fastq files, writing to two sets of 
#    fastq files based on stored, partitioned IDs.
# Updated 3/13/19 to require perfectly matched IDs within pairs
# Original concept based on code snippets from JMG 12/2018

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
  global vString
  global oFname
  
  try:
    opts, args = getopt.getopt(sys.argv[1:], "hvs:b:g:a:o:", ["help", "version", "sample=", "barcodes=","genome=","altgenome=","output="])
  except getopt.GetopError, err:
    print str(err)
    usage()
    sys.exit(2)
  sampleID = None
  barcodeFile = None
  sp1 = 'hg19' #default
  sp2 = 'mm10' #default
  vString = 'Version 0.95 - RPH 03-13-2019'
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
    elif o in ("-s", "--sample"):
      sampleID = a
    elif o in ("-b", "--barcodes"):
      barcodeFile = a
    elif o in ("-g", "--genome"):
      sp1 = a
    elif o in ("-a", "--altgenome"):
      sp2 = a
    elif o in ("-o", "--output"):
    	oFname = a
    else:
      assert False, "unhandled option"

def printhelp():
  print '''\
  The goal is to separate single-cell RNAseq fastq files into two 
  species specific subsets.  We import a list of barcodes that were 
  identified by clustering to a specific genome.  Iterate through the 
  mixed-genome bam file to select read IDs matching the extracted barcodes.  
  Finally, open paired-end fastq files and output barcode-matching records 
  to one genome, non-matching to an alternate genome. 
  '''

def usage():
  print '''\
  Usage: Usage: samtools view -h <BAM> | 
  splitByID2.py -s <sampleID> -b <barcode file>
  Options include:
  --version, -v    : Prints the version number
  --sample, -s     : Sample ID (required)
  --barcode, -b    : Tab-delimited text file of barcodes (required)
  --genome, -g     : Selected genome identifier (default hg19)
  --altgenome, -a  : Non-matching genome identifier (default mm10)
  --output, -o     : File to store output (redirected from stdout)
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

def loadIDs(idFile,sampStr):
  '''
  Open, read, and subset by sample string selected barcode ids in to a list
  '''
  ids=open(idFile,'r')
  hids=ids.readlines()
  hids= [s.rstrip() for s in hids]
  #split and store sample ID and barcodes separately
  hids_sample = [s.split('_',1)[0] for s in hids]
  hids_bc = [s.split('_',1)[1] for s in hids]  
  #return subset of list for this sample only
  return [hids_bc[i] for i in [i for i, s in enumerate(hids_sample) if sampStr in s]]

def main():
  '''Main.'''
  if(len(sys.argv)) < 2:
    usage()
    print 'No action specified.'
    sys.exit()
  else:
    checkOptions()
  if sampleID is None:
    print "Sample ID not specified, stopping..."
    sys.exit()
  if oFname:  #if -o flag gives output file, otherwise stdout
  	sys.stdout = open(oFname,'w')
  	
  #print settings to log this run
  print "splitByID2 run started " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
  print vString
  print "Barcodes will be loaded from " + barcodeFile
  print "Sample ID is " + sampleID
  print "Chosen species is " + sp1 + ", alternate species is " + sp2
  
  #import barcodes selected by sampleID
  bc_sp1 = loadIDs(barcodeFile,sampleID)
  bcCount = len(bc_sp1)
  
  # initialize counts
  ReadCount = 0 # total reads input
  GoodCount = 0 # total matching reads input
  sp1Count = 0 # reads matched to sp1 barcode
  sp2Count = 0 # reads not matched to sp1 barcode
  
  sp1IDs = dict() #dicts to store selected IDs for chosen genome
  
  # parse SAM (stdin), search barcodes for match in bc_sp1, store in dict
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

    #check barcode against selection list
    if getTag('CB', spl[11:]) is None: # record has no barcode; skip it
      continue
    bc = getTag('CB', spl[11:]).rstrip('-1') # assume each barcode appended with "-1"
    indices = [i for i, s in enumerate(bc_sp1) if s in bc] #should return only one
    if len(indices) > 1: #too many matches
      print("Record " + spl[0] + int(len(indices)) + " indices") # test to report n matches > 1
    if (len(indices) > 0): #one match, assume sp1
      sp1Count += 1
      sp1IDs[spl[0]] = aScore  #store id -> score 
    else: #no match, assume sp2, don't store
      sp2Count += 1
      continue
  #Write status of ID matching to stdout or redirected to oFname
  print("Barcodes found for sample " + sampleID + ": " + str(bcCount))
  print("Total " + str(ReadCount) + " reads, " + str(GoodCount) + " are aligned")
  print(str(sp1Count) + " reads matching " + sp1)
  print(str(sp2Count) + " reads not matching, assume " + sp2)

  f.close()

  # The read IDs are now stored in the dict.  Open paired-end fastq files.
  fastqs = [f for f in os.listdir(".") if f.endswith('.fastq')]
  if (len(fastqs) != 2):
    print("Fastq files not found in fastq sub-directory, exiting...")
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
  
  #should change to check if sampleID already in R1 and R2 before doing this
  if sampleID is not None: #if specify sample ID, insert it into file name (expected by cellranger)
    R1=R1.replace("L001",sampleID.upper()+"_L001")
    R2=R2.replace("L001",sampleID.upper()+"_L001")
    # Extract subdir name from sample ID
    subdirName = "SI-GA-A" + sampleID[-1:]
  
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
    else:  # assume record from alternate genotype
      sp2Ctr += 1
      SeqIO.write(rec1,sp2OutFileR1,'fastq')
      SeqIO.write(rec2,sp2OutFileR2,'fastq')
    
  
  #clean up
  fastq_iter1.close()
  fastq_iter2.close()
  sp1OutFileR1.close()
  sp1OutFileR2.close()
  sp2OutFileR1.close()
  sp2OutFileR2.close()
  
  #report counts
  print("Found " + str(ctr) + " records in fastq input")
  print("Output " + str(sp1Ctr) + " records to " + sp1 + " folder")
  print("Output " + str(sp2Ctr) + " records to " + sp2 + " folder")
  print "splitByID2 run ended " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M")

if __name__ == '__main__':
  main()

'''
This package contains helper sub-routines for vcfviz including argparse types and VCF parsing functions.
'''

import argparse as ap
import gzip
import os
import subprocess

#This dictionary contains the sequence length for Human chromosomes, build 37
SQLENS = {'1': 249250621,
          '2': 243199373,
          '3': 198022430,
          '4': 191154276,
          '5': 180915260,
          '6': 171115067,
          '7': 159138663,
          '8': 146364022,
          '9': 141213431,
          '10': 135534747,
          '11': 135006516,
          '12': 133851895,
          '13': 115169878,
          '14': 107349540,
          '15': 102531392,
          '16': 90354753,
          '17': 81195210,
          '18': 78077248,
          '19': 59128983,
          '20': 63025520,
          '21': 48129895,
          '22': 51304566,
          'X': 155270560,
          'Y': 59373566,
          'MT': 16569}

#generic functions that can be used by argparse as a variable type
def readableFile(fileName):
    '''
    A type that can be used by argparse for making sure a readable file exists and is actually readable by the user
    @param fileName - the filename to check
    @return - typically the filename; raises an ArgumentTypeError if there's an issue
    '''
    if os.path.isfile(fileName) and os.access(fileName, os.R_OK):
        return fileName
    else:
        raise ap.ArgumentTypeError("Cannot read file '%s'." % fileName)

def writableFile(fileName):
    '''
    A type that can be used by argparse for making sure a writable path is in fact writable
    @param fileName - the filename to check
    @return - typically the filename; raises an ArgumentTypeError if there's an issue
    '''
    if os.access(os.path.dirname(fileName), os.W_OK):
        return fileName
    else:        
        raise ap.ArgumentTypeError("Cannot write file '%s'." % fileName)

#VCF related functions
def getVCFChroms(filename):
    '''
    Gets a list of chromosomes in a VCF file
    @param filename - the vcf file to load, must be .vcf.gz and have a tabix index
    @return - list of chromosome names in the vcf file
    '''
    cmdFrags = ['tabix', '-l', filename]
    result = subprocess.Popen(cmdFrags, stdout=subprocess.PIPE).communicate()[0]
    pieces = result.rstrip().split('\n')
    return pieces

def getVCFColumns(filename):
    '''
    Gets a list of column names in a VCF file
    @param filename - the vcf file to load, must be .vcf.gz
    @return - list of columns in the vcf file
    '''
    #this used to be a tabix command, but that tends to be much slower and is unnecessary
    fp = gzip.open(filename, 'r')
    for l in fp:
        if l[0:6] == '#CHROM':
            retPieces = l[1:].rstrip().split('\t')
            return retPieces
        elif l[0] != '#':
            break
    
    raise Exception('VCF columns not found in "'+filename+'".')

def vcfIterator(fileName, chroms=None, sampleLabels=[]):
    '''
    This is an iterator for the vcf file
    @param fileName - the vcf filename, expected to be .vcf.bgz indicating a bgzip'ed VCF; the .tbi file will need to be there as well
    @param chroms - the list of chromosomes to pull out, if None it will iterate over ALL variants in the file in order corresponding to the order from "tabix -l <fileName>"
    @param sampleLabels - the list of sample labels required to be present
    '''
    columns = getVCFColumns(fileName)
    cInd = {}
    for i, c in enumerate(columns):
        cInd[c] = i
    
    requiredColumns = ['CHROM', 'POS', 'REF', 'ALT', 'FORMAT']+sampleLabels
    for rc in requiredColumns:
        if not cInd.has_key(rc):
            raise Exception('Missing required column "'+rc+'" in VCF file: '+fileName)
    
    if chroms == None:
        chroms = getVCFChroms(fileName)
    
    for c in chroms:
        cmdFrags = ['tabix', fileName, c]
        result = subprocess.Popen(cmdFrags, stdout=subprocess.PIPE).communicate()[0]
        varList = result.rstrip().split('\n')
        for l in varList:
            pieces = l.split('\t')
            d = {columns[i]: pieces[i] for i in xrange(0, len(columns))}
            yield d

def parseInfoField(formatCol, dataCol):
    '''
    This function converted the format and data column into a tuple with the following information
    @param formatCol - the format string from the VCF
    @param dataCol - the data column from the VCF
    @return - dict of format type to data value
    '''
    formats = formatCol.split(':')
    infoPieces = dataCol.split(':')
    ret = {}
    for x in xrange(0, len(infoPieces)):
        ret[formats[x]] = infoPieces[x]
    return ret
    
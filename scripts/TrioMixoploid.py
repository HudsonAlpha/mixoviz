#!/usr/bin/env python
'''
Usage: python TrioMixoploid.py -h

This script calculates the fraction of cells that are diploid under the assumption that it is a mixture of diploid/triploid cells and
the extra copy is inherited from the maternal line.
'''

import argparse as ap
import numpy as np
import os
import sys

import vcfutil

def getFields(infoDict):
    '''
    This function is just a helper to convert all the fields into an expected data type, and handle exceptions by setting to values
    @param infoDict - the dictionary created from vcfutil.parseInfoField(...)
    @return - tuple (gt, gq, dp , ad)
        gt - genotype (str)
        gq - genotype quality (int)
        dp - total depth, sum(ad) (int)
        ad - allele depth (list of int)
    '''
    gt = infoDict['GT']
    try:
        gq = int(infoDict['GQ'])
    except:
        gq = 0
    try:
        ad = [int(x) for x in infoDict['AD'].split(',')]
        dp = sum(ad)
    except:
        ad = []
        dp = 0
    return (gt, gq, dp, ad)

def calculateRatio(ratio0011, ratio1100):
    '''
    This function will calculate 'p' and 'e' for our two ratios.  Note: we assume that the maternal line is the source of any triploidy in
    the system of equations.
    @param ratio0011 - the ratio of ALT/TOTAL for alleles where the paternal genotype is 0/0 and the maternal genotype is 1/1
    @param ratio1100 - the ratio of ALT/TOTAL for alleles where the paternal genotype is 1/1 and the maternal genotype is 0/0
    @return - tuple (p, e)
        p - the fraction of cells that are diploid given the above ratios (1-p is the fraction triploid)
        e - the error term that allows the system of equation to be solved (i.e. reference/technical bias term)
    '''
    #logical form      ==>  standard linear alg. form for solving for p and e
    #f01+e = -p/6+2/3  ==>  p+6e = 4-6*f01
    #f10+e = p/6+1/3   ==>  -p+6e = 2-6*f10
    #this stores the constants in front of the unknown variable on the left-hand side
    systemLHS = [[1, 6],
                 [-1, 6]]
    
    #this stores the calculate constant on the right-hand side based on the observed average allelic ratios
    systemRHS = [4-6*ratio0011, 2-6*ratio1100]
    
    #calculate and return p, e
    result = np.linalg.solve(systemLHS, systemRHS)
    return result[0], result[1]

def calcTrioBiallelic(vcfFN, proband, father, mother, MIN_DEPTH, MIN_QUALITY):
    '''
    This function will scan the VCF, perform the calculations, and print a TSV output to STDOUT
    @param vcfFN - the .vcf.bgz file to parse
    @param proband - the label for the proband/child
    @param father - the label for the father to test
    @param mother - the label for the mother to test
    @param MIN_DEPTH - the minimum depth required by all variant calls to consider it
    @param MIN_QUALITY - the minimum quality required by all variant calls to consider it
    '''
    #get the chromosomes we plan to go through
    chromList = vcfutil.getVCFChroms(vcfFN)
    
    #go through each chromosome gathering the alleles with each GT combination
    totalRef0011 = []
    totalAlt0011 = []
    totalRef1100 = []
    totalAlt1100 = []
    
    #header for everything
    print '##COMMAND:'
    print '##  python '+' '.join(sys.argv)
    print '##PARAMETERS:'
    print '##  MIN_DEPTH = at least '+str(MIN_DEPTH)+' reads to include variant'
    print '##  MIN_QUALITY = at least '+str(MIN_QUALITY)+' quality score to include variant'
    print '##chrom - the chromosome tested'
    print '##diploid_frac - the fraction of cells that are diploid based on the mean statistics'
    print '##triploid_frac - the fraction of cells that are triploid based on the mean statistics'
    print '##e - the error value from the system using mean statistics, values greater than .01 may indicate an atypical sample'
    print '##diploid_frac_median - the fraction of cells that are diploid based on the median statistics'
    print '##triploid_frac_median - the fraction of cells that are triploid based on the median statistics'
    print '##e_median - the error value from the system using median statistics, values greater than .01 may indicate an atypical sample'
    print '#'+'\t'.join(['chrom', 'diploid_frac', 'triploid_frac', 'e', 'diploid_frac_median', 'triploid_frac_median', 'e_median'])
    
    #iterate through the VCF
    for chrom in chromList:
        c = chrom
        if c[0:3] == 'chr':
            c = c[3:]
        
        #go through each variants
        it = vcfutil.vcfIterator(vcfFN, chroms=[chrom], sampleLabels=[proband, father, mother])
        dv = {}
        for var in it:
            #we only check bi-allelic SNPs
            if len(var['REF']) == 1 and len(var['ALT']) == 1:
                proInfo = vcfutil.parseInfoField(var['FORMAT'], var[proband])
                patInfo = vcfutil.parseInfoField(var['FORMAT'], var[father])
                matInfo = vcfutil.parseInfoField(var['FORMAT'], var[mother])
                
                proGT, proGQ, proDP, proAD = getFields(proInfo)
                patGT, patGQ, patDP, patAD = getFields(patInfo)
                matGT, matGQ, matDP, matAD = getFields(matInfo)
                
                #make sure everything passes these user-set parameters
                if (patGQ >= MIN_QUALITY and
                    patDP >= MIN_DEPTH and
                    matGQ >= MIN_QUALITY and
                    matDP >= MIN_DEPTH and
                    proGQ >= MIN_QUALITY and
                    proDP >= MIN_DEPTH and
                    ((patGT == '0/0' and matGT == '1/1') or (patGT == '1/1' and matGT == '0/0'))):
                    
                    k = (patGT, matGT)
                    if not dv.has_key(k):
                        dv[k] = [[], []]
                    dv[k][0].append(proAD[0])
                    dv[k][1].append(proAD[1])
        
        try:
            #these are the counts we care about
            r0011, a0011 = dv[('0/0', '1/1')]
            r1100, a1100 = dv[('1/1', '0/0')]
        except:
            #one of these values doesn't exist, so we cannot perform the calculation
            print '\t'.join([str(x) for x in [c]+['--']*6])
            continue
        
        try:
            #this will raise an exception for non-autosomes
            cInt = int(c)
            
            #add these to the total for the final overall score
            totalRef0011 += r0011
            totalAlt0011 += a0011
            totalRef1100 += r1100
            totalAlt1100 += a1100
        
        except:
            #skip non-autosomes
            pass
        
        #convert for faster math
        r0011 = np.array(r0011)
        a0011 = np.array(a0011)
        r1100 = np.array(r1100)
        a1100 = np.array(a1100)
        
        #calculate the B-allele frequency for each call
        freq0011 = 1.0*a0011/(r0011+a0011)
        freq1100 = 1.0*a1100/(r1100+a1100)
        
        #calculate the ratios and then plug them into the system of equations
        p, e = calculateRatio(np.mean(freq0011), np.mean(freq1100))
        p2, e2 = calculateRatio(np.median(freq0011), np.median(freq1100))
        print '\t'.join([str(x) for x in [c, p, 1-p, e, p2, 1-p2, e2]])
    
    #calculate the overall ratios
    totalRef0011 = np.array(totalRef0011)
    totalAlt0011 = np.array(totalAlt0011)
    totalRef1100 = np.array(totalRef1100)
    totalAlt1100 = np.array(totalAlt1100)
    
    freq0011 = 1.0*totalAlt0011/(totalRef0011+totalAlt0011)
    freq1100 = 1.0*totalAlt1100/(totalRef1100+totalAlt1100)
    
    p, e = calculateRatio(np.mean(freq0011), np.mean(freq1100))
    p2, e2 = calculateRatio(np.median(freq0011), np.median(freq1100))
    print '\t'.join([str(x) for x in ['autosomes', p, 1-p, e, p2, 1-p2, e2]])
    
if __name__ == '__main__':
    #first set up the arg parser
    DESC = 'This script calculates the ratios of diploid/triploid cells in a trio under the assumption that the extra copy is inherited from the maternal line'
    p = ap.ArgumentParser(description=DESC, formatter_class=ap.RawTextHelpFormatter)
    
    #optional arguments with default
    DEFAULT_DEPTH = 20
    DEFAULT_QUAL = 20
    p.add_argument('-d', metavar='depth', dest='depth', type=int, default=DEFAULT_DEPTH, help='minimum read depth to consider a variant (default: '+str(DEFAULT_DEPTH)+')')
    p.add_argument('-q', metavar='quality', dest='quality', type=int, default=DEFAULT_QUAL, help='minimum quality to consider a variant (default: '+str(DEFAULT_QUAL)+')')
    
    #required main arguments
    p.add_argument('inputVCF', type=vcfutil.readableFile, help='a tabix-formatted VCF file to analyze (data.vcf.gz)')
    p.add_argument('proband', type=str, help='proband identifier in VCF')
    p.add_argument('father', type=str, help='father identifier in VCF')
    p.add_argument('mother', type=str, help='mother identifier in VCF')
    
    #parse the arguments
    args = p.parse_args()
    
    #run the trio B-allele plot script
    calcTrioBiallelic(args.inputVCF, args.proband, args.father, args.mother, args.depth, args.quality)
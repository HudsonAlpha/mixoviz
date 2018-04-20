#!/usr/bin/env python
'''
Usage: python BAllele.py -h

This script creates scatterplots of the B-allele frequencies, one per chromosome.
'''

import argparse as ap
import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import vcfutil
    
def plotChromosomeCalls(vcfFN, sampleLabel, outDir, MIN_DEPTH, MIN_QUAL):
    '''
    This is the primary plotting function
    @param vcfFN - the VCF filename, must be a bgzipped vcf (.vcf.gz) with a tabix index (.vcf.gz.tbi)
    @param sampleLabel - the column label in the VCF for the sample we care about
    @param outDir - the directory to save the images to, all files will be saved as <chrom>.png within that directory
    @param MIN_DEPTH - the minimum depth to include a variant in the plot
    @param MIN_QUAL - the minimum quality to include a variant in the plot
    '''
    
    #do this first just to make sure it's all good
    if not os.path.exists(outDir):
        os.makedirs(outDir)
    
    chromList = vcfutil.getVCFChroms(vcfFN)
    for chrom in chromList:
        #this is done on a per-chromosome basis
        xs = []
        ys = []
        it = vcfutil.vcfIterator(vcfFN, chroms=[chrom], sampleLabels=[sampleLabel])
        
        for var in it:
            #this function will split the format and data into an easy-to-use dictionary for us
            infoDict = vcfutil.parseInfoField(var['FORMAT'], var[sampleLabel])
            
            #we only care about bi-allelic variants
            if ',' in var['ALT']:
                continue
            
            try:
                #if any of these things fail, we don't want the variant to be included
                gq = (int(infoDict['GQ']) if infoDict['GQ'] != '.' else 0)
                refAD, altAD = [int(x) for x in infoDict['AD'].split(',')]
            except:
                #something failed above such as
                #either GQ or AD is absent
                #GQ is not an integer
                #AD is not a list of integers
                pass
            else:
                #make sure we pass the filter spec
                if (refAD+altAD >= MIN_DEPTH and
                    gq >= MIN_QUAL):
                    pos = int(var['POS'])
                    ratio = 100.0*altAD/(altAD+refAD)
                    xs.append(pos)
                    ys.append(ratio)
        
        #now we plot it
        outFN = outDir+'/'+chrom+'.png'
        chromLen = max(vcfutil.SQLENS.get(chrom, 0), (xs[-1] if len(xs) > 0 else 0))
        
        plt.figure()
        plt.scatter(xs, ys, alpha=.01)
        plt.xlabel('Position on Chromosome '+chrom)
        plt.ylabel('Call ratio (sum(AD) >= '+str(MIN_DEPTH)+' && qual >= '+str(MIN_QUAL)+')')
        plt.title(vcfFN.split('/')[-1]+' '+sampleLabel+'['+chrom+']')
        plt.xlim([0, chromLen])
        plt.ylim([0, 100])
        plt.grid()
        plt.savefig(outFN)
        plt.close()
                        
if __name__ == '__main__':
    #first set up the arg parser
    DESC = "This is a script for generating a B-allele frequencies plot per chromosome"
    p = ap.ArgumentParser(description=DESC, formatter_class=ap.RawTextHelpFormatter)
    
    #optional arguments with default
    p.add_argument('-d', metavar='depth', dest='depth', type=int, default=8, help='minimum read depth to consider a variant (default: 8)')
    p.add_argument('-q', metavar='quality', dest='quality', type=int, default=0, help='minimum quality to consider a variant (default: 0)')
    
    #required main arguments
    p.add_argument('inputVCF', type=vcfutil.readableFile, help='the input VCF files to analyze')
    p.add_argument('sample', type=str, help='the sample identifier in the vcf')
    p.add_argument('outputDir', type=str, help='the output .png file to write')
    
    #parse the arguments
    args = p.parse_args()
    
    #run the B-allele frequency script
    plotChromosomeCalls(args.inputVCF, args.sample, args.outputDir, args.depth, args.quality)

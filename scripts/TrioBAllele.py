#!/usr/bin/env python
'''
Usage: python TrioBAllele.py -h

This script creates a 3x3 multiplot of bi-allelic sites in a trio, one figure per chromosome.
'''

import argparse as ap
import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt
import numpy as np
import os

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
    #if infoDict['GT'] == './.' and infoDict.get('GQ', 0) == '.' and infoDict.get('AD', 0) == '.' and infoDict.get('DP', 0) == '.':
    #    return ('0/0', 20, 20, [20, 0])
    
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
        
def plotTrioBiallelic(vcfFN, proband, father, mother, outDir, MIN_DEPTH, MIN_QUALITY):
    '''
    This function will actually plot the figures per chromosome
    @param vcfFN - the .vcf.bgz file to parse
    @param proband - the label for the proband/child
    @param father - the label for the father to test
    @param mother - the label for the mother to test
    @param outDir - the directory to save all images to
    @param MIN_DEPTH - the minimum depth required by all variant calls to consider it
    @param MIN_QUALITY - the minimum quality required by all variant calls to consider it
    '''
    #make sure we can do this first
    if not os.path.exists(outDir):
        os.makedirs(outDir)
    
    chromList = vcfutil.getVCFChroms(vcfFN)
    
    #iterate through the VCF
    dataValues = {}
    for chrom in chromList:
        it = vcfutil.vcfIterator(vcfFN, chroms=[chrom], sampleLabels=[proband, father, mother])
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
                    proDP >= MIN_DEPTH):
                    
                    k = (chrom, patGT, matGT)
                    if not dataValues.has_key(k):
                        dataValues[k] = ([], [], [])
                    dataValues[k][0].append(int(var['POS']))
                    dataValues[k][1].append(proAD[0])
                    dataValues[k][2].append(proAD[1])
    
    #this is the order from top left to bottom right of the genotypes in the final figure
    typeOrder = ['0/0', '0/1', '1/1']
    
    #go through each chromosome gathering the alleles with each GT combination
    totalRef0011 = 0.0
    totalAlt0011 = 0.0
    totalRef1100 = 0.0
    totalAlt1100 = 0.0
    for chrom in chromList:
        c = chrom
        if c[0:3] == 'chr':
            c = c[3:]
        try:
            #this will raise an exception for non-autosomes
            cInt = int(c)
            
            #only allows autosomes to get added to these totals
            totalRef0011 += np.sum(dataValues[(chrom, '0/0', '1/1')][1])
            totalAlt0011 += np.sum(dataValues[(chrom, '0/0', '1/1')][2])
            totalRef1100 += np.sum(dataValues[(chrom, '1/1', '0/0')][1])
            totalAlt1100 += np.sum(dataValues[(chrom, '1/1', '0/0')][2])
        except Exception as e:
            pass
        
    #calculate the ratios so we can figure out what to plot
    if totalAlt0011 == 0.0:
        print 'WARNING: no 0/0 and 1/1 alleles detected'
        ratio0011 = 0.0
    else:
        ratio0011 = totalAlt0011/(totalAlt0011+totalRef0011)
    if totalAlt1100 == 0.0:
        print 'WARNING: no 1/1 and 0/0 alleles detected'
        ratio1100 = 1.0
    else:
        ratio1100 = 1-totalAlt1100/(totalAlt1100+totalRef1100)
    
    combinedRatio = .5*ratio0011+.5*ratio1100
    derivedRatio = 4-6*combinedRatio
    print 'Derived ratio=', derivedRatio
    
    #TODO: make these horizontal lines into an option
    plotHlines = [[0],
        [0, (1-derivedRatio)/3, (derivedRatio-4)/-6, 1-(derivedRatio-4)/-6],
        [(derivedRatio-4)/-6],
        [0, 1-(derivedRatio-4)/-6],
        [0, (1-derivedRatio)/3, (derivedRatio-4)/-6, 1-(derivedRatio-4)/-6, 1-(1-derivedRatio)/3, 1.0],
        [(derivedRatio-4)/-6, 1.0],
        [1-(derivedRatio-4)/-6],
        [(derivedRatio-4)/-6, 1-(derivedRatio-4)/-6, 1-(1-derivedRatio)/3, 1.0],
        [1.0]]
    
    #pltHlines =[[]]*9
    
    #this is just figuring out what to print to the screen in a tsv format 
    header = ['chrom']
    for pType in xrange(0, 3):
        for mType in xrange(0, 3):
            header += [typeOrder[pType]+'_'+typeOrder[mType], '', '']
    print '\t'.join(header)
    
    #now we can go through each chromosome and plot the results
    for chrom in chromList:
        #create a 3x3 figure
        wholefig = plt.figure()
        f, axarr = plt.subplots(3, 3, sharex=True, sharey=True)
        f.set_figheight(12)
        f.set_figwidth(12)
        plt.suptitle(vcfFN.split('/')[-1]+' '+proband+'['+chrom+']')
        
        #calculate the chromosome length
        chromLen = vcfutil.SQLENS.get(chrom, 0)
        for pType in xrange(0, 3):
            for mType in xrange(0, 3):
                k = (chrom, typeOrder[pType], typeOrder[mType])
                chromLen = max(chromLen, dataValues.get(k, ([0], [0], [0]))[0][-1])
        
        plt.xlim([0, chromLen])
        plt.ylim([0, 100])
        
        #row values stored what will eventually be printed to the screen for this chromosome
        rowValues = [chrom]
        
        for pType in xrange(0, 3):
            for mType in xrange(0, 3):
                #get the data for this chromosome and parental GTs
                k = (chrom, typeOrder[pType], typeOrder[mType])
                dv = dataValues.get(k, ([], [], []))
                
                #calculate the B-allele frequencies
                ratios = 100.0*np.array(dv[2])/(np.array(dv[1])+np.array(dv[2]))
                axarr[pType, mType].scatter(dv[0], ratios, alpha=.01)
                axarr[pType, mType].set_title(k[1]+' '+k[2])
                axarr[pType, mType].grid()
                
                try:
                    #if it is an autosome, plot the red ratio lines
                    c = chrom
                    if c[0:3] == 'chr':
                        c = c[3:]
                    c = int(c)
                    for hlineValue in plotHlines[pType*3+mType]:
                        axarr[pType, mType].axhline(100.0*hlineValue, color='red', linestyle='dashed', linewidth=2)
                except:
                    #this should only happen if int(c) fails, indicating non-autosome
                    pass
                
                #add the values to print for this chromosome
                if len(dv[0]) > 0:
                    #this method is weighting the variants by their coverage
                    refTot = np.sum(dv[1])
                    altTot = np.sum(dv[2])
                    rowValues += [refTot, altTot, 100.0*altTot/(refTot+altTot)]
                else:
                    rowValues += [0, 0, 'undefined']
        
        print '\t'.join([str(x) for x in rowValues])
                
        # hide tick and tick label of the big axes
        f.add_subplot(111, frameon=False)
        plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
        plt.xlabel("Position")
        plt.ylabel("B-allele frequency")
        
        #save and close the figure
        outFN = outDir+'/'+chrom+'.png'
        plt.savefig(outFN)
        plt.close(f)
        plt.close(wholefig)
    
if __name__ == '__main__':
    #first set up the arg parser
    DESC = 'This script creates a B-allele trio plot per chromosome for all non-biallelic SNP sites'
    p = ap.ArgumentParser(description=DESC, formatter_class=ap.RawTextHelpFormatter)
    
    #optional arguments with default
    DEFAULT_DEPTH = 20
    p.add_argument('-d', metavar='depth', dest='depth', type=int, default=DEFAULT_DEPTH, help='minimum read depth to consider a variant (default: '+str(DEFAULT_DEPTH)+')')
    p.add_argument('-q', metavar='quality', dest='quality', type=int, default=20, help='minimum quality to consider a variant (default: 20)')
    
    #required main arguments
    p.add_argument('inputVCF', type=vcfutil.readableFile, help='a tabix-formatted VCF file to analyze (data.vcf.gz)')
    p.add_argument('proband', type=str, help='proband identifier in VCF')
    p.add_argument('father', type=str, help='father identifier in VCF')
    p.add_argument('mother', type=str, help='mother identifier in VCF')
    p.add_argument('outputDir', type=str, help='the output directory')

    #parse the arguments
    args = p.parse_args()
    
    #run the trio B-allele plot script
    plotTrioBiallelic(args.inputVCF, args.proband, args.father, args.mother, args.outputDir, args.depth, args.quality)
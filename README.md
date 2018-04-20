# mixoviz
Detection and visualization of 2n/3n mixoploidy from a VCF.  This repo contains the scripts used to calculate ratios of diploid-triploid
mixoploidy within a Whole Genome Sequencing (WGS) sample.  The only input requirement is a tabix-indexed VCF file containing a trio (child, mother, father)
with alleles depths (AD tag) and call quality (GQ tag).  We have tested the scripts using Python 2.7 and VCF files from GATK3.

## Available Scripts

1. BAllele.py - Generates a B-allele plot for each chromosome for a single sample (trio not required).
2. TrioBAllele.py - Generates a 3x3 B-allele plot for each chromosome using trio information to deconvolute the variants.
3. TrioMixoploid.py - Calculates the percentage of diploid and triploid cells present in a sample under the assumption that the source of the extra haplotype is the mother.  Calculations are performed on each chromosome (i.e. mosaic trisomy) and across all autosomes (i.e. 2n/3n mixoploidy).
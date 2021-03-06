# mixoviz
Detection and visualization of 2n/3n mixoploidy from a VCF.  This repo contains the scripts used to calculate ratios of diploid-triploid
mixoploidy within a Whole Genome Sequencing (WGS) sample.  The only input requirement is a tabix-indexed VCF file containing a trio (child, mother, father)
with alleles depths (AD tag) and call quality (GQ tag).  We have tested the scripts using VCF files from GATK3 on human samples.

### Requirements

1. Python - tested with Python 3.6.1
2. matplotlib, numpy, PyVCF - available via `pip install matplotlib numpy pyvcf`
3. [tabix](http://www.htslib.org/doc/tabix.html) - part of htslib, all VCF inputs are expected as tabix-indexed VCF files

### Available Scripts
For specific instructions on each script, run `python3 <ScriptName.py> -h`.

1. BAllele.py - Generates a B-allele plot for each chromosome for a single sample (trio not required).
2. TrioBAllele.py - Generates a 3x3 B-allele plot for each chromosome using trio information to deconvolute the variants.
3. TrioMixoploid.py - Calculates the percentage of diploid and triploid cells present in a sample under the assumption that the source of the extra haplotype is the mother.  Calculations are performed on individual chromosomes (i.e. mosaic trisomy) and across all autosomes (i.e. 2n/3n mixoploidy).

### Reference

Pre-print on bioRxiv: [Holt, James M. et al. Programmatic Detection of Diploid-Triploid Mixoploidy via Whole Genome Sequencing.](https://www.biorxiv.org/content/early/2018/07/18/371468)

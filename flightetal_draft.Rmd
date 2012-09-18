# Affymetrix(r) Mismatch (MM) Probes: Useful After All

Robert M Flight, Abdallah M Eteleeb, Eric C Rouchka

## Abstract 

Affymetrix® GeneChip® microarray design define probe sets consisting of 11, 16, or 20 distinct 25 base pair (BP) probes
for determining mRNA expression for a specific gene, which may be covered by one or more probe sets.  Each probe has a
corresponding perfect match (PM) and mismatch (MM) set.  Traditional analytical techniques have either used the MM probes
to determine the level of cross-hybridization or reliability of the PM probe, or have been completely ignored.  Given the
availability of reference genome sequences, we have reanalyzed the mapping of both PM and MM probes to reference genomes
in transcript regions.  Our results suggest that depending of the species of interest, 66%-93% of the PM probes can be 
used reliably in terms of single unique matches to the genome, while a small number of the MM probes (typically less 
than 1%) could be incorporated into the analysis.  In addition, we have examined the mapping of PM and MM probes to five
different human genome projects, resulting in approximately a 70% overlap of uniquely mapping PM probes, and a subset of
51 uniquely mapping MM probes commonly found in all five projects, 24 of which are found within annotated exonic regions. 
These results suggest that individual variation in transcriptome regions provides an additional complexity to microarray
data analysis. Given these results, we conclude that the development of custom chip definition files (CDFs) should
include MM probe sequences to provide the most effective means of transcriptome analysis of Affymetrix® GeneChip® 
arrays.

**Categories and Subject Descriptors**: J.3 [Life and Medical Sciences]: Biology and Genetics

**General Terms**: Algorithms, Measurement, Theory

**Keywords**: Bioinformatics, microarray, probe set, custom definition files.

## Introduction

Oligonucleotide-based microarray technologies provide a methodology whereby a researcher can indirectly measure the
expression level of an mRNA molecule being actively transcribed under a set of conditions by labeling a cDNA fragment
that hybridizes to a complementary probe sequence specific to a particular transcript.  Since their first use on
customized cDNA arrays [1] in the mid-1990s, they have been used as the de-facto standard for measuring global 
transcriptional changes under differing conditions.  While RNA-Seq [2] may eventually supplant microarrays as the 
method of choice, a large number of microarray experiments exist that have been deposited into publicly available
repositories such as NCBI’s Gene Expression Omnibus (GEO) [3] and EBI’s ArrayExpress [4].  As a case in point, GEO
contains 32,471 series as of 9/6/2012.  The majority of the entries in GEO were performed on arrays designed by
companies such as Affymetrix®, Inc. (Santa Clara, CA), Agilent Technologies, Inc. (Santa Clara, CA), Illumina®, Inc.
(San Diego, CA), and GE Healthcare Lifesciences (Piscataway, NJ), with nearly half of the series (16,181) being
performed on various Affymetrix® arrays. 

The design of Affymetrix® GeneChip® arrays in particular provides for probe sets consisting of 11, 16, or 20 distinct
25 base pair (BP) probes, with each probe having a corresponding perfect match (PM) and mismatch (MM) probe.  The PM
and MM differ by the exchange of the complementary base at the 13th position in the probe.  While MM probes were
originally designed to account for signal in the PM resulting from non-specific cross-hybridization, they are often
underutilized or completely ignored. Mismatch probes have been explored for use in long oligonucleotide arrays as well
[5], but their utilization is limited to the Affymetrix® platform. 

Affymetrix® provides a default GeneChip® analysis package known as the Micro Array Suite 5.0 (MAS 5.0) [6] that measures the signal intensity for a particular probe pair as: 

signal = TukeyBiweight{log(PMj - MM*j)} (1)

Where MM`*` is a modified version of MM that is never bigger than the intensity value of the PM.  The motivation behind
the modified mismatch intensity MM* is to report all probe-level intensities as positive values, and to remove the
influence of the minority of probes where the MM intensity value is significantly higher than the corresponding PM
intensity.  In addition to the intensity signal, MAS 5.0 also produces a detection p-value which flags a transcript as
“P” (present), “M” (marginal), or “A” (absent) based on the reliability of the probe set based on differences between
PM and MM intensities.

Known issues in the use of PM and MM probe intensities to generate a single probe set intensity values led to the
development of other approaches, including RMA [7] and GCRMA [8]  which completely ignore the MM probes. 

With the availability of individual probe and reference genome sequences, it is possible to re-map probes based on new
sources of genome annotations. This allows custom Chip Description Files (CDFs) wherein probes are grouped into novel
probe sets based on exon, transcript, and gene level annotation [9-22]. Most notable is the effort of the BrainArray
group [10] which updates custom CDFs for a large number of Affymetrix® GeneChips® by creating probe sets based on
annotated features such as Entrez Gene [23], Ensembl transcript, Ensembl gene, and RefSeq Gene [24].  Using custom CDFs
has been shown to impact the reliability of expression analysis [10, 20-22, 25]. However, to the authors’ knowledge,
only the PM probe sequences are used when generating custom CDFs.

Based on the observation that a small, yet significant number of PM-MM probe pairs exist where the MM intensity is
significantly increased over the PM intensity, our initial inclination was that these differences in intensities were
not due to cross-hybridization or rogue probes alone.  Therefore, keeping in mind that Affymetrix® probes have been
designed according to continually evolving genome assemblies, we proceeded to analyze PM and MM probes across eight
commonly studied species (Table 1) by looking at PM and MM probes that uniquely map to the respective genome.

In addition to changing functional annotations, one potential problem area for microarray probe design is the presence
of single nucleotide polymorphisms (SNPs) within a population.  As the probes are designed using a reference genome or
transcriptome, a “one-size fits all” approach has been taken for the probes on a particular array.  However, SNPs are
known to occur relatively frequently throughout the genome, with build 137 of dbSNP [26] containing over 53.5 million
reference SNPs for the human genome. We have previously studied the effects of SNPs on Affymetrix® GeneChips® [27]
showing that a large number of SNPs lie in the areas where microarray probes have been designed.  This has been taken
into account in the BrainArray’s custom CDF files which incorporate SNP information.  To study the effects that
individual variation can play in microarray analysis, we looked at the mappings of PM and MM probes within five
distinct publicly available assemblies of human genomes.

## Methods

### Mapping of PM and MM Probes

Chromosomal-based genome assemblies were downloaded from the UCSC Goldenpath Genomes ftp server using an anonymous
login (ftp://hgdownload.cse.ucsc.edu/goldenPath/) [28] for eight commonly studied species, including C. elegans
(roundworm), D. melanogaster (fruit fly), S. cerevisiae (baker’s yeast), X. tropicalis (western clawed frog), D. rerio
(zebrafish), M. musculus (house mouse), R. norvegicus (brown Norway rat), and H. sapiens (human) (Table 1).  Genome
indices were created using bowtie-build version 0.12.8 [29] with the default parameters.  Perfect match (PM) probe
sequences for Affymetrix® GeneChips® were obtained from Bioconductor (v 2.10) probe packages, which are constructed
from data available in NetAffx [30] (Table 2) with each new Bioconductor release. Mismatch (MM) probe sequences were
constructed by replacing the 13th base in the supplied PM probe sequence with the complementary base. PM and MM probes
were aligned to the indexed genomes using bowtie version 0.12.8 [29] with the parameters –v 0 and –a which used
together will report all valid probes matching with 100% identity.

###	Generation of Exons and Overlap

Exon regions were obtained from the UCSC genome browser as BED files with an entry for each exon. mergeBed from the
bedTools suite was used to merge overlapping exons from multiple transcripts into single contiguous exons. These
merged exons were used when defining overlaps of probe alignments with an exon. Probe and exon overlaps were defined
as any type of overlap with at least 23 bases overlapping on the same strand. Overlaps were determined using Genomic
Ranges version 1.8.7 [31].

###	DNA Microarray Data

For each GeneChip®, CEL files were downloaded from GEO for 20 random samples (with the exception of S. cerevisiae (12)
and X. tropicalis (4), a the GSMs are listed in [gsmFiles.txt](gsmFiles.txt)). Probe intensities were background
corrected using the MAS background correction method implemented in Bioconductor. Depending on the application,
intensities were log (base 2), square root transformed, or used as is. 

###	Negative PM-MM Set

A PM-MM set of probes was considered to be negative if nine (two for S. cerevisiae and six for X. tropicalis) or more
samples had a negative value for the difference in the PM-MM intensities. For examination of intensity distribution,
any PM-MM pair with a negative different greater than 1000 in one or more samples was considered and examined. 

###	Probe Correlations

For each MM probe that uniquely overlapped one merged exon (designated as a true-match MM, (TMmm)), the correlation
with all other MM probes in the probe set (mm) and the correlation with all other TM probes that also mapped uniquely
to the same exon (if there were three or more other probes also mapped to the exon) was calculated (tm). 

###	Human Variation
To gain an understanding of individual variation and the unique mapping of microarray probes, five whole genome
assemblies were downloaded for the human genome [32-36] (Table 3).   Probes from the HGU133APlus2.0 Affymetrix®
GeneChip® were aligned to each of these genomes using the methods previously described for mapping PM and MM probes.

## RESULTS

### Probes Matching Genomic Locations 

Given the PM probe sequences and the inferred MM sequences, individual probes were mapped to the corresponding genome
assembly as outlined in Methods.  The percentage of perfect match probes mapping to the genome ranged from a low of
80% (X. tropicalis) to a high of 95% (D. melanogaster), with the exception of S. cerevisiae (Table 4).  It must be
noted that the lower percentage of S. cerevisiae matches (53%) is expected, as the Yeast Genome 2.0 GeneChip® contains
probes for two yeast species, S. cerevisiae and S. pombe.  

Affymetrix® probesets are given suffix definitions depending upon the uniqueness of the exemplar sequence used to
design a probe set.  A designation of `_`at indicates the probe set perfectly matches a single transcript; `_`a`_`at
probe sets only perfectly match transcripts of the same gene; _s_at perfectly match multiple transcripts for the same
gene family; and `_`x`_`at indicates the probe set is identical or highly similar to other genes.  One of the
difficulties with these designations is that it relies upon a set of annotations at a particular point in time.  

Analysis of the probes that map to the genome (Table 4) indicates that 82% to 98% of the mapped probes map uniquely to
a single genomic location.  The fact that a number of probes map to multiple locations is not to be unexpected due to
the restrictions placed on probe set design.  However, it is expected that those probes mapping to multiple locations
would not be from the `_`at class of probes.  

To determine the reliability of these probes with the fluctuation of unknown transcripts, those probes that map with
100% identity to two or more locations in the genome were considered.  While these probes typically represent less
than 10% of the total number of probes for a given GeneChip®, their classification could be important in detecting
cross hybridization.  One might expect that the greatest percentage of these would be within the `_`x`_`at and
`_`s`_`at classes.  However, as Table 5 shows, the larger genomes actually contain the greatest percentage in the
`_`at and `_`a`_`at classes, with anywhere from 18% (S. cerevisiae) to 91% (R. norvegicus) of the probes matching
multiple locations belonging to the `_`at class.  In addition, a small number of MM probes map to the genome as well. 
To better understand the effects this small set of MM probes might have on gene expression, we further reduced this to
a smaller subset where the mapping was within exon regions.  For these probes, we analyzed their signal intensities
from random samples compared to the overall distribution of PM and MM intensities, and the distribution of PM and MM
intensities within the corresponding exonic sequences (Figure 1).  As these plots indicate, MM probes mapping within
exonic regions closely follow the expression density of PM probes mapping within exonic regions, and are significantly 
shifted from the overall expression profiles of MM probes.  These results suggest that while the number of these
probes is small, they offer significant information that should not be ignored, and furthermore, can confound analyses
where MM data is incorporated.

As some of these MM probes may bind to transcripts, we further considered those MM probes that uniquely mapped to
exons (irrespective of whether the corresponding PM probe mapped zero, one or multiple times to exons or the full
genome), examining the differences in signal intensity between the MM and its associated PM.  In many cases there is a
significant negative difference in the expression level of the PM-MM pair (an example is shown in Figure 2). 

If these MM probes are grouped instead with the other probes within the transcriptional region for which they uniquely
match  (we have renamed these probes as “true match” (TM) probes since they truly match the region in the genome),
there is a much better association between the probe intensities, as shown in Figure 3.

Further analysis was performed to test the correlation of the TM probes with the expression levels of both the
annotated probe group MM probes and with the TM-mapped transcript probes (Figure 4).  The box plot in Figure 4 clearly
indicates that the TM intensities more closely correlate with those from the group based on mapping to the same exon.

To determine if the observed difference in the correlations from Figure 4 is due to measurement of different mRNA
entities, we considered the MM probe both within its annotated location as well as within the new mapped location
(TMmm). The intensity of the MM probe was compared with the intensity of its corresponding neighbor probes (MM probes
for the annotated location; PM probes for the new mapped location).  An average correlation value between intensities
was calculated on a per-exon probe basis.  The resulting correlations are summarized in Table 8.

What is interesting is that with few exceptions, the transcripts and genes being measured by the probe sets are the
same, implying that the TMmm probe aligns to the same gene as its complimentary PM probe. Upon further examination it
appears that many of the PM complements of the TMmm do not align to the genome at all (data not shown). One possibili
ty is that these probes are in regions that have seen changes in the reference sequence over the years, or that the
original sequencing of the ESTs used to design the probe sets was of poor quality. 

###	Effects of Individual Variation 

To gain an understanding of the effect of individual variation, the unique mapping of PM and MM probes to five
distinct human genome assemblies was analyzed (Table 6).  Four of the five projects have roughly the same number of
uniquely mapped PM probes (within 2% variation).  The fifth project (BGI) provides an exception to this trend.  While
there are a number of potential explanations for this (including sequence and assembly quality and coverage of the
sequencing), one potential feature to be considered is the fact that this sequencing project involves the sequencing
of an African individual, and it is the only project not to have a large component of the library consisting of
Caucasian individuals.

While there is a large agreement for the number of probes, we also checked if the probes represented were consistent
among all of the projects.  A Venn diagram depicting the number of overlapping perfectly matching probes is given in
Figure 5.  As can be seen from this figure, a total of 422,279 probes uniquely map for all five assemblies.  Of these,
422,119 are perfect match probes, indicating that 81% of the HGU133APlus2.0 perfect match probes are reliable in terms
of their mapping to the genome for these assemblies.  One of the interesting results is that there are 51 shared
perfectly matching mismatch probes.  Of these, 24 fall within RefSeq annotated exonic regions (results not shown),
with 22 of the 24 showing higher correlation in the TMmm assignments calculated in Table 6.

## Conclusion

MM probes are theoretically designed to capture background and non-specific binding. Alignment of the MM probes to the
genome shows that in a very small percentage of cases, MM probes align uniquely to the genome in transcribed regions.
Signal from these probes should be useful for quantifying true transcriptional events rather than for PM signal
adjustment.

In addition, current custom CDF generation workflows ignore the MM probes during the probe alignment process. Given
that some MM probes align to reference genomes, they should be considered for inclusion when creating custom CDFs. The
utility of the probes may be limited due to variation among individuals. 

## Acknowledgements

This work was partially funded by National Institutes of Health (NIH) grant 8P20GM103436-12.  Its contents are solely
the responsibility of the authors and do not represent the official views of NIH or the National Institute of General
Medical Sciences.

% Affymetrix MisMatch (MM) Probes: </br> Useful After All
% Robert M Flight; Abdallah Eteleeb; Eric C Rouchka
% 14/12/12











## Transcriptomics

<img src="figure/dna_processing.svg" width="750"/>

## Transcriptomics

<img src="figure/dna_processing_highlight.svg" width="750"/>

# Affymetrix GeneChips


## GeneChip

- measures abundance of RNA **transcripts**
- 25 mer oligonucleotides on solid support
- oligonucleotides organized into **probesets**
- **probeset** consists of 11, 16, 20, 25 probe pairs
- each probe pair consists of a perfect match (PM) and mis-match (MM) probe

## GeneChip & ProbeSet

<img src="figure/insertChipImage.png"  height="450"  alt="chip" title="chip" /> 


<img src="figure/insertProbeSet1.svg" width="300"   alt="probeset" title="probeset" /> 



## Perfect Match & MisMatch

- Perfect Match (PM)
    - supposed to perfectly match the sequence of interest
    - has exact complementarity
    - binds **perfectly**
- MisMatch (MM)
    - 13th base is reverse complement of **PM** sequence
    - supposed to account for non-specific binding in the **PM**
    - therefore should have lower signal
    - useful for **PM** signal correction


## 




<pre><code>
pm CACCCAGCTGGT<font color="#0000FF">C</font>CTGTGGATGGGA</br>mm CACCCAGCTGGT<font color="#0000FF">G</font>CTGTGGATGGGA

</pre></code>


## Perfect Match & MisMatch

- True Signal
    - should be PM - MM
- But ...
    - MM may have higher signal than PM
    - most modern summarization methods ignore it


## CDF: Chip Definition File

- Defines organization of probes into probesets
    - genes
    - gene families
    - transcripts

### Defined by Affymetrix based on available annotations


## 

<pre><code style="font-size:10pt">
[Unit43914_Block1]
Name=243114_at
BlockNumber=1
NumAtoms=11
NumCells=22
StartPosition=0
StopPosition=10
CellHeader=X	Y	PROBE	FEAT	QUAL	EXPOS	POS	CBASE	PBASE	TBASE	ATOM	INDEX	CODONIND	CODON	REGIONTYPE	REGION
Cell1=928	614	N	control	243114_at	0	13	C	C	C	0	715624	-1	-1	99	
Cell2=928	613	N	control	243114_at	0	13	C	G	C	0	714460	-1	-1	99	
Cell3=185	304	N	control	243114_at	1	13	A	A	A	1	354041	-1	-1	99	
Cell4=185	303	N	control	243114_at	1	13	A	T	A	1	352877	-1	-1	99	
Cell5=1102	372	N	control	243114_at	2	13	A	A	A	2	434110	-1	-1	99	
Cell6=1102	371	N	control	243114_at	2	13	A	T	A	2	432946	-1	-1	99	
Cell7=129	770	N	control	243114_at	3	13	C	C	C	3	896409	-1	-1	99	
Cell8=129	769	N	control	243114_at	3	13	C	G	C	3	895245	-1	-1	99	
Cell9=375	556	N	control	243114_at	4	13	A	A	A	4	647559	-1	-1	99	
Cell10=375	555	N	control	243114_at	4	13	A	T	A	4	646395	-1	-1	99	
Cell11=526	364	N	control	243114_at	5	13	C	C	C	5	424222	-1	-1	99	
Cell12=526	363	N	control	243114_at	5	13	C	G	C	5	423058	-1	-1	99	
Cell13=1059	92	N	control	243114_at	6	13	A	A	A	6	108147	-1	-1	99	
Cell14=1059	91	N	control	243114_at	6	13	A	T	A	6	106983	-1	-1	99	
Cell15=1098	186	N	control	243114_at	7	13	C	C	C	7	217602	-1	-1	99	
Cell16=1098	185	N	control	243114_at	7	13	C	G	C	7	216438	-1	-1	99	
Cell17=399	565	N	control	243114_at	8	13	G	C	G	8	658059	-1	-1	99	
Cell18=399	566	N	control	243114_at	8	13	G	G	G	8	659223	-1	-1	99	
Cell19=1062	242	N	control	243114_at	9	13	A	A	A	9	282750	-1	-1	99	
Cell20=1062	241	N	control	243114_at	9	13	A	T	A	9	281586	-1	-1	99	
Cell21=511	580	N	control	243114_at	10	13	A	A	A	10	675631	-1	-1	99	
Cell22=511	579	N	control	243114_at	10	13	A	T	A	10	674467	-1	-1	99	
</code></pre>

## Custom CDF

- Reorganize probes 
    - available probe sequences & chip locations
    - available target genome sequences
    - available genome annotations (genes, transcripts, exons, etc)
- Generate probe sets
    - have perfect alignment
    - all probes bind the same genomic element (exon, transcript, gene)
    - binding region does not encompass a lot of SNPs

## Better Results

- Reorganization gives results that are:
    - more reproducible
    - more consistent
    
- Dai et al, 2005, NAR 33(20):e175 
- Sandberg & Larsson, 2007, BMC Bioinf. 8:48
- 41 `BrainArray` platforms on GEO
    - http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/genomic_curated_CDF.asp


## But...

- only uses perfect match (PM) probes!
- what about mismatch probes (MM)?
- why not?
    - have sequences
    - matter of aligning to genome


# MM Probe Alignment


## Organisms




<!-- html table generated in R 2.15.1 by xtable 1.7-0 package -->
<!-- Mon Dec 10 15:49:33 2012 -->
<TABLE style="border-spacing:20px 5px;">
<TR> <TH> Organism </TH> <TH> Reference Assembly </TH> <TH> Build Date </TH>  </TR>
  <TR> <TD> *C elegans* </TD> <TD> ce6 </TD> <TD> May 2008 </TD> </TR>
  <TR> <TD> *D melanogaster* </TD> <TD> dm3 </TD> <TD> Apr. 2006 </TD> </TR>
  <TR> <TD> *S cerevisiae* </TD> <TD> sc3 </TD> <TD> Apr. 2011 </TD> </TR>
  <TR> <TD> *X tropicalis* </TD> <TD> xt3 </TD> <TD> Nov 2009 </TD> </TR>
  <TR> <TD> *D rerio* </TD> <TD> dr6 </TD> <TD> Dec. 2008 </TD> </TR>
  <TR> <TD> *M musculus* </TD> <TD> mm10 </TD> <TD> Dec. 2011 </TD> </TR>
  <TR> <TD> *R norvegicus* </TD> <TD> rn4 </TD> <TD> Nov 2004 </TD> </TR>
  <TR> <TD> *H sapiens* </TD> <TD> hg19 </TD> <TD> Feb. 2009 </TD> </TR>
   </TABLE>


## GeneChips




<!-- html table generated in R 2.15.1 by xtable 1.7-0 package -->
<!-- Mon Dec 10 15:49:33 2012 -->
<TABLE style="border-spacing:20px 5px;">
<TR> <TH> Organism </TH> <TH> GeneChip </TH>  </TR>
  <TR> <TD> *C elegans* </TD> <TD> C. elegans Genome </TD> </TR>
  <TR> <TD> *D melanogaster* </TD> <TD> Drosophila Genome 2.0 </TD> </TR>
  <TR> <TD> *S cerevisiae* </TD> <TD> Yeast Genome 2.0 </TD> </TR>
  <TR> <TD> *X tropicalis* </TD> <TD> X. tropicalus Genome </TD> </TR>
  <TR> <TD> *D rerio* </TD> <TD> Zebrafish Genome </TD> </TR>
  <TR> <TD> *M musculus* </TD> <TD> Mouse Genome 430 2.0 </TD> </TR>
  <TR> <TD> *R norvegicus* </TD> <TD> Rat Genome 230 2.0 </TD> </TR>
  <TR> <TD> *H sapiens* </TD> <TD> Human Genome U133 Plus 2.0 </TD> </TR>
   </TABLE>


## Microarray Data

- random data from gene expression omnibus (GEO)
- 20 random CEL files for each organism
    - Only 4 for *X. tropicalis*, 12 for Yeast
    
## Probe Sequences & Alignment

- Sequences:
    - from **Bioconductor** `probe` packages
    - MM sequences generated from PM sequences
- Alignments:
    - align both PM and MM to reference using `bowtie v0.12.8`
    - report **all** alignments with **0** mismatches
    
# Results

## Number of Alignments




<!-- html table generated in R 2.15.1 by xtable 1.7-0 package -->
<!-- Mon Dec 10 15:49:33 2012 -->
<TABLE style="font-size:70%; text-align:left; border-spacing:20px 5px;">
<TR> <TH> Organism </TH> <TH> Number of Probe Pairs </TH> <TH> PM Mapped to Reference </TH> <TH> MM Mapped to Reference </TH> <TH> PM Unique </TH> <TH> MM Unique </TH>  </TR>
  <TR> <TD> *Ce* </TD> <TD align="right"> 249165 </TD> <TD align="right"> 226856 </TD> <TD align="right"> 143 </TD> <TD align="right"> 213745 </TD> <TD align="right">  96 </TD> </TR>
  <TR> <TD> *Dm* </TD> <TD align="right"> 265400 </TD> <TD align="right"> 251602 </TD> <TD align="right">  89 </TD> <TD align="right"> 245712 </TD> <TD align="right">  54 </TD> </TR>
  <TR> <TD> *Dr* </TD> <TD align="right"> 249752 </TD> <TD align="right"> 200608 </TD> <TD align="right"> 1282 </TD> <TD align="right"> 171282 </TD> <TD align="right"> 726 </TD> </TR>
  <TR> <TD> *Hs* </TD> <TD align="right"> 604258 </TD> <TD align="right"> 562673 </TD> <TD align="right"> 1094 </TD> <TD align="right"> 521642 </TD> <TD align="right"> 608 </TD> </TR>
  <TR> <TD> *Mm* </TD> <TD align="right"> 496468 </TD> <TD align="right"> 456674 </TD> <TD align="right"> 557 </TD> <TD align="right"> 427920 </TD> <TD align="right"> 394 </TD> </TR>
  <TR> <TD> *Rn* </TD> <TD align="right"> 342410 </TD> <TD align="right"> 304646 </TD> <TD align="right"> 391 </TD> <TD align="right"> 286784 </TD> <TD align="right"> 282 </TD> </TR>
  <TR> <TD> *Sc* </TD> <TD align="right"> 120855 </TD> <TD align="right"> 63731 </TD> <TD align="right">   1 </TD> <TD align="right"> 61942 </TD> <TD align="right">   1 </TD> </TR>
  <TR> <TD> *Xt* </TD> <TD align="right"> 648548 </TD> <TD align="right"> 519177 </TD> <TD align="right"> 1884 </TD> <TD align="right"> 426237 </TD> <TD align="right"> 1014 </TD> </TR>
   </TABLE>



## Types of Alignments

Which types of PM probes align to multiple locations?




<!-- html table generated in R 2.15.1 by xtable 1.7-0 package -->
<!-- Mon Dec 10 15:49:36 2012 -->
<TABLE style="font-size:70%; text-align:left; border-spacing:20px 5px;">
<TR> <TH> Organism </TH> <TH> _x_at </TH> <TH> _s_at </TH> <TH> _a_at </TH> <TH> _at </TH> <TH> control </TH>  </TR>
  <TR> <TD> *Ce* </TD> <TD> 4040 (31%) </TD> <TD> 6416 (49%) </TD> <TD> 0 (0%) </TD> <TD> 2465 (19%) </TD> <TD> 190 (1.4%) </TD> </TR>
  <TR> <TD> *Dm* </TD> <TD> 361 (6.1%) </TD> <TD> 2742 (47%) </TD> <TD> 224 (3.8%) </TD> <TD> 2520 (43%) </TD> <TD> 43 (0.73%) </TD> </TR>
  <TR> <TD> *Dr* </TD> <TD> 1376 (4.7%) </TD> <TD> 374 (1.3%) </TD> <TD> 1203 (4.1%) </TD> <TD> 25879 (88%) </TD> <TD> 494 (1.7%) </TD> </TR>
  <TR> <TD> *Hs* </TD> <TD> 9402 (23%) </TD> <TD> 10634 (26%) </TD> <TD> 803 (2%) </TD> <TD> 19961 (49%) </TD> <TD> 231 (0.56%) </TD> </TR>
  <TR> <TD> *Mm* </TD> <TD> 4075 (14%) </TD> <TD> 2920 (10%) </TD> <TD> 4893 (17%) </TD> <TD> 16703 (58%) </TD> <TD> 163 (0.57%) </TD> </TR>
  <TR> <TD> *Rn* </TD> <TD> 440 (2.5%) </TD> <TD> 394 (2.2%) </TD> <TD> 805 (4.5%) </TD> <TD> 16127 (90%) </TD> <TD> 96 (0.54%) </TD> </TR>
  <TR> <TD> *Sc* </TD> <TD> 81 (4.5%) </TD> <TD> 1126 (63%) </TD> <TD> 0 (0%) </TD> <TD> 271 (15%) </TD> <TD> 311 (17%) </TD> </TR>
  <TR> <TD> *Xt* </TD> <TD> 14092 (15%) </TD> <TD> 12086 (13%) </TD> <TD> 34877 (38%) </TD> <TD> 31754 (34%) </TD> <TD> 131 (0.14%) </TD> </TR>
   </TABLE>


## Comparison of PM and MM Signals




Signal density for all probes and those that align to known exons.

<img src="figure/insertDensity.svg"  height="325"  alt="plot of chunk insertDensity" title="plot of chunk insertDensity" /> 


## Comparison of PM and MM Signals

- MM probes generally show lower signal intensity
- MM probes in exons tend toward the signal of PM probes in exons
- Although few in number, these would confound any analysis depending on PM - MM

## PM MM Negative Difference

<img src="figure/pmmmDensity.svg"  height="325"  alt="Probe Set Density" title="Probe Set Density" /> 


<img src="figure/pmmmIntensity.svg" width="300"   alt="Probe Set Intensity" title="Probe Set Intensity" /> 



## TM Signal 

- compare signal of exon matching MM with other probes in that exon
- matches much better than with the other MM probes
- call these "true match" probes, because they are based on alignment
- is this a general phenomenom?

<img src="figure/tmIntensity.svg" width="400"   alt="True Match Intensity Comparison" title="True Match Intensity Comparison" /> 




## MM vs TM Correlation

- compare the correlation of TM with other MM or with TM on same exon

<img src="figure/correlationBoxPlot.svg"  height="350"  alt="Correlation of TM with MM or other TM" title="Correlation of TM with MM or other TM" /> 



## Different Transcripts?




<!-- html table generated in R 2.15.1 by xtable 1.7-0 package -->
<!-- Mon Dec 10 15:50:11 2012 -->
<TABLE style="font-size:50%; text-align:left; border-spacing:20px 5px;">
<TR> <TH> tm </TH> <TH> mm </TH> <TH> Annotated RefSeq </TH> <TH> Exon RefSeq </TH> <TH> Annotated Symbol </TH> <TH> Exon Symbol </TH>  </TR>
  <TR> <TD align="right"> 0.87 </TD> <TD align="right"> 0.53 </TD> <TD> NM_001164750, NM_001164751, NM_001164752, NM_001164753, NM_001164754, NM_001164755, NM_001164756, NM_004318, NM_020164, NM_032466, NM_032467, NM_032468 </TD> <TD> NM_032466, NM_032468, NM_001164755, NM_001164754, NM_001164753, NM_001164752, NM_001164751 </TD> <TD> ASPH </TD> <TD> ASPH </TD> </TR>
  <TR> <TD align="right"> 0.86 </TD> <TD align="right"> 0.61 </TD> <TD> NM_000898 </TD> <TD> NM_000898 </TD> <TD> MAOB </TD> <TD> MAOB </TD> </TR>
  <TR> <TD align="right"> 0.84 </TD> <TD align="right"> 0.60 </TD> <TD> NM_005328 </TD> <TD> NM_005328 </TD> <TD> HAS2 </TD> <TD> HAS2 </TD> </TR>
  <TR> <TD align="right"> 0.81 </TD> <TD align="right"> 0.55 </TD> <TD> NM_001173487, NM_001173488, NM_017544 </TD> <TD> NM_001173488, NM_001173487, NM_017544 </TD> <TD> NKRF </TD> <TD> NKRF </TD> </TR>
  <TR> <TD align="right"> 0.81 </TD> <TD align="right"> 0.78 </TD> <TD> NM_014390 </TD> <TD> NM_014390 </TD> <TD> SND1 </TD> <TD> SND1 </TD> </TR>
  <TR> <TD align="right"> 0.79 </TD> <TD align="right"> 0.82 </TD> <TD> NM_213600, NR_033151 </TD> <TD> NM_213600, NR_033151 </TD> <TD> PLA2G4F </TD> <TD> PLA2G4F </TD> </TR>
   </TABLE>


- Is the MM annotated transcript different than the mapped TM transcript (hg19)?
    - MM and TM map to the same transcript!
    - MM actually perfectly matches its transcript!
    

# Variation

## Data

- Five human genome assemblies

<!-- html table generated in R 2.15.1 by xtable 1.7-0 package -->
<!-- Mon Dec 10 15:50:12 2012 -->
<TABLE style="font-size:90%; text-align:left; border-spacing:20px 5px;">
<TR> <TH> Name </TH> <TH> Abbr </TH> <TH> Assembly </TH> <TH> Bioproject </TH> <TH> Race </TH>  </TR>
  <TR> <TD> GRCh37 </TD> <TD> Hg19 </TD> <TD> 420368 </TD> <TD> 31257 </TD> <TD> Mixed </TD> </TR>
  <TR> <TD> HS_Celera_WGSA </TD> <TD> Celera </TD> <TD> 281338 </TD> <TD> 1431 </TD> <TD> Mixed </TD> </TR>
  <TR> <TD> HuRefPrime </TD> <TD> JCVI </TD> <TD> 281188 </TD> <TD> 19621 </TD> <TD> Caucasian </TD> </TR>
  <TR> <TD> BIGAF </TD> <TD> BGI </TD> <TD> 165398 </TD> <TD> 42201 </TD> <TD> African </TD> </TR>
  <TR> <TD> HsapALLPATHS1 </TD> <TD> HSAP1 </TD> <TD> 238948 </TD> <TD> 59877 </TD> <TD> Caucasian </TD> </TR>
   </TABLE>


## Mapping Differences

<!-- html table generated in R 2.15.1 by xtable 1.7-0 package -->
<!-- Mon Dec 10 15:50:39 2012 -->
<TABLE style="border-spacing:20px 5px;">
<TR> <TH> Assembly </TH> <TH> Total </TH> <TH> Perfect Match </TH> <TH> Mismatch </TH>  </TR>
  <TR> <TD> Hg19 </TD> <TD align="right"> 522250 </TD> <TD align="right"> 521642 </TD> <TD align="right"> 608 </TD> </TR>
  <TR> <TD> Celera </TD> <TD align="right"> 515111 </TD> <TD align="right"> 514518 </TD> <TD align="right"> 593 </TD> </TR>
  <TR> <TD> JCVI </TD> <TD align="right"> 530213 </TD> <TD align="right"> 529569 </TD> <TD align="right"> 644 </TD> </TR>
  <TR> <TD> BGI </TD> <TD align="right"> 469973 </TD> <TD align="right"> 469714 </TD> <TD align="right"> 259 </TD> </TR>
  <TR> <TD> HSAP1 </TD> <TD align="right"> 522480 </TD> <TD align="right"> 521922 </TD> <TD align="right"> 558 </TD> </TR>
   </TABLE>



## Mapping Differences

<img src="figure/vennFigureDiffs.svg"  height="350"  alt="plot of chunk vennFigureDiffs" title="plot of chunk vennFigureDiffs" /> 





- 161 shared MM probes
- 24 fall within RefSeq annotated exonic regions
- 16 show higher correlation with TM than MM


## Useful!!

- MM probes **theoretically** capture non-specific binding
- Small percentage appear able to capture true transcriptional events
- Therefore, custom CDF workflows should include the MM probe sequences
    - `Bowtie` makes this lightning fast compared to `Blat`

## Tools Used

- Produced using:
    - `knitr`
    - `pandoc`
    - `bzslides`
    
- Source code on Github
    - github.com/rmflight/affyMM

- Draft publication on Github
    - rmflight.github.com/affyMM
    
    
## Questions??
    
## Session Info

R version 2.15.1 (2012-06-22)
Platform: x86_64-pc-mingw32/x64 (64-bit)

locale:
[1] LC_COLLATE=English_United States.1252 
[2] LC_CTYPE=English_United States.1252   
[3] LC_MONETARY=English_United States.1252
[4] LC_NUMERIC=C                          
[5] LC_TIME=English_United States.1252    

attached base packages:
[1] grid      stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] GenomicRanges_1.8.13    xtable_1.7-0           
 [3] Biostrings_2.24.1       IRanges_1.14.4         
 [5] affy_1.34.0             hgu133plus2probe_2.10.0
 [7] hgu133plus2cdf_2.10.0   hgu133plus2.db_2.7.1   
 [9] org.Hs.eg.db_2.7.1      RSQLite_0.11.2         
[11] DBI_0.2-5               AnnotationDbi_1.18.3   
[13] Biobase_2.16.0          BiocGenerics_0.2.0     
[15] VennDiagram_1.5.1       plyr_1.7.1             
[17] ggplot2_0.9.2.1         knitr_0.8              

loaded via a namespace (and not attached):
 [1] affyio_1.24.0         BiocInstaller_1.4.9   colorspace_1.1-1     
 [4] dichromat_1.2-4       digest_0.5.2          evaluate_0.4.2       
 [7] formatR_0.6           gtable_0.1.1          labeling_0.1         
[10] MASS_7.3-18           memoise_0.1           munsell_0.4          
[13] preprocessCore_1.18.0 proto_0.3-9.2         RColorBrewer_1.0-5   
[16] reshape2_1.2.1        scales_0.2.2          stats4_2.15.1        
[19] stringr_0.6.1         tools_2.15.1          zlibbioc_1.2.0       


## 

require(knitr)
knit("rmflight_biomedCom_2012.Rmd")
system("pandoc -s -S -t dzslides --slide-level=2 --mathjax rmflight_biomedCom_2012.md -o rmflight_biomedCom_2012.html")

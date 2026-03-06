# Changelog

## v2.0.0

* Finally removed the dependency to the R program all together! This means
  now Stargazer can be run without R and is a 100% Python program (except
  for the Beagle program which was written in Java). As the replacement of
  R, Stargazer relies on four of the most basic Python modules: numpy,
  pandas, scipy and matplotlib.
  [#13](https://github.com/NickersonGenomeSciUW/Stargazer/issues/13)
* From now on, by default Stargazer will not output profiles for read depth,
  copy number, and allele fraction. Use the `-p/--include_profiles` argument
  to include those profiles.
  [#19](https://github.com/NickersonGenomeSciUW/Stargazer/issues/19)
* Extended Stargazer to CYP4A11, CYP17A1, PTGIS (formerly CYP8A1),
  CYP4A22, ABCB1, VKORC1, 2C_CLUSTER, XPC.
  [#23](https://github.com/NickersonGenomeSciUW/Stargazer/issues/23)
* Added a new command line which will make Stargazer run in GDF-creation
  mode. Stargazer switches to gdf-creation mode with'-G/--create-gdf-file' 
  argument. You must provide both '-c/--control-gene' and '-B/--bam-file'.
* Added a new reporting feature: Stargazer  creates a report that shows
  the phased haplotypes, but also the other possible candidates to consider 
  such as unphased, multiple star alleles in one haplotype, warns against 
  BEAGLE imputed missing genotypes.

## v1.0.8

* Official release on February 25, 2020.
* Extended to CACNA1S, CFTR, IFNL3, RYR1, SLCO1B3, SULT1A1, and UGT1A4.
* Implemented a new system for parsing command-line arguments and providing
  help messages. As a result, some tools and their arguments have been
  changed. Enter 'python3 stargazer.py --help' for details.
* Added a new subtool called 'sdf2gdf' under the 'setup' tool. This subtool
  takes in a SDF (samtools depth format) file and outputs a GDF
  (GATK-DepthOfCoverage format) file. This subtool was designed for the
  users of GATK4 who rely on 'samtools depth' instead of DepthOfCoverage
  from GATK3 to obtain per-base read depth because DepthOfCoverage is not
  supported by GATK4 currently. Enter 'python3 stargazer.py setup --help'
  for details.
* Added a new subtool called 'sgea' under the 'pipeline' tool. This subtool
  allows per-project genotyping with SGE and Anaconda (GATK4-based). Enter
  'python3 stargazer.py pipeline --help' for details.
* Added a new subtool called 'changepoint' under the 'view' tool. This
  subtool can detect changepoints (breakpoints) in copy number. Enter
  'python3 stargazer.py view --help' for details.
* Updated to plot read depth profiles, in addition to usual copy number
  and allele fraction profiles.
* Reorganized the code so that SV detection, including GDF file processing
  for CN analysis, is now performed solely with R, instead of the
  combination of Python and R. One major benefit of this change is the
  significant increase in computing speed. This version is on average
  about 25% faster than the previous version.
* Updated the phase-by-extension algorithm to handle multiallelic sites
  in addition to usual biallelic sites. As an example, consider three star
  alleles commonly found in the UGT1A1 gene: \*28, \*36, and \*37. The \*28
  allele has one core variant (234668879:CAT>CATAT) and three tag variants,
  the \*36 allele has one core variant (234668879:CAT>C) and four tag
  variants, and finally the \*37 allele has one core variant
  (234668879:CAT>CATATAT) and four tag variants. Notice how all three core
  variants belong to the same multiallelic site. In order to phase an input
  VCF file containing these variants, the reference panel must contain
  phase information for all the variants; however, the 1000 Genomes Project
  has phase information only for 234668879:CAT>CATAT. The previous version
  of Stargazer could not phase this site either because the
  phase-by-extension algorithm was restricted to biallelic sites. But this
  restriction is now gone, and Stargazer can call UGT1A1 genotypes such
  as \*28/\*37.
* Updated the SNP table and the star allele table for several existing genes.
* Fixed some minor bugs.

## v1.0.7

* Official release on November 20, 2019.
* Extended to CYP1B1, CYP4B1, CYP19A1, and TBXAS1 (formerly CYP5A1).
* This version uses a new algorithm for haplotype phasing rare variants
  that are not present in reference VCF -- i.e. cannot be phased
  statistically: 1) The algorithm is termed as "phase-by-extension".
  2) The algorithm does not replace statistical phasing; it's only
  supplementary. 3) The algorithm utilizes haplotype information obtained
  by statistical phasing and a scoring system to determine which of the
  two haplotypes in a given sample is more likely to carry the rare variant
  based on the total number of "tagging" SNPs relevant to that rare variant
  and matching the observed haplotype -- hence the algorithm's name
  "phasing by extension." 4) Take the non-functional CYP2D6*21 allele as an
  example, which is defined by 2580_2581insC (core), 2851C>T (tagging),
  and 4181G>C (tagging). 2851C>T and 4181G>C are present in the 1KGP panel
  and thus statistically phasable, while 2580_2581insC is not. In order to
  call a sample with 2580_2581insC as having CYP2D6*21, Stargazer will first
  check which of the haplotypes contains 2851C>T and 4181G>C and then
  assign 2580_2581insC to that haplotype.
* Added a new tool named 'view'. Briefly, this tool allows you to view,
  for example, the list of SNPs and their information of a given star
  allele for specific samples (e.g. genotype call, allelic depth,
  functional annotation). See the related example in Stargazer and also the
  Documentation page for more details.
* Added a new tool named 'sges'. The major difference between this tool and
  the existing 'sgep' tool (formerly called 'sge') is that the latter was
  designed for running per-project genotyping pipeline (i.e. single gene,
  multiple samples) while the former is designed for running per-sample
  genotyping pipeline (i.e. multiple genes, single sample). Here, the
  letter p in "sgep" means SGE pipeline for per-[p]roject while the
  letter s in "sges" means SGE pipeline for per-[s]ample. The "sgep" tool
  is more suited for projects involving a large volume of samples (e.g.
  cohort-based research study) while the "sges" tool is more suited for
  situations in which the return of results is expected to occur (e.g.
  individual PGx report). See the Documentation page for more details.

## v1.0.6

* Official release on September 20, 2019.
* Extended to CYP2A13, CYP2F1, CYP2J2, CYP2R1, CYP2S1, CYP2W1, CYP3A7,
  CYP3A43, CYP26A1, POR.
* We're very excited to introduce the Database of Pharmacogenomic
  Structural Variants (DPSV)! The objective of DPSV is to provide a
  comprehensive summary of PGx SVs detected by Stargazer from real NGS
  samples (see examples in the DPSV page).
* This version makes allele fraction profiles as well as copy number profiles.
* This version uses phased allele fractions to determine relative gene
  copy numbers in samples with CN>2 (e.g., CYP2D6\*1/\*2x2 vs. \*1x2/\*2).
* This version uses more advanced SV detection algorithms, including the
  use of copy number-stable region or CNSR.
* This version uses improved systems for handling input VCF files created
  from various tools/sources.
* This version uses expanded haplotype reference panels for increased
  phasing accuracy (+/- 100kb instead of 3kb).
* This version updates Beagle v5.0 -> v5.1 for increased accuracy and speed.
* This version uses expanded target regions for more accurate SV detection.
* This version produces logging messages that are more user-friendly.
* Some command line arguments have been changed. See the Documentation page.

## v1.0.5

* Official release on July 23, 2019.
* Extended to G6PD and NUDT15.
* Additional tools have been added to Stargazer and, because of this, the
  command line has been changed.
* This version uses Beagle v5.0 instead of Beagle v4.1 for phasing SNVs/indels.
* Stargazer now supports "VCF only" mode for both NGS data and SNP chip data.

## v1.0.4

* Official release on March 3, 2019.
* Stargazer has been extended to call star alleles in 28 genes.
* Many useful optional arguments have been added.
* This version is described in Lee et al., 2019.

## v1.0.3

* Official release on July 9, 2018.
* Extended to CYP2A6/CYP2A7.

## v1.0.2

* Official release on June 14, 2018.
* To determine the duplicated star allele in samples with three gene copies
  or more (e.g., CYP2D6\*1x2/\*4 vs. \*1/\*4x2), Stargazer computes allele
  fractions from sequence reads that carry the corresponding variant.
  Previous versions of Stargazer test if the observed allele fraction
  from a sample with three gene copies or more is greater than the mean
  of allele fractions from all samples within the same sequencing project
  that are heterozygous for the variant of interest and do not have any
  structural variation. This version instead uses an optimal decision
  boundary found with Bayesian updating for two main reasons. First, the
  empirical mean is not always obtainable (i.e., there is only one sample
  with the variant) or the mean value might not be accurate if not many
  samples have this variant. Second, the approach allows utilization of an
  informative prior that says allele fractions should be centered at 0.5
  if heterozygous samples without structural variation are used.

## v1.0.1

* Official release on April 11, 2018.
* For detection of structural variation, this version no longer filters
  out loci based on the variance in read depth across the samples. Instead,
  it filters out pre-selected regions that have been empirically shown to
  produce high noise (e.g., regions in which reads are mapping multiply).
* In order to call structural variants, this version fits every pairwise
  combination of known sequence structures (one for each chromosome) against
  the sample's observed copy number profile and then selects the combination
  that produces the least deviance. This combinatorial testing is used in
  Stargazer_v1.0.0 as well but only for samples with more than one
  structural variation (abnormal structure for the first chromosome and
  abnormal structure for the second chromosome). Basically, this version
  generalizes the combinatorial testing to be applied to even samples
  without any structural variation (normal structure and normal structure)
  and samples with only one structural variation (normal structure and
  abnormal structure). As a result, the copy number plot now displays the
  two best sequence structures for the sample in addition to the sample's
  original copy number.

## v1.0.0

* Official release on March 11, 2018.
* This version is described in Lee et al., 2018.

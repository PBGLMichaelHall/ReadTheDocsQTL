===============
QTL_BSA-Sorghum
===============

:Author: Michael Hall
:Date:   4/13/2022

QTLSorghum
==========

QTLseqr is an R package for QTL mapping using NGS Bulk Segregant
Analysis.

QTLseqr is still under development and is offered with out any
guarantee.

**For more detailed instructions please read the vignette**\ `here <https://github.com/bmansfeld/QTLseqr/raw/master/vignettes/QTLseqr.pdf>`__
---------------------------------------------------------------------------------------------------------------------------------------------

For updates read the `NEWS.md <https://github.com/bmansfeld/QTLseqr/blob/master/NEWS.md>`__
-------------------------------------------------------------------------------------------

Installation
============

.. raw:: html

   <!-- You can install and update QTLseqr by using our [drat](http://dirk.eddelbuettel.com/code/drat.html) repository hosted on our github page: -->

.. raw:: html

   <!-- ```{r drat-install, eval = FALSE} -->

.. raw:: html

   <!-- install.packages("QTLseqr", repos = "http://bmansfeld.github.io/drat") -->

.. raw:: html

   <!-- ``` -->

.. raw:: html

   <!-- OR You can install QTLseqr from github with: -->

You can install QTLseqr from github with:

.. code:: r

   # install devtools first to download packages from github
   install.packages("devtools")

   # use devtools to install QTLseqr
   devtools::install_github("PBGLMichaelHall/QTLseqr")

**Note:** Apart from regular package dependencies, there are some
Bioconductor tools that we use as well, as such you will be prompted to
install support for Bioconductor, if you haven’t already. QTLseqr makes
use of C++ to make some tasks significantly faster (like counting SNPs).
Because of this, in order to install QTLseqr from github you will be
required to install some compiling tools (Rtools and Xcode, for Windows
and Mac, respectively).

**If you use QTLseqr in published research, please cite:**

   Mansfeld B.N. and Grumet R, QTLseqr: An R package for bulk segregant
   analysis with next-generation sequencing *The Plant Genome*
   `doi:10.3835/plantgenome2018.01.0006 <https://dl.sciencesocieties.org/publications/tpg/abstracts/11/2/180006>`__

We also recommend citing the paper for the corresponding method you work
with.

QTL-seq method:

   Takagi, H., Abe, A., Yoshida, K., Kosugi, S., Natsume, S., Mitsuoka,
   C., Uemura, A., Utsushi, H., Tamiru, M., Takuno, S., Innan, H., Cano,
   L. M., Kamoun, S. and Terauchi, R. (2013), QTL-seq: rapid mapping of
   quantitative trait loci in rice by whole genome resequencing of DNA
   from two bulked populations. *Plant J*, 74: 174–183.
   `doi:10.1111/tpj.12105 <https://onlinelibrary.wiley.com/doi/full/10.1111/tpj.12105>`__

G prime method:

   Magwene PM, Willis JH, Kelly JK (2011) The Statistics of Bulk
   Segregant Analysis Using Next Generation Sequencing. *PLOS
   Computational Biology* 7(11): e1002255.
   `doi.org/10.1371/journal.pcbi.1002255 <http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002255>`__

Abstract
--------

Next Generation Sequencing Bulk Segregant Analysis (NGS-BSA) is
efficient in detecting quantitative trait loci (QTL). Despite the
popularity of NGS-BSA and the R statistical platform, no R packages are
currently available for NGS-BSA. We present QTLseqr, an R package for
NGS-BSA that identifies QTL using two statistical approaches: QTL-seq
and G’. These approaches use a simulation method and a tricube smoothed
G statistic, respectively, to identify and assess statistical
significance of QTL. QTLseqr, can import and filter SNP data, calculate
SNP distributions, relative allele frequencies, G’ values, and
log10(p-values), enabling identification and plotting of QTL.

Examples:
=========

Load/install libraries
======================

\``\` r {r libraries} # install.packages(“tinytex”) #
install.packages(“vcfR”) # install.packages(“tidyr”) #
install.packages(“ggplot2”)
devtools::install_github(“PBGLMichaelHall/QTLseqr”,force = TRUE)
library(QTLseqr) library(tinytex) library(vcfR) library(tidyr)
library(ggplot2)

::

   # Set the Working Directory to where VCF file is stored in file system

   ``` r 
   setwd("/home/michael/Desktop/QTLseqr/extdata")

Vcf file must only contain bialleleic variants. (filter upstream, e.g., with bcftools view -m2 -M2), also the QTLseqR functions will only take SNPS, ie, length of REF and ALT== 1
==================================================================================================================================================================================

.. code:: r

   vcf <- read.vcfR(file = "freebayes_D2.filtered.vcf")

.. figure:: vertopal_3a51f92c0e3d4ee6abdbf942c924322c/69a7aa01ca2cde4201af39a38111309465617484.png
   :alt: Screenshot from 2022-03-30 15-12-19

   Screenshot from 2022-03-30 15-12-19

.. code:: r


   #Convert to tidy data frame
   VCF_TIDY <- vcfR2tidy(vcf)

.. figure:: vertopal_3a51f92c0e3d4ee6abdbf942c924322c/917999698fdd9e11d6f58a39baf048f9c0352229.png
   :alt: Screenshot from 2022-03-30 15-18-37

   Screenshot from 2022-03-30 15-18-37

Call the Parser
===============

.. code:: r

   QTLParser_1_MH(vcf = VCF_TIDY, HighBulk = "D2_F2_tt",LowBulk = "D2_F2_TT", filename = Hall)

.. figure:: vertopal_3a51f92c0e3d4ee6abdbf942c924322c/30037cf81667cbec9fdee3a2b653e9be0d2c5ac2.png
   :alt: Screenshot from 2022-03-30 15-10-38

   Screenshot from 2022-03-30 15-10-38

Preview the CSV file
====================

.. figure:: vertopal_3a51f92c0e3d4ee6abdbf942c924322c/7dc8731e7f30813c3fb29eec408af5c55c9860a2.png
   :alt: mcsv

   mcsv

Invoke unique command to extract Sample names reverse comapatible to the VCF
============================================================================

.. code:: r

   unique(VCF_TIDY$gt$Indiv)

.. figure:: vertopal_3a51f92c0e3d4ee6abdbf942c924322c/2a3484a525f9512ccc36ad623e56294731dab167.png
   :alt: Screenshot from 2022-03-30 15-33-29

   Screenshot from 2022-03-30 15-33-29

.. code:: r

   #Set High bulk and Low bulk sample names and parser generated file name
   #The file name is generated from the QTLParser_1_MH function in line 119

   HighBulk <- "D2_F2_tt"
   LowBulk <- "D2_F2_TT"
   file <- "Hall.csv"

   #Choose which chromosomes/contigs will be included in the analysis,

   Chroms <- c("Chr01","Chr02","Chr03","Chr04","Chr05","Chr06","Chr07","Chr08","Chr09","Chr10")

   df <-
     importFromTable(
       file = file,
       highBulk = HighBulk,
       lowBulk = LowBulk,
       chromList = Chroms
     ) 

.. figure:: vertopal_3a51f92c0e3d4ee6abdbf942c924322c/cc0678ad3d686e5bd31c0cb034647f45bbeb676e.png
   :alt: Screenshot from 2022-03-30 15-38-14

   Screenshot from 2022-03-30 15-38-14

Inspect the head of the df object
=================================

.. figure:: vertopal_3a51f92c0e3d4ee6abdbf942c924322c/be41d91010309b9265fc6ab00ddf0af65653ede8.png
   :alt: Screenshot from 2022-03-30 15-41-45

   Screenshot from 2022-03-30 15-41-45

.. code:: r

   #plot histograms associated with filtering arguments such as mamximum and minumum Total Depths and reference Allele Frequency to determine cut off values 

   ggplot(data =df) + geom_histogram(aes(x = DP.LOW + DP.HIGH)) + xlim(0,400)

   ggsave(filename = "Depth_Histogram.png",plot=last_plot())

.. figure:: vertopal_3a51f92c0e3d4ee6abdbf942c924322c/222b5a91a7323706dade57cf1212c0a3f1878523.png
   :alt: hist34

   hist34

.. code:: r


   ggplot(data = df) + geom_histogram(aes(x = REF_FRQ))

   ggsave(filename = "Ref_Freq_Histogram.png",plot = last_plot())

.. figure:: vertopal_3a51f92c0e3d4ee6abdbf942c924322c/a7a2aa48d18b933b3dc9482543946e052faff3ee.png
   :alt: hist4

   hist4

\``\` r {r Filtering, warning = FALSE}

#Filter SNPs based on some criteria df_filt <- filterSNPs( SNPset = df,
refAlleleFreq = 0.20, minTotalDepth = 100, maxTotalDepth = 400,
minSampleDepth = 40, # minGQ = 0 )

::

   ![Screenshot from 2022-04-01 08-29-37](https://user-images.githubusercontent.com/93121277/161207708-a50c7061-94ce-4729-b764-417bdb488d32.png)


   ```r

   #Run G' analysis
   df_filt<-runGprimeAnalysis_MH(
     SNPset = df_filt,
     windowSize = 5000000,
     outlierFilter = "deltaSNP",
     filterThreshold = 0.1)

.. figure:: vertopal_3a51f92c0e3d4ee6abdbf942c924322c/4bd4d836e4dc548ac16e5ac7ec0f0ed15032729e.png
   :alt: Screenshot from 2022-04-01 08-31-03

   Screenshot from 2022-04-01 08-31-03

G’ Distribution Plot
====================

.. code:: r

   #The plot reveals a skewed G Prime statistic with a really small variance. Perhaps it is due to the small number of variants called.
   #In addition, Hampels outlier filter in the second argument, can also be changed to "deltaSNP"
   plotGprimeDist(SNPset = df_filt, outlierFilter = "Hampel")

.. figure:: vertopal_3a51f92c0e3d4ee6abdbf942c924322c/e6315b7cd1928f2fdaff703e4fa32c2b91bb2247.png
   :alt: Screenshot from 2022-04-01 08-53-51

   Screenshot from 2022-04-01 08-53-51

.. code:: r

   #We can see raw data before and after our filtering step
   plotGprimeDist_MH(SNPset = df_filt, outlierFilter = "deltaSNP",filterThreshold = 0.1)

.. figure:: vertopal_3a51f92c0e3d4ee6abdbf942c924322c/5f15599f35d24769fa2bea8d5181ef975a00554c.png
   :alt: Screenshot from 2022-04-01 08-59-04

   Screenshot from 2022-04-01 08-59-04

::





   ``` r {r QTLSEQ, warning = FALSE}

   #Run QTLseq analysis
   df_filt2 <- runQTLseqAnalysis_MH(
     SNPset = df_filt,
     windowSize = 5000000,
     popStruc = "F2",
     bulkSize = c(45, 38),
     replications = 10000,
     intervals = c(95, 99)
   )

.. figure:: vertopal_3a51f92c0e3d4ee6abdbf942c924322c/549438872bc377baf4c36a4e7c9009f0c9da8488.png
   :alt: Screenshot from 2022-04-01 08-39-28

   Screenshot from 2022-04-01 08-39-28

Plot G Statistic Distribution
=============================

.. code:: r

   hist(df_filt2$G,breaks = 950,xlim = c(0,10),xlab = "G Distribution",main = "Histogram of G Values")

.. figure:: vertopal_3a51f92c0e3d4ee6abdbf942c924322c/5cae3910801aa3c0cb631ff7e10cd02dd3773fe2.png
   :alt: Screenshot from 2022-04-01 08-44-38

   Screenshot from 2022-04-01 08-44-38

.. code:: r

   #Plot Snps as a function of chromosome and position values
   plotQTLStats(SNPset = df_filt2, var = "nSNPs")
   ggsave(filename = "nSNPs.png",plot = last_plot())

.. figure:: vertopal_3a51f92c0e3d4ee6abdbf942c924322c/195b48fa6761e56d448532daa6f3695fc85cfce6.png
   :alt: Screenshot from 2022-04-01 15-20-34

   Screenshot from 2022-04-01 15-20-34

.. code:: r

   #Using QTLStats funciton plot Gprime Statistic with False Discovery Rate Threhshold as a third argument boolean operator as TRUE. The q value is used as FDR threshold null value is 0.05%.
   plotQTLStats(SNPset = df_filt, var = "Gprime", plotThreshold = TRUE, q = 0.01)
   ggsave(filename = "GPrime.png",plot = last_plot())

.. figure:: vertopal_3a51f92c0e3d4ee6abdbf942c924322c/206370151dcb11c9cc1848cc4c270fa1555907cc.png
   :alt: Screenshot from 2022-04-01 15-21-27

   Screenshot from 2022-04-01 15-21-27

.. code:: r

   #Again using plotQTLStats change second argument varaible to deltaSNP and plot.
   plotQTLStats(SNPset = df_filt2, var = "deltaSNP", plotIntervals  = TRUE)
   ggsave(filename = "DeltaSNPInterval.png",plot = last_plot())

.. figure:: vertopal_3a51f92c0e3d4ee6abdbf942c924322c/fb5a0b6910e3bf180b25d721e26651d0f32b060c.png
   :alt: Screenshot from 2022-04-01 15-22-05

   Screenshot from 2022-04-01 15-22-05

.. code:: r

   #Finally with plotQTLStats plot negLog10Pval
   plotQTLStats(SNPset = df_filt2, var = "negLog10Pval",plotThreshold = TRUE,q=0.01)
   ggsave(filename = "negLog10Pval.png",plot = last_plot())

.. figure:: vertopal_3a51f92c0e3d4ee6abdbf942c924322c/44d0c3a9b17264c41f6bd39e7dbb9f0c3ebd7735.png
   :alt: Screenshot from 2022-04-01 15-22-41

   Screenshot from 2022-04-01 15-22-41

.. code:: r

   #Add subset argument to focus on particular chromosomes one, three, four, and six.
   #The reason is due to signficant QTL regions
   plotQTLStats(SNPset = df_filt2, var = "Gprime",plotThreshold = TRUE,q=0.01,subset = c("Chr01","Chr03","Chr04","Chr06"))

.. figure:: vertopal_3a51f92c0e3d4ee6abdbf942c924322c/c523e8d291fea8b1c15d9801fafc03a37c71d838.png
   :alt: Screenshot from 2022-04-01 15-24-01

   Screenshot from 2022-04-01 15-24-01

Use RMVP package to view SNPs on chromosomes/contigs
====================================================

.. figure:: vertopal_3a51f92c0e3d4ee6abdbf942c924322c/833c68789aa223179cf0553e99f9b8de4a80a44f.png
   :alt: Screenshot from 2022-04-13 09-07-31

   Screenshot from 2022-04-13 09-07-31

Export summary CSV
==================

.. code:: r

   QTLTable(SNPset = df_filt, alpha = 0.01, export = TRUE, fileName = "my_BSA_QTL.csv")

Preview the Summary QTL
=======================

.. figure:: vertopal_3a51f92c0e3d4ee6abdbf942c924322c/40960888b938144206a456d12f1fb3baf1938cab.png
   :alt: Screenshot from 2022-04-01 09-18-54

   Screenshot from 2022-04-01 09-18-54

.. code:: r

   #Use the function to plot allele frequencies per chromosome
   #Second argument size specifes size of scalar factor on nSNPs and if you have a relatively small SNP set .001 is a good startin point otherwise set to 1
   Obs_Allele_Freq(SNPSet = df_filt, size = .001)

.. figure:: vertopal_3a51f92c0e3d4ee6abdbf942c924322c/f081b550d85fca66d4a358e7e4739330100f912d.png
   :alt: Screenshot from 2022-04-01 12-38-33

   Screenshot from 2022-04-01 12-38-33

.. code:: r

   ##Use the function to investigate chromosomal region of interest
   Obs_Allele_Freq2(SNPSet = df_filt, ChromosomeValue = "Chr04", threshold = .90)

.. figure:: vertopal_3a51f92c0e3d4ee6abdbf942c924322c/869d5d9c05ceb5d29acbd02cdd81638a8e953938.png
   :alt: Screenshot from 2022-04-01 15-45-06

   Screenshot from 2022-04-01 15-45-06

.. figure:: vertopal_3a51f92c0e3d4ee6abdbf942c924322c/db050b186cfd6e202e36d62ebd5065182a6596a1.png
   :alt: Screenshot from 2022-04-01 15-44-36

   Screenshot from 2022-04-01 15-44-36

.. code:: r

   setwd("/home/michael/Desktop/QTLseqr/extdata")
   # Theory and Analytical Framework of Sampling from BSA
   par(mfrow=c(1,1))
   # Define Ranges of Success
   success <- 0:90
   # The Difference between realized and Expected Frequencies 
   # ns : Sample Size taken from Low Bulk
   # 2(ns)p1_star ~ Binomial(2(ns),p1)
   # p1 Expected Frequencies
   # Expected Frequencies:
   # E(n1) = E(n2) = E(n3) = E(n4) = C/2 = 110
   # We prefer for accuracy to have ns >> C >> 1
   plot(success, dbinom(success, size = 90, prob = .50), type = "h",main="Binomial Sampling from Diploid Orgainism from High Bulk",xlab="2(ns)(p1_STAR)",ylab="Density")

.. figure:: vertopal_3a51f92c0e3d4ee6abdbf942c924322c/66008bdbf28e215b4e9ef9211aa869180233a89e.png
   :alt: Screenshot from 2022-04-12 13-59-19

   Screenshot from 2022-04-12 13-59-19

.. code:: r


   # ns : Sample Size from High Bulk
   # 2(ns)p2_star ~ Binomial(2(ns),p2)
   # p2 Expected Frequencies
   success <- 0:76
   plot(success, dbinom(success, size = 76, prob = 0.5), type = "h",main="Binomial Sampling from Diploid Organism from Low Bulk",xlab="2(n2)(p2_STAR)",ylab="Density")

.. figure:: vertopal_3a51f92c0e3d4ee6abdbf942c924322c/84dc9baeefd7e73372e54ef2e5d7c135f24f6220.png
   :alt: Screenshot from 2022-04-12 14-00-01

   Screenshot from 2022-04-12 14-00-01

.. code:: r



   # Read in the csv file from High bulk tt
   tt<-read.table(file = "D2_F2_tt.csv",header = TRUE,sep = ",")
   # Calculate average Coverage per SNP site
   mean(tt$DP)
   # Find REalized frequencies
   p1_STAR <- sum(tt$AD_ALT.) / sum(tt$DP)

   # Read in the csv file from Low Bulk TT
   TT<-read.table(file ="D2_F2_TT.csv",header = TRUE,sep=",")
   # Calculate average Coverage per SNP sit
   mean(TT$DP)
   # Find Realized frequencies
   p2_STAR <- sum(TT$AD_ALT.) / sum(TT$DP)
   # Take the average of the Averages
   C <-(mean(tt$DP)+mean(TT$DP))/2
   C<-round(C,0)
   # Find realized frequencies

   par(mfrow=c(1,1))
   #Define Ranges of Success (Allele Frequencies High and Low)
   success <- 0:100
   #n1|p1_star ~ Poisson(lambda)
   plot(success, dpois(success, lambda = C*(1-p1_STAR)), type = 'h',main="n1|p1_STAR ~ Poisson(C[1-p1_STAR])",xlab="n1|(n3/n1+n3)",ylab="Prob")

.. figure:: vertopal_3a51f92c0e3d4ee6abdbf942c924322c/7a5066025f47b1d9838a7b85b8ada85edd76fbb1.png
   :alt: Screenshot from 2022-04-12 14-00-53

   Screenshot from 2022-04-12 14-00-53

.. code:: r

   hist(TT$AD_REF., probability = TRUE,main="Histogram of Actually Realized n1 Values",xlab="n1")

.. figure:: vertopal_3a51f92c0e3d4ee6abdbf942c924322c/3428e831a768ea01144cb350cf10d6881440f663.png
   :alt: Screenshot from 2022-04-12 14-01-21

   Screenshot from 2022-04-12 14-01-21

.. code:: r

   #n2|p2_star ~ Poisson(lambda)
   plot(success, dpois(success, lambda = C*(1-p2_STAR)), type='h', main="n2|p2_STAR ~ Poisson(C[[1-p2_STAR])",xlab="n2|(n4/n2+n4)",ylab="Prob")

.. figure:: vertopal_3a51f92c0e3d4ee6abdbf942c924322c/bc0f7e6b430ffd2b5130be3b8c9241b1177b87f5.png
   :alt: Screenshot from 2022-04-12 14-01-49

   Screenshot from 2022-04-12 14-01-49

.. code:: r

   hist(tt$AD_REF., probability = TRUE, main = "Histogram of Actually Realized n2 Values",xlab="n2")

.. figure:: vertopal_3a51f92c0e3d4ee6abdbf942c924322c/9a11484175e6b38592b7f08857ba20c1a42cf063.png
   :alt: Screenshot from 2022-04-12 14-02-17

   Screenshot from 2022-04-12 14-02-17

.. code:: r

   #n3|p1_star ~ Poisson(lambda)
   plot(success, dpois(success, lambda = C*p1_STAR),type='h',main="n3|p1_STAR ~ Poisson(C[1-p1_STAR])",xlab="n3|(n3/n1+n3)",ylab="Prob")

.. figure:: vertopal_3a51f92c0e3d4ee6abdbf942c924322c/9e252e823d00f6181f8e06d8d86d6e55f0763d0a.png
   :alt: Screenshot from 2022-04-12 14-02-47

   Screenshot from 2022-04-12 14-02-47

.. code:: r

   hist(TT$AD_ALT., probability = TRUE, main="Histogram of Acutally Realized n3 Values",xlab="n3")

.. figure:: vertopal_3a51f92c0e3d4ee6abdbf942c924322c/6b53af63a62a6f2fdbc0fca552dee81558e11a5f.png
   :alt: Screenshot from 2022-04-12 14-03-11

   Screenshot from 2022-04-12 14-03-11

.. code:: r

   #n4|p2_star ~ Poisson(lambda)
   plot(success, dpois(success, lambda = C*p2_STAR), type = 'h',main="n4|p2_STAR ~ Poisson(C[1-p2_STAR])",xlab="n4|n4/(n2+n4)",ylab="Prob")

.. figure:: vertopal_3a51f92c0e3d4ee6abdbf942c924322c/158e623f9111a0b361cd99b7718d95adc26333c8.png
   :alt: Screenshot from 2022-04-12 14-03-41

   Screenshot from 2022-04-12 14-03-41

.. code:: r

   hist(tt$AD_ALT., probability = TRUE, main="Histogram of Acutally Realized n4 Values",xlab="n4")

.. figure:: vertopal_3a51f92c0e3d4ee6abdbf942c924322c/e9d156e83df7ec0d091aef884fc4f0406f524f1a.png
   :alt: Screenshot from 2022-04-12 14-04-02

   Screenshot from 2022-04-12 14-04-02

#Assuming average sequencing coverage expected values for n1,n2,n3,n4
C/2

p2 >> p1 QTL is present
=======================

However, ns >> C >> 1 is NOT TRUE
=================================

::






   ################################################################################################################################################

   # Rice QTL Analysis: High Bulk sample size of 430 Tolerant to cold environments and Low Bulk sample size of 385 Suseptilble to cold environments
   ```r


   #Set Working Directory
   setwd("/home/michael/Desktop/RiceCold2")

   #vcf file must only contain bialleleic variants. (filter upstream, e.g., with bcftools view -m2 -M2), also the QTL-Rice-Cold functions will only take SNPS, ie, length of REF and ALT== 1
   vcf <- read.vcfR(file = "wGQ-Filt-freebayes~bwa~IRGSP-1.0~both-segregant_bulks~filtered-default.vcf.gz")

.. figure:: vertopal_3a51f92c0e3d4ee6abdbf942c924322c/38f751c6ec44e430aab8f9ec2a9977d885465444.png
   :alt: Screenshot from 2022-04-01 09-27-24

   Screenshot from 2022-04-01 09-27-24

.. code:: r

   #Convert to tidy data frame
   VCF_TIDY <- vcfR2tidy(vcf)

.. figure:: vertopal_3a51f92c0e3d4ee6abdbf942c924322c/1cbe0915e4562eafbbdda84aa13d5c072c28734c.png
   :alt: Screenshot from 2022-04-01 09-37-35

   Screenshot from 2022-04-01 09-37-35

.. code:: r

   #Call the Parser
   QTLParser_1_MH(vcf = VCF_TIDY, HighBulk = "ET-pool-385",LowBulk = "ES-pool-430", filename = "Hall")

Standard Parser Output makes a list of Chromosome as named in VCF file
======================================================================

|Screenshot from 2022-04-01 09-41-00| # Also find unique sample names
reverse compatible with VCF file

.. code:: r

   unique(VCF_TIDY$gt$Indiv)

.. figure:: vertopal_3a51f92c0e3d4ee6abdbf942c924322c/04636c37909340f178f7765e1b3909804ed11034.png
   :alt: Screenshot from 2022-04-01 09-42-58

   Screenshot from 2022-04-01 09-42-58

.. code:: r

   #Set High bulk and Low bulk sample names and parser generated file name

   HighBulk <- "ET-pool-385"
   LowBulk <- "ES-pool-430"
   file <- "Hall.csv"

.. _preview-the-csv-file-1:

Preview the CSV file
====================

.. figure:: vertopal_3a51f92c0e3d4ee6abdbf942c924322c/ebb97f4dce8489b7644013df512e95b31d3512f1.png
   :alt: Screenshot from 2022-04-01 10-33-41

   Screenshot from 2022-04-01 10-33-41

.. code:: r

   #Input chromosomes values which will be included in the analysis,
   Chroms <- c("NC_029256.1","NC_029257.1","NC_029258.1","NC_029259.1","NC_029260.1","NC_029261.1","NC_029262.1","NC_029263.1","NC_029264.1","NC_029265.1","NC_029266.1","NC_029267.1")



   df <-
     importFromTable(
       file = file,
       highBulk = HighBulk,
       lowBulk = LowBulk,
       chromList = Chroms
     ) 

.. figure:: vertopal_3a51f92c0e3d4ee6abdbf942c924322c/ca1a5941413908c25760e96d76ee91c6a7cb7350.png
   :alt: Screenshot from 2022-04-01 10-35-10

   Screenshot from 2022-04-01 10-35-10

Plot histograms associated with filtering arguments to determine if cut off values are appropriate
==================================================================================================

.. code:: r


   ggplot(data =df) +geom_histogram(aes(x = DP.LOW + DP.HIGH)) + xlim(0,400)
   ggsave(filename = "Depth_Histogram.png",plot=last_plot())

.. figure:: vertopal_3a51f92c0e3d4ee6abdbf942c924322c/5d224fec26805da9230bd7f80e3ac999ca0b63eb.png
   :alt: dplowhigh

   dplowhigh

.. code:: r


   ggplot(data = df) +geom_histogram(aes(x = REF_FRQ))
   ggsave(filename = "Ref_Freq_Histogram.png",plot = last_plot())

.. figure:: vertopal_3a51f92c0e3d4ee6abdbf942c924322c/8309baa6f1ec32be1144a10173c2e0a7086691a8.png
   :alt: reffreq

   reffreq

.. code:: r



   #Filter SNPs based on some criteria
   df_filt <-
     filterSNPs(
       SNPset = df,
       refAlleleFreq = 0.20,
       minTotalDepth = 100,
       maxTotalDepth = 400,
       minSampleDepth = 40,
       #    minGQ = 0
     )
     

.. figure:: vertopal_3a51f92c0e3d4ee6abdbf942c924322c/cac0927fa9421d897a9a79dd8c99931aeb74344a.png
   :alt: Screenshot from 2022-04-01 09-54-07

   Screenshot from 2022-04-01 09-54-07

.. code:: r

   #Run G' analysis
   df_filt<-runGprimeAnalysis_MH(
     SNPset = df_filt,
     windowSize = 1e6,
     outlierFilter = "deltaSNP",
     filterThreshold = 0.1)

.. figure:: vertopal_3a51f92c0e3d4ee6abdbf942c924322c/7123d7e92d073f85ae9fccdaa6c2f613743b70d2.png
   :alt: Screenshot from 2022-04-01 09-54-50

   Screenshot from 2022-04-01 09-54-50

.. code:: r



   #Run QTLseq analysis
   df_filt2 <- runQTLseqAnalysis_MH(
     SNPset = df_filt,
     windowSize = 1e6,
     popStruc = "F2",
     bulkSize = c(430, 385),
     replications = 10000,
     intervals = c(95, 99)
   )

.. figure:: vertopal_3a51f92c0e3d4ee6abdbf942c924322c/9f88cc692a8fb4aa986e1039643eeaa1eea7ca3e.png
   :alt: Screenshot from 2022-04-01 09-56-09

   Screenshot from 2022-04-01 09-56-09

.. _plot-g-statistic-distribution-1:

Plot G Statistic Distribution
=============================

.. code:: r

   hist(df_filt2$G,breaks = 950,xlim = c(0,10),xlab = "G Distribution",main = "Histogram of G Values")

.. figure:: vertopal_3a51f92c0e3d4ee6abdbf942c924322c/eaa2c527519040e11c8c67e4c2ae330a1f54fad1.png
   :alt: gstat

   gstat

.. code:: r


   # G' Distribution Plot
   plotGprimeDist_MH(SNPset = df_filt2, outlierFilter = "Hampel")
   ggsave(filename = "Hampel_GPrime.png",plot = last_plot())

.. figure:: vertopal_3a51f92c0e3d4ee6abdbf942c924322c/dcf618107bbe6a42e15745d12e6a7856ffb9fe29.png
   :alt: gprime

   gprime

.. code:: r



   plotGprimeDist_MH(SNPset = df_filt2, outlierFilter = "deltaSNP",filterThreshold = 0.1)
   ggsave(filename = "DeltaSNP.png",plot = last_plot())

.. figure:: vertopal_3a51f92c0e3d4ee6abdbf942c924322c/bf83cea9ea886ebcf0818f42a11fb99b25d23438.png
   :alt: deltaSNP

   deltaSNP

.. code:: r



   #make the Plot
   plotQTLStats(SNPset = df_filt2, var = "nSNPs")
   ggsave(filename = "nSNPs.png",plot = last_plot())

.. figure:: vertopal_3a51f92c0e3d4ee6abdbf942c924322c/c4bcbb10de269549361d77c2f8a4e611d2787535.png
   :alt: Screenshot from 2022-04-01 10-36-37

   Screenshot from 2022-04-01 10-36-37

::



   ```r
   plotQTLStats(SNPset = df_filt, var = "Gprime", plotThreshold = TRUE, q = 0.01)
   ggsave(filename = "GPrime.png",plot = last_plot())

.. figure:: vertopal_3a51f92c0e3d4ee6abdbf942c924322c/de52a93ff150b69a0e3a7018a4cdfe2a0829a0bf.png
   :alt: Screenshot from 2022-04-01 10-37-50

   Screenshot from 2022-04-01 10-37-50

.. code:: r

   plotQTLStats(SNPset = df_filt2, var = "deltaSNP", plotIntervals  = TRUE)
   ggsave(filename = "DeltaSNPInterval.png",plot = last_plot())

.. figure:: vertopal_3a51f92c0e3d4ee6abdbf942c924322c/8e43311a8673a6b6c7d6d4623074c2c31b40dd65.png
   :alt: Screenshot from 2022-04-01 10-39-24

   Screenshot from 2022-04-01 10-39-24

.. code:: r

   plotQTLStats(SNPset = df_filt2, var = "negLog10Pval",plotThreshold = TRUE,q=0.01,subset = c("NC_029256.1","NC_029257.1","NC_029263.1","NC_029265.1"))
   ggsave(filename = "negLog10Pval.png",plot = last_plot())

.. figure:: vertopal_3a51f92c0e3d4ee6abdbf942c924322c/d66b6486773f1b82de7056cda9f0cd0c8da7910b.png
   :alt: Screenshot from 2022-04-01 10-05-35

   Screenshot from 2022-04-01 10-05-35

.. code:: r


   plotQTLStats(SNPset = df_filt2, var = "Gprime",plotThreshold = TRUE,q=0.01,subset = c("NC_029256.1","NC_029257.1","NC_029263.1","NC_029265.1"))

.. figure:: vertopal_3a51f92c0e3d4ee6abdbf942c924322c/1d13bf78fba89721d0288c9d309a34596ab8ed0d.png
   :alt: Screenshot from 2022-04-01 10-06-19

   Screenshot from 2022-04-01 10-06-19

Lets take a look at SNPs per chromosome using rMVP Package
==========================================================

.. code:: r

   sample<-"RiceColdTolerance"
   pathtosample <- "/home/michael/Desktop/RiceCold2/wGQ-Filt-freebayes~bwa~IRGSP-1.0~both-segregant_bulks~filtered-default.vcf"
   out<- paste0("mvp.",sample,".vcf")
   memo<-paste0(sample)
   dffile<-paste0("mvp.",sample,".vcf.geno.map")

   message("Making MVP data S1")
   MVP.Data(fileVCF=pathtosample,
            #filePhe="Phenotype.txt",
            fileKin=FALSE,
            filePC=FALSE,
            out=out
   )
   message("Reading MVP Data S1")
   df <- read.table(file = dffile, header=TRUE)
   message("Making SNP Density Plots")
   MVP.Report.Density(df[,c(1:3)], bin.size = 1000000, col = c("blue", "yellow", "red"), memo = memo, file.type = "jpg", dpi=300)

.. figure:: vertopal_3a51f92c0e3d4ee6abdbf942c924322c/4cef037fdad54f47151f5124e0673fe43107ef00.png
   :alt: Screenshot from 2022-04-12 14-34-01

   Screenshot from 2022-04-12 14-34-01

.. _export-summary-csv-1:

Export summary CSV
==================

.. code:: r


   getQTLTable(SNPset = df_filt2, alpha = 0.01, export = TRUE, fileName = "my_BSA_QTL.csv")

Preview the QTL Summary
=======================

.. figure:: vertopal_3a51f92c0e3d4ee6abdbf942c924322c/819707c1c9e934534655362fa263ab4a62d5604c.png
   :alt: Screenshot from 2022-04-01 09-58-43

   Screenshot from 2022-04-01 09-58-43

.. code:: r


   #Use the function to plot allele frequencies per chromosome
   Obs_Allele_Freq(SNPSet = df_filt)

Looks Dense
===========

|LowB| |HighB|

Filter Low Allelic Depth Frequencies
====================================

.. code:: r

   ##Use the function to investigate chromosomal region of interest
   Obs_Allele_Freq2(SNPSet = df_filt, ChromosomeValue = 8, threshold = .75)

Preview the plot with idenitfied SNP positions
==============================================

|LB3| |HB3|

Investigate SNP@POS 24525659
============================

|header| |values|

.. |Screenshot from 2022-04-01 09-41-00| image:: vertopal_3a51f92c0e3d4ee6abdbf942c924322c/58b25cc4e2eff7318f295e59578277ac460f1e33.png
.. |LowB| image:: vertopal_3a51f92c0e3d4ee6abdbf942c924322c/c6f682bd1feee9cbe40b7e6e4f8b7560d7e3f9b3.png
.. |HighB| image:: vertopal_3a51f92c0e3d4ee6abdbf942c924322c/82f46c812413eef22b05013da743d8c6bb8c11f0.png
.. |LB3| image:: vertopal_3a51f92c0e3d4ee6abdbf942c924322c/f6785e6c19f60f31c25ce8cdbde6796765141888.png
.. |HB3| image:: vertopal_3a51f92c0e3d4ee6abdbf942c924322c/730b4f313191eca1c17c24eaf067da7a06a34fba.png
.. |header| image:: vertopal_3a51f92c0e3d4ee6abdbf942c924322c/1b4581665ed4569b7405a11a78eef3caf4a70813.png
.. |values| image:: vertopal_3a51f92c0e3d4ee6abdbf942c924322c/2b62cd5e52b70a956d02300074e82c5d67e3cdc2.png

=======================
QTL_Rice_Cold_Tolerance
=======================

:Author: Michael Hall
:Date:   4/13/2022

QTL-Rice-Cold-Tolerance
=======================




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
   
Package Dependencies
--------------------


**Note:** Apart from regular package dependencies, there are some
Bioconductor tools that we use as well, as such you will be prompted to
install support for Bioconductor, if you haven’t already. QTLseqr makes
use of C++ to make some tasks significantly faster (like counting SNPs).
Because of this, in order to install QTLseqr from github you will be
required to install some compiling tools (Rtools and Xcode, for Windows
and Mac, respectively).

Citation
========

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
----------------------

.. code:: r 
   
   install.packages(“vcfR”) 
   install.packages(“tidyr”) 
   install.packages(“ggplot2”)
   devtools::install_github(“PBGLMichaelHall/QTLseqr”,force = TRUE)   
   library(QTLseqr) 
   library(vcfR) 
   library(tidyr)
   library(ggplot2)
   library(dplyr)

::

Set the Working Directory
-------------------------

.. code:: r 

   setwd("/home/michael/Desktop/RiceCold2")

Pre-Filtering Rules
===================

.. code:: r

   Vcf file can contain bialleleic variants before parsing, however, out of a principal investigators preference, the user can (filter upstream, e.g., with bcftools view -m2 -M2), also the QTLseqR functions will only call SNPS, so filter out **INDELS** with the following command line.


.. figure:: ../images/upstreamfilter.png
   :alt: 


The Lonely Parser
=================

   Calling my Parser **QTLParser_1_MH**
   This method requires 4 arguments, a **vcf**, **highBulk**, **lowBulk**, and **filename**.
   Proceeding this Call you must invoke **importFromTable** before Filtering.

.. code:: r 

    df <- QTLParser_1_MH(vcf = "freebayes_D2.filtered.vcf", highBulk = "D2_F2_tt", lowBulk = "D2_F2_TT", filename = "Hall.csv")

Import Data
===========

  **Method 1 (Biased due to parser configuration)**
  
  Calling **importFromTable** on Hall.csv file 
  This method requires 5 inputs to 5 arguments, **file**, **highBulk**, **lowBulk**, **chromList** and **sep**.
  


importFromTable
---------------

.. code:: r

    Chroms <- c("NC_029256.1","NC_029257.1","NC_029258.1","NC_029259.1","NC_029260.1","NC_029261.1","NC_029262.1","NC_029263.1","NC_029264.1","NC_029265.1","NC_029266.1","NC_029267.1")
   
   df <- importFromTable(file = "Hall.csv", highBulk = "D2_F2_tt", lowBulk = "D2_F2_TT", chromList = Chroms, sep = ",")

   **Method 2 (Most convienent)**

   Calling **importFromVCF**
   This method requires 5 arguments, a vcf **file**, **highBulk**, **lowBulk**, **chromList**, **filename**, and **filter.**
   **The filtering argument is a Boolean accepting only TRUE or FALSE. If TRUE then it filters out all SNPs that did not "PASS" in that INFO field.**
   **If it is FALSE then there is no filter applied at all.** 


importFromVCF
-------------

.. code:: r

   Chroms <- c("NC_029256.1","NC_029257.1","NC_029258.1","NC_029259.1","NC_029260.1","NC_029261.1","NC_029262.1","NC_029263.1","NC_029264.1","NC_029265.1","NC_029266.1","NC_029267.1")
   
   df <- importFromVCF(file = "freebayes_D2.filtered.vcf",highBulk = "D2_F2_tt",lowBulk = "D2_F2_TT",chromList = Chroms,filename = "Hall",filter = FALSE)
   
GATK
----

   Method 3 (Best in my opinion)

   Calling **importFromGATK**
   This method requires 4 arguments, **a vcf file**, **highBulk**, **lowBulk**, and **chromlist**.
   If you do not have the software on your machine, first visit this website.
   https://gatk.broadinstitute.org/hc/en-us/articles/360036194592-Getting-started-with-GATK4
   Go to section 4 and click the first from left to right **here** hyperlink
.. code:: r   
   
.. figure:: ../images/GATKDownload.png
   :alt: 
   
   **Download the lastest version by clicking on gatk-4.2.6.1.zip**
   **Extract the contents into your Downloads Folder**
   **What open source opeating system are you running?**
   
   .. code:: r   
   
.. figure:: ../images/WhichVersionUbuntu.png
   :alt: 
   
   
   
.. code:: r   
   
.. figure:: ../images/GATKRelease.png
   :alt: 
  
   **Navigate to the folder containg gatk executable python script**
   
.. code:: r

.. figure:: ../images/gatk.png
   
   Call **VariantsToTable** sub executable program with all appropriate flags
   
   
.. code:: r

.. figure:: ../images/gatkcommand.png

   This should produce a file called **Hall.table**

    Chroms <- c("NC_029256.1","NC_029257.1","NC_029258.1","NC_029259.1","NC_029260.1","NC_029261.1","NC_029262.1","NC_029263.1","NC_029264.1","NC_029265.1","NC_029266.1","NC_029267.1")
   
   df <- importFromGATK(file = "Hall.table", highBulk = "D2_F2_tt", lowBulk = "D2_F2_TT", chromlist = Chroms)

   **Method 1 is the most biased and therefore cuts out more SNPs than Methods 2 & 3 which produce nearly identical SNP sets.**


.. code:: r




Input Fields
============

.. code:: r

   #Set High bulk and Low bulk sample names and parser generated file name
   #The file name is generated from the QTLParser_1_MH function in line 119

   HighBulk <- "ET-pool-385"
   LowBulk <- "ES-pool-430"
   file <- "Hall.csv"

   #Choose which chromosomes/contigs will be included in the analysis,

   Chroms <- c("NC_029256.1","NC_029257.1","NC_029258.1","NC_029259.1","NC_029260.1","NC_029261.1","NC_029262.1","NC_029263.1","NC_029264.1","NC_029265.1","NC_029266.1","NC_029267.1")


importFromTable
===============
.. code:: r

   df <-
     importFromTable(
       file = file,
       highBulk = HighBulk,
       lowBulk = LowBulk,
       chromList = Chroms
     ) 

Histograms
----------

.. code:: r

   #plot histograms associated with filtering arguments such as mamximum and minumum Total Depths and reference Allele Frequency to determine cut off        values 
   ggplot(data =df) + geom_histogram(aes(x = DP.LOW + DP.HIGH)) + xlim(0,400)
   ggsave(filename = "Depth_Histogram.png",plot=last_plot())

.. figure:: ../images/65.png
   :alt: 

.. code:: r

   ggplot(data = df) + geom_histogram(aes(x = REF_FRQ))
   ggsave(filename = "Ref_Freq_Histogram.png",plot = last_plot())

.. figure:: ../images/66.png
   :alt: 

filterSNPs
==========

.. code:: r

   #Filter SNPs based on some criteria 
   df_filt <- filterSNPs( SNPset = df,
   refAlleleFreq = 0.20, minTotalDepth = 100, maxTotalDepth = 400,
   minSampleDepth = 40, 
   # minGQ = 0 )

.. figure:: ../images/67.png
   :alt: 



runGprimeAnalysis_MH
====================

.. code:: r

   #Run G' analysis
   df_filt<-runGprimeAnalysis_MH(
     SNPset = df_filt,
     windowSize = 1e6,
     outlierFilter = "deltaSNP",
     filterThreshold = 0.1)

.. figure:: ../images/68.png
   :alt: 

 

plotGprimeDist_MH
==================

.. code:: r

   #The plot reveals a skewed G Prime statistic with a really small variance. Perhaps it is due to the small number of variants called.
   #In addition, Hampels outlier filter in the second argument, can also be changed to "deltaSNP"
   
   plotGprimeDist(SNPset = df_filt, outlierFilter = "Hampel")

.. figure:: ../images/69.png
   :alt: 


.. code:: r

   #We can see raw data before and after our filtering step
   
   plotGprimeDist_MH(SNPset = df_filt, outlierFilter = "deltaSNP",filterThreshold = 0.1)

.. figure:: ../images/70.png
   :alt: 
   

runQTLseqAnalysis_MH
====================

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

.. figure:: ../images/71.png
   :alt: 



Plot G Statistic Distribution as a Histogram
--------------------------------------------

.. code:: r

   hist(df_filt2$G,breaks = 950,xlim = c(0,10),xlab = "G Distribution",main = "Histogram of G Values")

.. figure:: ../images/72.png
   :alt:

plotQTLStats
============



nSNPs
-----

.. code:: r

   #Plot Snps as a function of chromosome and position values
   
   plotQTLStats(SNPset = df_filt2, var = "nSNPs")
   ggsave(filename = "nSNPs.png",plot = last_plot())

.. figure:: ../images/73.png
   :alt: 

Gprime
------

.. code:: r

   #Using QTLStats funciton plot Gprime Statistic with False Discovery Rate Threhshold as a third argument boolean operator as TRUE. The q value is used as FDR threshold null value is 0.05%.
   
   plotQTLStats(SNPset = df_filt, var = "Gprime", plotThreshold = TRUE, q = 0.01)
   ggsave(filename = "GPrime.png",plot = last_plot())

.. figure:: ../images/74.png
   :alt: 


deltaSNP
--------

.. code:: r

   #Again using plotQTLStats change second argument varaible to deltaSNP and plot.
   
   plotQTLStats(SNPset = df_filt2, var = "deltaSNP", plotIntervals  = TRUE)
   ggsave(filename = "DeltaSNPInterval.png",plot = last_plot())

.. figure:: ../images/75.png
   :alt: 

negLog10Pval
------------
 
.. code:: r

   #Finally with plotQTLStats plot negLog10Pval
   
   plotQTLStats(SNPset = df_filt, var = "negLog10Pval",plotThreshold = TRUE,q=0.15)
   ggsave(filename = "negLog10Pval.png",plot = last_plot())

.. figure:: ../images/76.png
   :alt: 

   
Gprime Subset
-------------

.. code:: r

   #Add subset argument to focus on particular chromosomes one, three, four, and six.
   #The reason is due to signficant QTL regions
   plotQTLStats(SNPset = df_filt, var = "Gprime",plotThreshold = TRUE,q=0.05,subset = c("NC_029256.1","NC_029258.1","NC_029259.1","NC_029261.1"))

.. figure:: ../images/77.png
   :alt:



rMVP Package
============

SNP Densities
--------------

.. code:: r

   #install.packages("rMVP")
   library(rMVP)
   sample<-"Semi_Dwarfism_in_Sorghum"
   pathtosample <- "/home/michael/Desktop/QTLseqr/extdata/subset_freebayes_D2.filtered.vcf.gz"
   out<- paste0("mvp.",sample,".vcf")
   memo<-paste0(sample)
   dffile<-paste0("mvp.",sample,".vcf.geno.map")

   message("Making MVP data S1")
   MVP.Data(fileVCF=pathtosample,
         #filePhe="Phenotype.txt",
         fileKin=FALSE,
         filePC=FALSE,
         out=out)
         
   message("Reading MVP Data S1")
   df <- read.table(file = dffile, header=TRUE)
   message("Making SNP Density Plots")
   MVP.Report.Density(df[,c(1:3)], bin.size = 1000000, col = c("blue", "yellow", "red"), memo = memo, file.type = "jpg", dpi=300)


.. figure:: ../images/78.png
   :alt: 

 

Export summary CSV
==================

.. code:: r

   QTLTable(SNPset = df_filt, alpha = 0.01, export = TRUE, fileName = "my_BSA_QTL.csv")

Preview the Summary QTL
-----------------------

.. figure:: ../images/79.png
   :alt: 

Theory
======

Contigency Table
----------------

.. figure:: ../images/contingency.png
   :alt: 

Obs_Allel_Freq
--------------

.. code:: r

   #Use the function to plot allele frequencies per chromosome
   #Second argument size specifes size of scalar factor on nSNPs and if you have a relatively small SNP set .001 is a good startin point otherwise set to 1
   Obs_Allele_Freq(SNPSet = df_filt, size = 1)


   
Obs_Allele_Freq2
----------------

.. code:: r

   #Use the function to plot allele frequencies per chromosome
   #Second argument size specifes size of scalar factor on nSNPs and if you have a relatively small SNP set .001 is a good startin point otherwise set to 1
   ##Use the function to investigate chromosomal region of interest
   Obs_Allele_Freq2(SNPSet = df_filt, ChromosomeValue = "NC_029263.1", threshold = .85)

.. figure:: ../images/80.png
   :alt:
   
.. figure:: ../images/233.png
   :alt:



Total Coverage and Expected Allelic Frequencies
-----------------------------------------------

.. code:: r

   #Assuming average sequencing coverage (C) expected values for n1,n2,n3,n4
   E(n1) = E(n2) = E(n3) = E(n4) = C/2 = 35

   # Read in the csv file from High bulk tt
   tt<-read.table(file = "ET-pool-385.csv",header = TRUE,sep = ",")
   # Calculate average Coverage per SNP site
   mean(tt$DP)
   # Find REalized frequencies
   p1_STAR <- sum(tt$AD_ALT.) / sum(tt$DP)

   # Read in the csv file from Low Bulk TT
   TT<-read.table(file ="ES-pool-430.csv",header = TRUE,sep=",")
   # Calculate average Coverage per SNP sit
   mean(TT$DP)
   # Find Realized frequencies
   p2_STAR <- sum(TT$AD_ALT.) / sum(TT$DP)
   # Take the average of the Averages
   C <-(mean(tt$DP)+mean(TT$DP))/2
   C<-round(C,0)
   #Average Coverage
   70
   C/2 = 35
  
   p2 >> p1 QTL is present
   However, ns >> C >> 1 is TRUE 


   
Theory and Analytical Framework of Sampling from BSA
====================================================


Binomial Sampling
-----------------
   
Low Bulk
---------
   
.. code:: r

   setwd("/home/michael/Desktop/QTLseqr/extdata")
   # Theory and Analytical Framework of Sampling from BSA
   par(mfrow=c(1,1))
   # Define Ranges of Success
   # Sample Size from High Bulk sn = 385
   success <- 0:770
   # The Difference between realized and Expected Frequencies 
   # ns : Sample Size taken from Low Bulk
   # 2(ns)p1_star ~ Binomial(2(ns),p1)
   # p1 Expected Frequencies
   # Expected Frequencies:
   # E(n1) = E(n2) = E(n3) = E(n4) = C/2 = 110
   # We prefer for accuracy to have ns >> C >> 1
   plot(success, dbinom(success, size = 770, prob = .50), type = "h",main="Binomial Sampling from Diploid Orgainism from Low Bulk",xlab="2(ns)(   p1_STAR)",ylab="Density")

.. figure:: ../images/LB.png
   :alt: 

High Bulk
---------

.. code:: r


   # ns : Sample Size from High Bulk
   # 2(ns)p2_star ~ Binomial(2(ns),p2)
   # p2 Expected Frequencies
   success <- 0:860
   plot(success, dbinom(success, size = 860, prob = 0.5), type = "h",main="Binomial Sampling from Diploid Organism from High Bulk",xlab="2(n2)(p2_STAR)",ylab="Density")

.. figure:: ../images/HB.png
   :alt: 

  
Conditional Distribution of n1 given realized average frequency
---------------------------------------------------------------


.. code:: r



   par(mfrow=c(1,1))
   #Define Ranges of Success (Allele Frequencies High and Low)
   success <- 0:100
   #n1|p1_star ~ Poisson(lambda)
   plot(success, dpois(success, lambda = C*(1-p1_STAR)), type = 'h',main="n1|p1_STAR ~ Poisson(C[1-p1_STAR])",xlab="n1|(n3/n1+n3)",ylab="Prob")

.. figure:: ../images/n1Gp1.png
   :alt: 

 
Observed n1
-----------

.. code:: r

   # Filter outliers
   TT <- TT %>% filter(AD_REF. <= 500)

   hist(TT$AD_REF., probability = FALSE,main="Histogram of Actually Realized n1 Values",xlab="n1",breaks = "Sturges")



.. figure:: ../images/N1.png
   :alt: 

 
Conditional Distribution of n2 given realized average frequency
--------------------------------------------------------------- 

.. code:: r

   #n2|p2_star ~ Poisson(lambda)
   plot(success, dpois(success, lambda = C*(1-p2_STAR)), type='h', main="n2|p2_STAR ~ Poisson(C[[1-p2_STAR])",xlab="n2|(n4/n2+n4)",ylab="Prob")

.. figure:: ../images/85.png
   :alt: 

Observed n2
-----------

.. code:: r

   tt <- tt %>% filter(AD_REF. <= 500)
   hist(tt$AD_REF., probability = TRUE, main = "Histogram of Actually Realized n2 Values",xlab="n2")

.. figure:: ../images/86.png
   :alt: 

 
Conditional Distribution of n3 given realized average frequency
--------------------------------------------------------------- 

.. code:: r

   #n3|p1_star ~ Poisson(lambda)
   plot(success, dpois(success, lambda = C*p1_STAR),type='h',main="n3|p1_STAR ~ Poisson(C[1-p1_STAR])",xlab="n3|(n3/n1+n3)",ylab="Prob")

.. figure:: ../images/59.png
   :alt: 

Observed n3
-----------

.. code:: r


   TT <- TT %>% filter(AD_ALT. <= 300)
   hist(TT$AD_ALT., probability = TRUE, main="Histogram of Acutally Realized n3 Values",xlab="n3")


.. figure:: ../images/60.png
   :alt:

Conditional Distribution of n4 given realized average frequency
--------------------------------------------------------------- 

.. code:: r

   #n4|p2_star ~ Poisson(lambda)
   plot(success, dpois(success, lambda = C*p2_STAR), type = 'h',main="n4|p2_STAR ~ Poisson(C[1-p2_STAR])",xlab="n4|n4/(n2+n4)",ylab="Prob")

.. figure:: ../images/n4Gp2.png
   :alt: 


Observed n4
-----------

.. code:: r

   hist(tt$AD_ALT., probability = TRUE, main="Histogram of Acutally Realized n4 Values",xlab="n4")

.. figure:: ../images/62.png
   :alt: 







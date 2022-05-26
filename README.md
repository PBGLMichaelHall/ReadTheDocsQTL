---
title: "QTL_BSA-Sorghum"
author: "Michael Hall"
date: "4/13/2022"
output: html_document
---
# Semi Dwarfism in Sorghum Mutant Lines 
QTLseqr is an R package for QTL mapping using NGS Bulk Segregant
Analysis.

QTLseqr is still under development and is offered with out any
guarantee.

I have forked the repository to my github account and added a generic parser in R to do a Quantitative Trait Locus Bulk Segregant analysis. The hyperlink to Read the Docs is provided below. The reason I forked his repository is because the analysis was limited to GATK parsed vcf files. To make it more robust I added a generic R Parser that does a similar parsing. The parser function takes 4 arguments, a VCF Tidy Data frame which first was a VCF file then read by read.vcfR from vcfR package and finally converted to VCF Tidy Data Frame by vcfR2tidy function in vcfR package. The second and third arguments are the respective names given for High Bulk Sample and Low Bulk Sample which is reverse compatible with original vcf file. Lastly, the fourth and final argument is a user specified file name given to CSV input file used at the beginning of downstream analysis utilitzing importFromTable Function.


# Read the Docs
https://qtl-bsa.readthedocs.io/en/latest/Sorghum.html

# Why this tool is useful in the Sorghum case is because it takes 7861 Called Variants and reduces it to 92 total variants dispersed on Chromosomes 1,3,4,6,7, and 9 . You can also filter by allelic frequency differences however there is not really a valid statistic to determine significance. But, how much of the phenotypic variation is explained by these identified QTLs?

# I have added the following functions to the QTLseqr package:
- it provides cut off values for G-Prime Statistic and neglog
- Prints a csv file of ONLY G-prime Significant Variants
- In this case, Rice has called 32,771 significant variants with respect to G Prime Statistic

# plotQTLStats_MH



```r

QTLseqr::plotQTLStats_MH(SNPset = df_filt2, var = "Gprime", plotThreshold = TRUE, q = 0.01)

```

![hidden_Variables](https://user-images.githubusercontent.com/93121277/168254018-6869143b-1ecd-4029-8fd3-fa1c3771231b.png)

![print](https://user-images.githubusercontent.com/93121277/168260947-c35711c0-2c7a-4555-8976-20d666ad4afe.png)

# tricube_Smooth 
Uses Local Regression, Likelihood and Density Estimation and a Local Polynomial Model Term to filter out noise in the data.
**stats::predict(locfit::locfit(Stat~locfit::lp(POS, h = windowsize, deg= deg, nn=nn)),POS)**

Where Stat is G Statistic, however other variables can be used given they are in the data frame like frequencies etc. This allows the user to choose optimal parameter values to smooth the model basd on degree of polynomial deg or Nearest neighbor component of the smoothing parameter. Default value for KNN is 0.7 and degree of polynomial is zero.


# Obs_Allele_Freq
Plot allelic frequences from High and Low Bulk using nSNPs in window as size argument. Hoever, the size argument is actually a scalar quantity multiplied against nSNPs and so depending on the magnitude of the variants it helps to scale it down or up. Please start with a unit scale '1' until a satisfactory plot is made. 
i.e.

```r

QTLseqr::Obs_Allele_Freq(SNPSet = df_filt, size = .001)

```
![obs](https://user-images.githubusercontent.com/93121277/168235788-1e5dbbd7-ce25-4373-921c-5a6677058642.png)


# Obs_Allele_Freq2
Plot Independent Allelic Frequencies from High and Low Bulk with SNP position used as identifiers. Also returns a sorted data frame with Gprime statistic decending. The idea is to compare Allelic Frequencies from each bulk with Gprime Statistic for SNP called. The threhold argument is the magnitude of difference between Allelic frequencies inherited from the high parent in the High Bulk and Allelic Frequencies inherited from the high parent in the Low Bulk. It is curretly set at 66% due to the expectation that the High Bulk will have an extremely enriched allele content in a significant Quantitative Trait Locus where 100% of the Alleles are called to be homozygous recessive. While assuming heterzygous alleles have not been filtered out the expectation will tend towards 33% in mendellian fashion, and so the difference is 66% approximately. However, in this case heterzygous calls have been filtered and the trait is determined to be recessive due to ratio patterns in offspring from the F2 generation. So, we can expect to see the highbulk trending towards 100% and the low bulk trending towards 0%.

i.e.

```r

QTLseqr::Obs_Allele_Freq2(SNPSet = df_filt, ChromosomeValue = "Chr04", threshold = .66)

```
![obsfreq2](https://user-images.githubusercontent.com/93121277/168234078-7ed348c3-bcf1-41ba-86f1-efb7d125340a.png)


![chart](https://user-images.githubusercontent.com/93121277/168235561-f3e518d5-7d9b-494f-a8aa-96a4cf6bb00a.png)



# Correlation
This function takes a VCF file and a chromosome list and plots correlation related plots. 
Required libraries:

```r

**base::library(Hmisc)** 
**base::library(PerformanceAnalytics)** 
**base::library(corrplot)**

```

```r

table <- QTLseqr::Correlation(vcffile = "SNPS_ONLY.freebayes_D2.filtered.vcf.gz", chromlist =  c("Chr01","Chr02","Chr03","Chr04","Chr05","Chr06","Chr07","Chr08","Chr09","Chr10"), p1 = FALSE, p2 = FALSE, p3 = FALSE, p4 = FALSE, p5 = TRUE)


```


# Due to limitation on R Memory I recommend to use one plot at a time.
i.e. 

p1 = TRUE, p2 = FALSE,......etc.
For example, plot 4, p4 = TRUE gives this plot, distributions and correlations of certain info fields.

![correlation](https://user-images.githubusercontent.com/93121277/168233558-f5e15dd8-e962-4c81-b16d-b938490428a0.png)

Or if you like heatmaps do plot 5, p5 = TRUE.

![heatmap](https://user-images.githubusercontent.com/93121277/168233848-c7e67d6a-2812-4da5-933b-fcb9a1438d28.png)




# AlleleFreqSlidingWindow
i.e.

```r

QTLseqr::AlleleFreqSlidingWindow(vcf = "SNPS_ONLY.freebayes_D2.filtered.vcf.gz",chromList = c("Chr01","Chr02","Chr03","Chr04","Chr05","Chr06","Chr07","Chr08","Chr09","Chr10") , windowSize = 1000000, highBulk ="D2_F2_tt", lowBulk = "D2_F2_TT", filename = "Sorghum101", threshold = .66)
dev.off()

```
![plot](https://user-images.githubusercontent.com/93121277/168237187-1960c835-acee-4e12-bffd-1c96612681eb.png)



### **For more detailed instructions please read the vignette [here](https://github.com/bmansfeld/QTLseqr/raw/master/vignettes/QTLseqr.pdf)**

### For updates read the [NEWS.md](https://github.com/bmansfeld/QTLseqr/blob/master/NEWS.md)

# Installation

<!-- You can install and update QTLseqr by using our [drat](http://dirk.eddelbuettel.com/code/drat.html) repository hosted on our github page: -->

<!-- ```{r drat-install, eval = FALSE} -->

<!-- install.packages("QTLseqr", repos = "http://bmansfeld.github.io/drat") -->

<!-- ``` -->

<!-- OR You can install QTLseqr from github with: -->

You can install QTLseqr from github with:

``` r
# install devtools first to download packages from github
utils::install.packages("devtools")

# use devtools to install QTLseqr
devtools::install_github("PBGLMichaelHall/QTLseqr")
```


**Note:** Apart from regular package dependencies, there are some
Bioconductor tools that we use as well, as such you will be prompted to
install support for Bioconductor, if you haven’t already. QTLseqr makes
use of C++ to make some tasks significantly faster (like counting SNPs).
Because of this, in order to install QTLseqr from github you will be
required to install some compiling tools (Rtools and Xcode, for Windows
and Mac, respectively).

**If you use QTLseqr in published research, please cite:**

> Mansfeld B.N. and Grumet R, QTLseqr: An R package for bulk segregant
> analysis with next-generation sequencing *The Plant Genome*
> [doi:10.3835/plantgenome2018.01.0006](https://dl.sciencesocieties.org/publications/tpg/abstracts/11/2/180006)

We also recommend citing the paper for the corresponding method you work
with.

QTL-seq method:

> Takagi, H., Abe, A., Yoshida, K., Kosugi, S., Natsume, S., Mitsuoka,
> C., Uemura, A., Utsushi, H., Tamiru, M., Takuno, S., Innan, H., Cano,
> L. M., Kamoun, S. and Terauchi, R. (2013), QTL-seq: rapid mapping of
> quantitative trait loci in rice by whole genome resequencing of DNA
> from two bulked populations. *Plant J*, 74: 174–183.
> [doi:10.1111/tpj.12105](https://onlinelibrary.wiley.com/doi/full/10.1111/tpj.12105)

G prime method:

> Magwene PM, Willis JH, Kelly JK (2011) The Statistics of Bulk
> Segregant Analysis Using Next Generation Sequencing. *PLOS
> Computational Biology* 7(11): e1002255.
> [doi.org/10.1371/journal.pcbi.1002255](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002255)

## Abstract

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

# Examples:

# Load/install libraries


``` r {r libraries}
devtools::install_github("PBGLMichaelHall/QTLseqr",force = TRUE)
# utils::install.packages("vcfR")
# utils::install.packages("tidyr")
# utils::install.packages("ggplot2")
# utils::install.packages("dplyr")
# utils::install.packages("data.table")
base::library(QTLseqr)
base::library(data.table)
base::library(dplyr)
base::library(tidyr)
base::library(vcfR)
base::library(ggplot2)

```
# Set the Working Directory to where VCF file is stored in file system

``` r 
base::setwd("/home/michael/Desktop/QTLseqr/extdata")
```
# Vcf file must only contain bialleleic variants and Pure SNPs. (filter upstream, e.g., with bcftools view -v snps -m2 -M2), also the QTLseqR functions will only take SNPS, ie, length of REF and ALT== 1. In addition, importFromVCF has another filter boolean value of TRUE or FALSE. If TRUE it filters to include only variants that have PASS in the FILTER INFO field. 

```r
# Invoke importFromVCF function and produce a .CSV file

QTLseqr::importFromVCF(file = "freebayes_D2.filtered.vcf",highBulk = "D2_F2_tt",lowBulk = "D2_F2_TT",chromList = c("Chr01","Chr02","Chr03","Chr04","Chr05","Chr06","Chr07","Chr08","Chr09","Chr10"),filename = "Hall",filter=TRUE)

```
![importVCF](https://user-images.githubusercontent.com/93121277/165956202-7ea2d997-6506-4e66-a1da-4e0510736c7e.png)

# The header of the .CSV file reveals 7,763 variant entries and a total of 16 columns. This file is used in the next step when importFromTable function is invoked.

![header](https://user-images.githubusercontent.com/93121277/168533725-5f34b628-928d-480a-bfac-481b23918f27.png)



```r
#Set High bulk and Low bulk sample names and parser generated file name
#The file name is generated from the QTLParser_1_MH function in line 119

HighBulk <- "D2_F2_tt"
LowBulk <- "D2_F2_TT"
file <- "Hall.csv"

#Choose which chromosomes/contigs will be included in the analysis,

Chroms <- c("Chr01","Chr02","Chr03","Chr04","Chr05","Chr06","Chr07","Chr08","Chr09","Chr10")

df <-
  QTLseqr::importFromTable(
    file = file,
    highBulk = HighBulk,
    lowBulk = LowBulk,
    chromList = chromList
    sep = ","
  ) 

```

![removing](https://user-images.githubusercontent.com/93121277/165956524-43fb8aca-c6c1-4f35-ba99-863bbcee0e50.png)



```r
#plot histograms associated with filtering arguments such as mamximum and minumum Total Depths and reference Allele Frequency to determine cut off values 

ggplot2::ggplot(data =df) + geom_histogram(aes(x = DP.LOW + DP.HIGH)) + xlim(0,400)

ggplot2::ggsave(filename = "Depth_Histogram.png",plot=last_plot())

```

![hist34](https://user-images.githubusercontent.com/93121277/156784666-7abcc556-e0ef-4b4b-b981-c97074fccb0c.png)

```r

ggplot2::ggplot(data = df) + geom_histogram(aes(x = REF_FRQ))

ggplot2::ggsave(filename = "Ref_Freq_Histogram.png",plot = last_plot())

```

![hist4](https://user-images.githubusercontent.com/93121277/156784371-a96fcac0-fe21-43de-ad6c-568c69e3d25a.png)

```r

ggplot2::ggplot(data = df) + geom_histogram(aes(x = DP.LOW))

```

![dpLow](https://user-images.githubusercontent.com/93121277/170448999-976f0cd6-2fa7-4401-8a36-f5059565bfa1.png)

```r

ggplot2::ggplot(data = df) + geom_histogram(aes(x = DP.HIGH))

```

![dpHigh](https://user-images.githubusercontent.com/93121277/170449550-d4fbe07f-aff3-4b5e-b514-96afb13d9892.png)


```r

ggplot2::ggplot(data = df) + geom_histogram(aes(x = GQ.LOW))

```


![GQlow](https://user-images.githubusercontent.com/93121277/170449603-9ce8d582-f2b6-4550-8a86-7fed616fc962.png)


```r

ggplot2::ggplot(data = df) + geom_histogram(aes(x = GQ.HIGH))

```

![GQHigh](https://user-images.githubusercontent.com/93121277/170449502-f0833ddc-85cb-465c-928b-5a6013ddeee8.png)



``` r 

#Filter SNPs based on some criteria
df_filt <-
  QTLseqr::filterSNPs(
    SNPset = df,
    refAlleleFreq = 0.20,
    minTotalDepth = 100,
    maxTotalDepth = 400,
    minSampleDepth = 40,
    #    minGQ = 0
  )
```
![dffilt](https://user-images.githubusercontent.com/93121277/165956705-0faa10da-4fc2-4cfe-bc2c-b1bf8b51cdb2.png)



```r

#Run G' analysis
df_filt<-QTLseqr::runGprimeAnalysis(
  SNPset = df_filt,
  windowSize = 5000000,
  outlierFilter = "deltaSNP",
  filterThreshold = 0.1)
```
![Screenshot from 2022-04-01 08-31-03](https://user-images.githubusercontent.com/93121277/161207885-45119458-fa34-4259-80be-a8cd2f30fb3c.png)

# G' Distribution Plot
```r
#The plot reveals a skewed G Prime statistic with a really small variance. Perhaps it is due to the small number of variants called.
#In addition, Hampels outlier filter in the second argument, can also be changed to "deltaSNP"
QTLseqr::plotGprimeDist(SNPset = df_filt, outlierFilter = "Hampel")
```
![Screenshot from 2022-04-01 08-53-51](https://user-images.githubusercontent.com/93121277/161211346-780ed554-5fa6-4d22-bada-d446714c06aa.png)

```r
#We can see raw data before and after our filtering step
QTLseqr::plotGprimeDist(SNPset = df_filt, outlierFilter = "deltaSNP",filterThreshold = 0.1)
```
![Screenshot from 2022-04-01 08-59-04](https://user-images.githubusercontent.com/93121277/161212125-9daf552a-ddb1-4a50-8603-08adbf21a869.png)


```




``` r 

#Run QTLseq analysis
df_filt2 <- QTLseqr::runQTLseqAnalysis(
  SNPset = df_filt,
  windowSize = 5000000,
  popStruc = "F2",
  bulkSize = c(45, 38),
  replications = 10000,
  intervals = c(95, 99)
)

```
![Screenshot from 2022-04-01 08-39-28](https://user-images.githubusercontent.com/93121277/161209279-3bd76a3c-97d4-44f0-801d-c5c3e7ba1942.png)


# Plot G Statistic Distribution

```r
graphics::hist(df_filt2$G,breaks = 950,xlim = c(0,10),xlab = "G Distribution",main = "Histogram of G Values")
```

![Screenshot from 2022-04-01 08-44-38](https://user-images.githubusercontent.com/93121277/161209927-59cd31d5-5663-4d0b-9d7d-ae6ee64a180d.png)



```r
#Plot Snps as a function of chromosome and position values
QTLseqr::plotQTLStats(SNPset = df_filt2, var = "nSNPs")
ggplot2::ggsave(filename = "nSNPs.png",plot = last_plot())
```
![Screenshot from 2022-04-01 15-20-34](https://user-images.githubusercontent.com/93121277/161271810-6709f553-1361-4f4d-a6f4-3f11eca4a8c9.png)


```r
#Using QTLStats funciton plot Gprime Statistic with False Discovery Rate Threhshold as a third argument boolean operator as TRUE. The q value is used as FDR threshold null value is 0.05%.
QTLseqr::plotQTLStats(SNPset = df_filt, var = "Gprime", plotThreshold = TRUE, q = 0.01)
ggplot2::ggsave(filename = "GPrime.png",plot = last_plot())
```
![Screenshot from 2022-04-01 15-21-27](https://user-images.githubusercontent.com/93121277/161271924-e1e6319b-be94-440f-8898-76849c583edd.png)


```r
#Again using plotQTLStats change second argument varaible to deltaSNP and plot.
QTLseqr::plotQTLStats(SNPset = df_filt2, var = "deltaSNP", plotIntervals  = TRUE)
ggplot2::ggsave(filename = "DeltaSNPInterval.png",plot = last_plot())
```

![Screenshot from 2022-04-01 15-22-05](https://user-images.githubusercontent.com/93121277/161272031-d384dad1-c471-43bf-96d5-2e6cc05d2cfd.png)

```r
#Finally with plotQTLStats plot negLog10Pval
QTLseqr::plotQTLStats(SNPset = df_filt2, var = "negLog10Pval",plotThreshold = TRUE,q=0.01)
ggplot2::ggsave(filename = "negLog10Pval.png",plot = last_plot())
```
![Screenshot from 2022-04-01 15-22-41](https://user-images.githubusercontent.com/93121277/161272116-144fb1c9-78a1-4c7c-9a7c-96b1a1e0b664.png)


```r
#Add subset argument to focus on particular chromosomes one, three, four, and six.
#The reason is due to signficant QTL regions
QTLseqr::plotQTLStats(SNPset = df_filt2, var = "Gprime",plotThreshold = TRUE,q=0.01,subset = c("Chr01","Chr03","Chr04","Chr06"))

```
![Screenshot from 2022-04-01 15-24-01](https://user-images.githubusercontent.com/93121277/161272349-3c9bcf3e-553c-43a5-af5d-6b1e80e99658.png)

# Now to find some hidden variables

```r

QTLseqr::plotQTLStats_MH(SNPset = df_filt, var = "Gprime", plotThreshold = TRUE, q = 0.01)

```
# Even though there appears to be a signal on the whole chromosome 4, it actually only calls a total of 92 G-Prime Significant Variants.

![1](https://user-images.githubusercontent.com/93121277/168262569-345e5f13-2078-41b6-b5f4-2bb3c208db8e.png)
![2](https://user-images.githubusercontent.com/93121277/168262589-3f42e166-de1b-4f7c-8891-83096b08ee23.png)


# Use RMVP package to view SNPs on chromosomes/contigs

```r

utils::install.packages("rMVP")
base::library(rMVP)
sample<-"Semi_Dwarfism_in_Sorghum"
pathtosample <- "/home/michael/Desktop/QTLseqr/extdata/subset_freebayes_D2.filtered.vcf.gz"
out<- base::paste0("mvp.",sample,".vcf")
memo<-base::paste0(sample)
dffile<-paste0("mvp.",sample,".vcf.geno.map")

base::message("Making MVP data S1")
rMVP::MVP.Data(fileVCF=pathtosample,
      #filePhe="Phenotype.txt",
      fileKin=FALSE,
      filePC=FALSE,
      out=out)

base::message("Reading MVP Data S1")
df <- utils::read.table(file = dffile, header=TRUE)
base::message("Making SNP Density Plots")
rMVP::MVP.Report.Density(df[,c(1:3)], bin.size = 5000000, col = c("blue", "yellow", "red"), memo = memo, file.type = "jpg", dpi=300)

```

![Screenshot from 2022-04-13 09-07-31](https://user-images.githubusercontent.com/93121277/163119749-97a84175-a79b-4366-8d21-0d8a3aa7d1a3.png)



# Export summary CSV



```r
QTLseqr::getQTLTable(SNPset = df_filt, alpha = 0.01, export = TRUE, fileName = "my_BSA_QTL.csv")
```
# Preview the Summary QTL

![Screenshot from 2022-04-01 09-18-54](https://user-images.githubusercontent.com/93121277/161214771-b0ded302-846c-43c9-a8a4-a87b97423f63.png)

``` r
#Use the function to plot allele frequencies per chromosome
#Second argument size specifes size of scalar factor on nSNPs and if you have a relatively small SNP set .001 is a good startin point otherwise set to 1
QTLseqr::Obs_Allele_Freq(SNPSet = df_filt, size = .001)
```

![Screenshot from 2022-04-01 12-38-33](https://user-images.githubusercontent.com/93121277/161247921-2fab6d12-6b11-433a-af3b-a5c969aabc8a.png)
![slice](https://user-images.githubusercontent.com/93121277/170481262-38cd9ba6-3c6f-4a6c-b36e-474992033fc4.png)


# Filter the SNP set for significant P-values with respect to G Prime Statistics

```r

ggplot2::ggplot(data = df_filt, mapping = ggplot2::aes(x = pvalue)) + ggplot2::geom_histogram(bins = 100)
ggplot2::ggplot(data = df_filt, mapping= ggplot2::aes(x = pvalue)) + ggplot2::geom_density()
df_filt2 <- df_filt %>% dplyr::filter(pvalue < 0.05)
# After filtering only 96 Variants remain
base::plot(df_filt2$pvalue, pch = 20, col = "blue", xlab = "index", ylab = "pvalue")

```



![hist](https://user-images.githubusercontent.com/93121277/170462211-f2848b2e-40a0-412b-8b28-dc07f94df1aa.png)

**After pvalue filter only 96 variants remain**
![density](https://user-images.githubusercontent.com/93121277/170462221-aa24e5f0-4838-4738-8054-eb797ef452d0.png)
![pvalue](https://user-images.githubusercontent.com/93121277/170462231-2b7cba53-3ac7-4d3a-9d5d-353628be1196.png)


``` r
##Use the function to investigate chromosomal region of interest
QTLseqr::Obs_Allele_Freq2(SNPSet = df_filt, ChromosomeValue = c("Chr01","Chr02","Chr03","Chr04","Chr05", "Chr06", "Chr07", "Chr08", "Chr09", "Chr10", threshold = 0, pvalueThresh = 0.05)
```

![pvalueSIG](https://user-images.githubusercontent.com/93121277/170471126-2f31504b-e067-4992-8e6d-48b802a2387c.png)


![slice](https://user-images.githubusercontent.com/93121277/170484889-031979b7-3ff9-4845-842f-415f8639fbd0.png)




# Assuming average sequencing coverage expected values for n1,n2,n3,n4
# C/2


# p2 >> p1 QTL is present
# However, ns >> C >> 1 is NOT TRUE



################################################################################################################################################

# Rice QTL Analysis: High Bulk sample size of 430 Tolerant to cold environments and Low Bulk sample size of 385 Suseptilble to cold environments


# Lets take a look at SNPs per chromosome using rMVP Package
```r
sample<-"RiceColdTolerance"
pathtosample <- "/home/michael/Desktop/RiceCold2/wGQ-Filt-freebayes~bwa~IRGSP-1.0~both-segregant_bulks~filtered-default.vcf"
out<- base::paste0("mvp.",sample,".vcf")
memo<-base::paste0(sample)
dffile<-base::paste0("mvp.",sample,".vcf.geno.map")

base::message("Making MVP data S1")
rMVP::MVP.Data(fileVCF=pathtosample,
         #filePhe="Phenotype.txt",
         fileKin=FALSE,
         filePC=FALSE,
         out=out
)
base::message("Reading MVP Data S1")
df <- utils::read.table(file = dffile, header=TRUE)
base::message("Making SNP Density Plots")
rMVP::MVP.Report.Density(df[,c(1:3)], bin.size = 1000000, col = c("blue", "yellow", "red"), memo = memo, file.type = "jpg", dpi=300)

```
![Screenshot from 2022-04-12 14-34-01](https://user-images.githubusercontent.com/93121277/162963811-7fc215c6-ad65-4a01-b4b9-94aebed3e112.png)


```r


#Set Working Directory
base::setwd("/home/michael/Desktop/RiceCold2")

# vcf file must only contain bialleleic variants. (filter upstream, e.g., with bcftools view -m2 -M2), also the QTL-Rice-Cold functions will only take SNPS, ie, length of REF and ALT== 1

```

# Assuming Gatk is available and downloaded on your machine invoke this command to produce the input file to the downstream analysis.

```r

python gatk VariantsToTable --variant freebayes_D2.filtered.vcf --fields CHROM --fields POS --fields REF --fields ALT --genotyp-fields AD --genotype-fields DP --genotype-fields GQ --genotype-fields PL --output Hall.table

```
# Here is a preview of .Table file, it is a olain text file and so very difficult to comprehend due to it being so unorganized, etc.

![header2](https://user-images.githubusercontent.com/93121277/168534436-46f9361e-c967-468e-ba6b-beab16f48c6a.png)


```r
#Set High bulk and Low bulk sample names and parser generated file name

HighBulk <- "ET-pool-385"
LowBulk <- "ES-pool-430"
file <- "Hall.table"
```


```r
#Input chromosomes values which will be included in the analysis,
Chroms <- c("NC_029256.1","NC_029257.1","NC_029258.1","NC_029259.1","NC_029260.1","NC_029261.1","NC_029262.1","NC_029263.1","NC_029264.1","NC_029265.1","NC_029266.1","NC_029267.1")



df <-
  QTLseqr::importFromGATK(
    file = file,
    highBulk = HighBulk,
    lowBulk = LowBulk,
    chromList = Chroms
  ) 

```

![Screenshot from 2022-04-01 10-35-10](https://user-images.githubusercontent.com/93121277/161226938-8f50aa1b-27b4-4fe1-9b77-b06d14ec8b31.png)



# Plot histograms associated with filtering arguments to determine if cut off values are appropriate




```r

ggplot2::ggplot(data =df) +geom_histogram(aes(x = DP.LOW + DP.HIGH)) + xlim(0,400)
ggplot2::ggsave(filename = "Depth_Histogram.png",plot=last_plot())

```
![dplowhigh](https://user-images.githubusercontent.com/93121277/158780363-0939b60a-4a19-4104-b435-e352d715f5df.png)

```r

ggplot2::ggplot(data = df) +geom_histogram(aes(x = REF_FRQ))
ggplot2::ggsave(filename = "Ref_Freq_Histogram.png",plot = last_plot())

```

![reffreq](https://user-images.githubusercontent.com/93121277/158780481-24a6992e-e239-41d4-9868-3b6da831e757.png)



```r


#Filter SNPs based on some criteria
df_filt <-
  QTLseqr::filterSNPs(
    SNPset = df,
    refAlleleFreq = 0.20,
    minTotalDepth = 100,
    maxTotalDepth = 400,
    minSampleDepth = 40,
    #    minGQ = 0
  )
  
  ```
  
  ![Screenshot from 2022-04-01 09-54-07](https://user-images.githubusercontent.com/93121277/161220139-0c079197-99db-4292-949a-cd23d0e2e7c7.png)


```r
#Run G' analysis
df_filt<- QTLseqr::runGprimeAnalysis(
  SNPset = df_filt,
  windowSize = 1e6,
  outlierFilter = "deltaSNP",
  filterThreshold = 0.1)


```
![Screenshot from 2022-04-01 09-54-50](https://user-images.githubusercontent.com/93121277/161220373-00be3d7e-b67d-44ca-a1af-5c4670751f39.png)


```r


#Run QTLseq analysis
df_filt2 <- QTLseqr::runQTLseqAnalysis(
  SNPset = df_filt,
  windowSize = 1e6,
  popStruc = "F2",
  bulkSize = c(430, 385),
  replications = 10000,
  intervals = c(95, 99)
)

```
![Screenshot from 2022-04-01 09-56-09](https://user-images.githubusercontent.com/93121277/161220504-ba5f8e90-126f-4ee7-9a1d-06d5e09e6a94.png)

# Plot G Statistic Distribution
```r
graphics::hist(df_filt2$G,breaks = 950,xlim = c(0,10),xlab = "G Distribution",main = "Histogram of G Values")
```
![gstat](https://user-images.githubusercontent.com/93121277/158780626-0dd9efaa-8c2b-448e-8e22-ce94a1ff6fbf.png)


```r

# G' Distribution Plot
QTLseqr::plotGprimeDist(SNPset = df_filt2, outlierFilter = "Hampel")
ggplot2::ggsave(filename = "Hampel_GPrime.png",plot = last_plot())

```

![gprime](https://user-images.githubusercontent.com/93121277/158780745-ce684a8b-5267-42f2-aab1-9038b66490a5.png)


```r


QTLseqr::plotGprimeDist(SNPset = df_filt2, outlierFilter = "deltaSNP",filterThreshold = 0.1)
ggplot2::ggsave(filename = "DeltaSNP.png",plot = last_plot())

```

![deltaSNP](https://user-images.githubusercontent.com/93121277/158780846-58095997-5814-4440-8594-e2c560412eee.png)


```r


#make the Plot
QTLseqr::plotQTLStats(SNPset = df_filt2, var = "nSNPs")
ggplot2::ggsave(filename = "nSNPs.png",plot = last_plot())


```
![Screenshot from 2022-04-01 10-36-37](https://user-images.githubusercontent.com/93121277/161227187-9a3d3138-e996-4965-9659-6a879d0be70a.png)



```


```r
QTLseqr::plotQTLStats(SNPset = df_filt, var = "Gprime", plotThreshold = TRUE, q = 0.01)
ggplot2::ggsave(filename = "GPrime.png",plot = last_plot())

```

![Screenshot from 2022-04-01 10-37-50](https://user-images.githubusercontent.com/93121277/161227403-aa546c4a-c6f0-4c3c-8c25-485bae654739.png)


```r
QTLseqr::plotQTLStats(SNPset = df_filt2, var = "deltaSNP", plotIntervals  = TRUE)
ggplot2::ggsave(filename = "DeltaSNPInterval.png",plot = last_plot())

```
![Screenshot from 2022-04-01 10-39-24](https://user-images.githubusercontent.com/93121277/161227673-447d11f9-1245-4042-affd-56ca43c39799.png)





```r
QTLseqr::plotQTLStats(SNPset = df_filt2, var = "negLog10Pval",plotThreshold = TRUE,q=0.01,subset = c("NC_029256.1","NC_029257.1","NC_029263.1","NC_029265.1"))
ggplot2::ggsave(filename = "negLog10Pval.png",plot = last_plot())

```
![Screenshot from 2022-04-01 10-05-35](https://user-images.githubusercontent.com/93121277/161222325-69f4ad58-8fb8-4eb7-bc91-fc900be47416.png)


```r

QTLseqr::plotQTLStats(SNPset = df_filt2, var = "Gprime",plotThreshold = TRUE,q=0.01,subset = c("NC_029256.1","NC_029257.1","NC_029263.1","NC_029265.1"))

```
![Screenshot from 2022-04-01 10-06-19](https://user-images.githubusercontent.com/93121277/161222339-dbbafc12-0631-4d22-b477-5f038256ee6e.png)

# Lets take a look at SNPs per chromosome using rMVP Package
```r
sample<-"RiceColdTolerance"
pathtosample <- "/home/michael/Desktop/RiceCold2/wGQ-Filt-freebayes~bwa~IRGSP-1.0~both-segregant_bulks~filtered-default.vcf"
out<- base::paste0("mvp.",sample,".vcf")
memo<-base::paste0(sample)
dffile<-base::paste0("mvp.",sample,".vcf.geno.map")

base::message("Making MVP data S1")
rMVP::MVP.Data(fileVCF=pathtosample,
         #filePhe="Phenotype.txt",
         fileKin=FALSE,
         filePC=FALSE,
         out=out
)
base::message("Reading MVP Data S1")
df <- utils::read.table(file = dffile, header=TRUE)
base::message("Making SNP Density Plots")
rMVP::MVP.Report.Density(df[,c(1:3)], bin.size = 1000000, col = c("blue", "yellow", "red"), memo = memo, file.type = "jpg", dpi=300)

```
![Screenshot from 2022-04-12 14-34-01](https://user-images.githubusercontent.com/93121277/162963811-7fc215c6-ad65-4a01-b4b9-94aebed3e112.png)



# Export summary CSV
```r

QTLseqr::getQTLTable(SNPset = df_filt2, alpha = 0.01, export = TRUE, fileName = "my_BSA_QTL.csv")
```
# Preview the QTL Summary
![Screenshot from 2022-04-01 09-58-43](https://user-images.githubusercontent.com/93121277/161220947-979d5bbf-8438-4110-a950-a33224878a01.png)



# Filter Allelic Depth Frequencies, prefer 85% or greater.
```r

QTLseqr::Obs_Allele_Freq2(SNPSet = df_filt2, ChromosomeValue = "NC_029263.1", threshold = .85)


```

# Preview the plot with idenitfied SNP positions

![obs9](https://user-images.githubusercontent.com/93121277/168242146-7b92a7aa-c8b6-4d49-8dec-d70eba1e4f53.png)






![new](https://user-images.githubusercontent.com/93121277/170455852-49ff3597-182f-4877-846d-78571c201ea4.png)

![save](https://user-images.githubusercontent.com/93121277/170457107-f8f9b2ef-5396-483f-9a13-a69c86f1fe25.png)



![screen](https://user-images.githubusercontent.com/93121277/170471567-cbb92c0d-8d96-4777-afa8-d4184ffa2350.png)


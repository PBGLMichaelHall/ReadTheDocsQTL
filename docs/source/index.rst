.. QTL_BSA-Sorghum documentation master file, created by
   sphinx-quickstart on Wed Apr 13 10:30:05 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to QTL_BSA-Sorghum's documentation!
===========================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:


QTLSorghum
----------

QTLseqr is an R package for QTL mapping using NGS Bulk Segregant Analysis.

QTLseqr is still under development and is offered with out any guarantee.
For more detailed instructions please read the vignette here
For updates read the NEWS.md
Installation

You can install QTLseqr from github with:


# install devtools first to download packages from github
install.packages("devtools")

# use devtools to install QTLseqr
devtools::install_github("PBGLMichaelHall/QTLseqr")

Note: Apart from regular package dependencies, there are some Bioconductor tools that we use as well, as such you will be prompted to install support for Bioconductor, if you haven’t already. QTLseqr makes use of C++ to make some tasks significantly faster (like counting SNPs). Because of this, in order to install QTLseqr from github you will be required to install some compiling tools (Rtools and Xcode, for Windows and Mac, respectively).

If you use QTLseqr in published research, please cite:

    Mansfeld B.N. and Grumet R, QTLseqr: An R package for bulk segregant analysis with next-generation sequencing The Plant Genome doi:10.3835/plantgenome2018.01.0006

We also recommend citing the paper for the corresponding method you work with.

QTL-seq method:

    Takagi, H., Abe, A., Yoshida, K., Kosugi, S., Natsume, S., Mitsuoka, C., Uemura, A., Utsushi, H., Tamiru, M., Takuno, S., Innan, H., Cano, L. M., Kamoun, S. and Terauchi, R. (2013), QTL-seq: rapid mapping of quantitative trait loci in rice by whole genome resequencing of DNA from two bulked populations. Plant J, 74: 174–183. doi:10.1111/tpj.12105

G prime method:

    Magwene PM, Willis JH, Kelly JK (2011) The Statistics of Bulk Segregant Analysis Using Next Generation Sequencing. PLOS Computational Biology 7(11): e1002255. doi.org/10.1371/journal.pcbi.1002255

Abstract

Next Generation Sequencing Bulk Segregant Analysis (NGS-BSA) is efficient in detecting quantitative trait loci (QTL). Despite the popularity of NGS-BSA and the R statistical platform, no R packages are currently available for NGS-BSA. We present QTLseqr, an R package for NGS-BSA that identifies QTL using two statistical approaches: QTL-seq and G’. These approaches use a simulation method and a tricube smoothed G statistic, respectively, to identify and assess statistical significance of QTL. QTLseqr, can import and filter SNP data, calculate SNP distributions, relative allele frequencies, G’ values, and log10(p-values), enabling identification and plotting of QTL.

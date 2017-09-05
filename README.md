# rnaseqcomp: A benchmark for RNA-seq quantification pipelines

### Introduction
RNA sequencing (RNA-seq) has been utilized as the standard technology for
measuring the expression abundance of genes, transcripts, exons or splicing
junctions. Numerous quantification methods were proposed to quantify such
abundances with/without combination of RNA-seq read aligners. 
It is currently difficult to evaluate the performance of the best method, due
in part to the high costs of running assessment experiments as well as the
computational requirements of running these algorithms. *rnaseqcomp* package
provides a series of statistical summaries and data visualization techniques
to evaluate the performance of these metohods.

### Installation

*rnaseqcomp* is an R/bioconductor package, which can be installed with
source code documented
in [GitHub](https://github.com/tengmx/rnaseqcomp) or
through [Bioconductor](https://bioconductor.org/packages/rnaseqcomp).

Install *rnaseqcomp* through GitHub with following code.
```s
library(devtools)
install_github("tengmx/rnaseqcomp")
```

Or, install through Bioconductor as follows.
```s
source("https://bioconductor.org/biocLite.R")
biocLite("rnaseqcomp")
```

### Using *rnaseqcomp*

First, load the package into R.
```s
library(rnaseqcomp)
```

Details of how to use this package, please see the 
[vignette](https://github.com/tengmx/rnaseqcomp/blob/master/vignettes/rnaseqcomp.pdf).


### Help

You are very welcome to leave any questions/bug messages at
[GitHub issues](https://github.com/tengmx/rnaseqcomp/issues).

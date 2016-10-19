# rnaseqcomp: A benchmark for RNA-seq quantification pipelines

### Introduction
RNA sequencing (RNA-seq) has been utilized as the standard technology for
measuring the expression abundance of genes, transcripts, exons or splicing
junctions. Numerous quantification methods were proposed to quantify such
abundances with/without combination of RNA-seq read aligners. 
It is currently difficult to evaluate the performance of the best method, due
in part to the high costs of running assessment experiments as well as the
computational requirements of running these algorithms. **rnaseqcomp** package
provides a series of statistical summaries and data visualization techniques
to evaluate the performance of these metohods.

### Installation

R-package **rnaseqcomp** can be installed:
```s
library(devtools)
install_github("tengmx/rnaseqcomp")
```
After installation, the package can be loaded into R.

```s
library(rnaseqcomp)
```

### Using rnaseqcomp

Details of how to use this package, please see the 
[vignette](https://github.com/tengmx/rnaseqcomp/blob/master/vignettes/rnaseqcomp.pdf).
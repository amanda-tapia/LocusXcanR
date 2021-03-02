# LocusXcanR: An R Shiny application for integrating TWAS results with other 'omics data
LocusXcanR is an R package that creates an R Shiny application to assist with TWAS fine-mapping. It allows researchers to integrate information from multiple sources (GWAS and TWAS) and interactively visualize TWAS results in-context, one genomic locus at a time. LocusXcanR aids TWAS interpretation, can be extended to other ‘omics data, and highlights R Shiny’s effectiveness at presenting results in an approachable, interactive, and visual format.

For additional details and an application of LocusXcanR, please reference our manuscript currently available on bioRxiv: https://www.biorxiv.org/content/10.1101/2021.02.23.432444v1

Source code available in <code> app.R </code>

## Prerequisites
The LocusXcanR package depends on R (>= 4.0.3) and imports the following R packages:
<ul>
  <li>data.table</li>
  <li>dplyr</li>
  <li>DT</li>
  <li>ggplot2</li>
  <li>Gviz</li>
  <li>magrittr</li>
  <li>plotly</li>
  <li>RColorBrewer</li>
  <li>shiny</li>
  <li>stats</li>
  <li>tidyr</li>
  <li>visNetwork</li>
</ul>

## How to install the package from GitHub
First, install the devtools package. You can do this from CRAN. Invoke R and then type:

<code> install.packages("devtools") </code>

Load the devtools package:

<code>library(devtools)</code>

Install the package from GitHub using install_github("author/package") as follows:

<code>install_github("amanda-tapia/LocusXcanR")</code>

## Accessing LocusXcanR vignette
A vignette documenting the use of LocusXcanR is contained within the package. After installing the package, run:

<code>vignette("LocusXcanR")</code>

This should bring up the vignette in the help window of the R console. If this doesn't work, you may need to try the following:

<pre>
  <code> 
  devtools::install(build_vignettes = TRUE)
  vignette("LocusXcanR")
  </code>
</pre>

## Example scripts and data
See <code>R/example/example_script.R</code> to render an example R Shiny application.

Example datasets can be found in <code>inst/extdata</code>.

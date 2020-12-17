# LocusXcanR: An R Shiny application for integrating TWAS results with other 'omics data
This R Shiny application was developed to assist with TWAS fine-mapping. It allows researchers to integrate information from multiple sources (GWAS and TWAS) and interactively visualize TWAS results in-context, one genomic locus at a time.

Source code available in <code> app.R </code>

## Prerequisites
The LocusXcanR package depends on R (>= 3.5.0) and imports the following R packages:
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

<code>install_github("amanda-tapia/TWAS_RShiny")</code>

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

### Hosting your Shiny app locally
### Hosting your Shiny app on the Shiny server
### Hosting your Shiny app on another server


# R Shiny application for integrated TWAS results
This R Shiny application was developed to assist with TWAS fine-mapping. It allows researchers to integrate information from multiple sources (GWAS and TWAS) and interactively visualize TWAS results in-context, one genomic locus at a time.

Source code available in <code> app.R </code>

## Prerequisites
The following R packages are needed for plotting figures and reporting tables:
<ul>
  <li>shiny</li>
  <li>tidyverse</li>
  <li>DT</li>
  <li>plotly</li>
  <li>RColorBrewer</li>
  <li>data.table</li>
  <li>visNetwork</li>
</ul>

## app.R input data
<strong>TWAS results</strong> <br>
The main results (p-values required) from TWAS. Option to convert p-values to -log10(p-values) and option to include results (p-values required) from TWAS conditional analysis. Other variables in the data set are optional (can be included, excluded, or left blank).

TO-DO: add an example table here, and include variable definitions.

<strong>GWAS results</strong> <br>
The main results (p-values required) from GWAS. Option to convert p-values to -log10(p-values). Other variables in the data set are optional (can be included, excluded, or left blank).

TO-DO: add an example table here, and include variable definitions.

<strong>Gene expression prediction weights</strong><br>
The weights, or effect sizes, of model variants, trained in your gene expression reference panel. Required variables are gene identifier, variant identifier, model weight.

TO-DO: add an example table here, and include variable definitions.

<strong>Gene-gene correlations</strong>

TO-DO: add an example table here, and include variable definitions.

<strong>LD among GWAS variants</strong>

TO-DO: add an example table here, and include variable definitions.


## Setup and usage example
### Hosting your Shiny app locally
### Hosting your Shiny app on the Shiny server
### Hosting your Shiny app on another server


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

## How to install the package from GitHub
First, install the devtools package. You can do this from CRAN. Invoke R and then type:

<code> install.packages("devtools") </code>

Load the devtools package:

<code>library(devtools)</code>

Install the package from GitHub using install_github("author/package") as follows:

<code>install_github("amanda-tapia/TWAS_RShiny")</code>

### Hosting your Shiny app locally
### Hosting your Shiny app on the Shiny server
### Hosting your Shiny app on another server


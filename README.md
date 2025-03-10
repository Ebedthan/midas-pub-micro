## Analysis Code of MIDAS scientific article #1

### Overview

We used R to analyse microbiome data. To facilitate reproducibility, `renv` package was used to bundle everything needed to reproduce results.

### How to reproduce the results ?

1.  Download the data folder

The data folder was not added to GitHub (because of its size) but is publicly available on Zenodo at this link <https://doi.org/10.5281/zenodo.14989905>.

2.  Check you have the exact version of R used

We used **R version 4.4.3 (2025-02-28)** on a `x86_64-pc-linux-gnu` PC.

3.  Download the repository

Using the CLI you can use:

``` bash
git clone https://github.com/Ebedthan/midas-pub-micro.git
```

or directly download the zipped file from GitHub.

After downloading, decompress the folder if needed and navigate into the folder.

4.  Reproduce the results

First open `R` and then install `renv`:

``` r
install.packages("renv")
```

Then restore the project state as defined in the lockfile by using

``` r
renv::restore()
```

Then you can run the analysis in the `diversity_analysis.R` and get the results.

Enjoy!

Bioinformatics pipeline for Global Spore Sampling Project
================

<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->
<!-- badges: end -->

This pipeline is a pre-release version of the OptimOTU pipeline, which
is implemented in R using the
[{targets}](https://books.ropensci.org/targets/) workflow management
package. It is probably difficult to run it outside of Puhti, the
[CSC](https://www.csc.fi) SLURM computing cluster where it was developed
and executed for this project. Code is included primarily for reference
purposes to accompany the GSSP data release.

## How does the pipeline work?

`{targets}` is a workflow manager which takes a plan (an R list of
`tar_target` objects), constructs a directed acyclic graph (DAG)
representing the dependencies between the individual targets, and
executes the commands either one at a time or in parallel, while
ensuring that any prerequisites for each target are available before
starting to execute its command.

Unlike many other workflow managers (e.g., Snakemake, NextFlow, GNU
make), `{targets}` does not require that the prerequisites for each
target are explicitly listed. Instead, the names of targets can be used
in the commands as though they were normal R variables, and `{targets}`
uses static code analysis to determine the dependencies between
different targets.

The code to define the GSSP pipeline is contained in the `scripts`
directory. The convention which has been used is that files beginning
with “`0`” define functions and load settings which are required to run
the plan. Files beginning with “`1`” define the actual targets. The
logical flow of the pipeline can be approximately followed by looking at
the files in numerical order, although `{targets}` does not necessarily
execute the targets in the order they are defined.

A simple target definition looks like this:

``` r
#### name_of_target ####
# type:
#  - including explanation of each column
#  - if it is a data.frame or tibble
#
# further explanation if needed
tar_target( # or alternate versions like tar_fst_tbl and tar_file_fast
  name_of_target, # the name of the target :)
  f1(x, y, z) |> # command to generate the target
    f2(a, c),
  deployment = "main" # or other optional arguments to tar_target
)
```

`tar_target()` is the basic `{targets}` function which defines a target.
`tar_fst_tbl()` and `tar_fast_file()` which are also frequently used in
the pipeline, are functions from the
[`tarchetypes`](https://docs.ropensci.org/tarchetypes/) package which
also define targets, but the resulting target is stored differently on
disk.

`f1` and `f2` could be functions defined either in base R, an R package,
or in one of the `0xx`-level script files. `x`, `y`, `z`, `a`, and `c`
could be either explicit constants, constants defined in one of the
`0xx` level scripts, or other targets defined in the `1xx`-level
scripts. `|>` is the base R [pipe
operator](https://search.r-project.org/R/refmans/base/html/pipeOp.html);
the pipeline uses it semi-interchangeably with the [`magrittr` pipe
operator](https://magrittr.tidyverse.org/), `%>%`, due to some of the
early code being written before `|>` was implemented in R version 4.1.

The multiline comment above the target is a convention used in this
pipeline for readability. If the code is opened in Rstudio, the name of
each target will show up in the quick navigation bar at the bottom of
the code window.

`{targets}` also includes powerful features for *branching*, i.e., using
one target definition multiple times for different data. These are too
complex to fully describe here, but anyone interested in understanding
the pipeline is referred to the [`{targets}` user
manual](https://books.ropensci.org/targets). In particular, the
`tar_map()` function, which typically encompasses several target
definitions, is used for for [static
branching](https://books.ropensci.org/targets/static.html), and the
`pattern` and `iteration` arguments to `tar_target()` (and similar
functions) are used for [dynamic
branching](https://books.ropensci.org/targets/dynamic.html)

## Execution

Execution of the pipeline on a single Linux computer using the included
[conda](https://conda.io) environment file is probably feasible.
Adaption to a different SLURM cluster would take additional effort, and
a different HPC or cloud system (SGE, AWS, etc.) even more. No tests
have been attempted on OSX or Windows.

In principle the steps to execute would be:

1)  Unzip the archive (if you are reading this, you have probably
    already accomplished this step.)

``` sh
tar -xzf GSSP_pipeline.tar.gz
cd GSSP
```

2)  Install all prerequisites listed in `conda/conda.yaml`, either using
    [conda](https://conda.io) or by other means. If using conda,
    activate the conda environment.

``` sh
conda env create -f GSSP/conda/GSSP.yaml
conda activate GSSP
```

3)  Download ProtaxFungi and unzip it in the main GSSP directory (or a
    sister directory).

``` sh
tar -xzf protaxFungi.tar.gz
```

4)  Download the raw sequence files from ENA project
    [PRJEB65748](https://www.ebi.ac.uk/ena/browser/view/PRJEB65748) into
    `sequences/01_raw/`.

5)  Download [USEARCH](https://www.drive5.com/usearch/) version 11.0.667
    and unzip in `bin/`.

``` sh
wget -nd -P bin/ https://www.drive5.com/downloads/usearch11.0.667_i86linux32.gz
gunzip bin/usearch11.0.667_i86linux32.gz
chmod +x bin/usearch11.0.667_i86linux32
```

6)  download [reference data for the Unite sh_matching
    pipeline](https://files.plutof.ut.ee/public/orig/9C/FD/9CFD7C58956E5331F1497853359E874DEB639B17B04DB264C8828D04FA964A8F.zip)
    and unzip in `data/sh_matching_data`.

``` sh
wget https://files.plutof.ut.ee/public/orig/9C/FD/9CFD7C58956E5331F1497853359E874DEB639B17B04DB264C8828D04FA964A8F.zip
unzip -j 9CFD7C58956E5331F1497853359E874DEB639B17B04DB264C8828D04FA964A8F.zip data/shs_out.txt data/sanger_refs_sh.fasta -d data/sh_matching_data
rm 9CFD7C58956E5331F1497853359E874DEB639B17B04DB264C8828D04FA964A8F.zip
```

6)  In R (for execution on a single computer):

``` r
library(targets)
tar_make()
```

Alternatively, `sbatch run_node.sh` submits the pipeline as a single job
using 8 cores to a SLURM cluster. `sbatch run_clustermq.sh` submits a
“master” job using a single core, which will itself submit an array job
requesting one 8-core worker per sequencing run (i.e., 33 workers).
Adapting these to run on a different cluster would involve modifying
`run_node.sh` and/or `run_clustermq.sh`, as well as potentially loading
environment modules, a conda environment, etc. For ClusterMQ execution,
a new template like `slurm/puhti_clustermq.tmpl` would also be required,
as well as modifications to `run_clustermq.R`.

Bioinformatics pipeline for "Wind is a primary driver of fungal dispersal across a mainland-island system"
================

This pipeline is a pre-release version of the OptimOTU pipeline, which is
implemented in R using the [{targets}](https://books.ropensci.org/targets/)
pipeline tool.
It is probably difficult to run it outside of the [CSC](https://www.csc.fi)
Puhti computing cluster where it was developed and executed for this project,
and it is included primarily for reference purposes to accompany the paper.

In principle the steps to execute would be:

1)  Unzip the archive (if you are reading this, you have probably
    already accomplished this step.)

```sh
unzip konnevesi_pipeline.zip
cd konnevesi
```

2)  Install all prerequisites listed in `conda/GSSP.yaml`, either using
    [conda](https://conda.io) or by other means. If using conda,
    activate the conda environment.

```sh
conda env create -f GSSP/conda/GSSP.yaml
conda activate GSSP
```

3)  Download ProtaxFungi and unzip it in the main konnevesi directory (or a
    sister directory).

```sh
tar -xzf protaxFungi.tar.gz
```

4)  Download the raw sequence files from ENA project
    [PRJEB76596](https://www.ebi.ac.uk/ena/browser/view/PRJEB76596) into
    `sequences/01_raw/KN_22`.

5)  Download [USEARCH](https://www.drive5.com/usearch/) version 11.0.667
    and unzip in `bin/`.

```sh
wget -nd -P bin/ https://www.drive5.com/downloads/usearch11.0.667_i86linux32.gz
gunzip bin/usearch11.0.667_i86linux32.gz
chmod +x bin/usearch11.0.667_i86linux32
```

6)  download [reference data for the Unite sh_matching
    pipeline](https://files.plutof.ut.ee/public/orig/9C/FD/9CFD7C58956E5331F1497853359E874DEB639B17B04DB264C8828D04FA964A8F.zip)
    and unzip in `data/sh_matching_data`.

```sh
wget https://files.plutof.ut.ee/public/orig/9C/FD/9CFD7C58956E5331F1497853359E874DEB639B17B04DB264C8828D04FA964A8F.zip
unzip -j 9CFD7C58956E5331F1497853359E874DEB639B17B04DB264C8828D04FA964A8F.zip data/shs_out.txt data/sanger_refs_sh.fasta -d data/sh_matching_data
rm 9CFD7C58956E5331F1497853359E874DEB639B17B04DB264C8828D04FA964A8F.zip
```

6)  In R (for execution on a single computer):

```r
library(targets)
tar_make()
```

Alternatively, `sbatch run_node.sh` submits the pipeline as a single job
using 8 cores to a SLURM cluster. `sbatch run_clustermq.sh` and
`sbatch run_crew.sh` submit a “master” job using a single core, which
will itself submit additional "worker" jobs using 8 cores each. The
clustermq backend submits a fixed number of workers as an array job,
while the crew backend dynamically submits workers depending on the
pipeline's needs.  Adapting these to run on a different cluster would
involve modifying `run_node.sh`, `run_clustermq.sh`, and/or `run_crew.sh`,
as well as potentially loading environment modules, a conda environment,
etc. For ClusterMQ or crew execution, a new template like
`slurm/puhti_clustermq.tmpl` and `slurm/puhti_crew.tmpl` would also be
required, as well as modifications to `run_clustermq.R` or `run_crew.R`.


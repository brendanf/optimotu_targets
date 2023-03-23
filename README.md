ITS2 taxonomy-first metabarcoding pipeline
================

<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->
<!-- badges: end -->

### Installation

- [ ] clone the repository from github;
  `git clone git@github.com:brendanf/sonja_spruce_logs.git`  
  This will by default create a directory called `sonja_spruce_logs` for
  your project, but you can put it in a different directory using
  `git clone git@github.com:brendanf/sonja_spruce_logs.git {name_of_directory}`

- [ ] (recommended) start a new branch; if you will do several projects,
  it is recommended to put each branch in a separate directory using
  `git worktree`. This allows them to share one copy of the `.git`
  directory, as well as `protaxFungi`.

  ``` sh
  cd sonja_spruce_logs #or another directory name if you chose one
  git worktree add ../{name_of_project}
  ```

  **Or**, if you only plan on doing one project, you can just create the
  branch in the main directory:

  ``` sh
  cd sonja_spruce_logs #or another directory name if you chose one
  git checkout -b {name_of_project}
  ```

- [ ] download [protaxFungi](https://github.com/psomervuo/protaxfungi)
  and unzip into a sister directory (i.e., `../protaxFungi`). **This
  currently does not work because the github version of protaxfungi is
  not the lastest version.**

  ``` sh
  cd ..
  wget https://github.com/psomervuo/protaxfungi/raw/master/protaxfungi.tgz
  tar -xf protaxfungi.tgz
  ```

  To get the current version, on puhti (Finnish computing cluster)
  *only*:

  ``` sh
  cd ..
  cp -r /scratch/project_2003156/protaxFungi .
  ```

- [ ] download [usearch](https://drive5.com/usearch/) into `bin/`, and
  rename or link the executable to be called just `bin/usearch`. Make
  sure it is executable.

  ``` sh
  cd {name_of_project}/bin
  wget https://www.drive5.com/downloads/usearch11.0.667_i86linux32.gz
  gunzip usearch11.0.667_i86linux32.gz
  ln -s usearch11.0.667_i86linux32 usearch
  chmod +x usearch
  cd ..
  ```

- [ ] download [reference data for the Unite sh_matching
  pipeline](https://files.plutof.ut.ee/public/orig/9C/FD/9CFD7C58956E5331F1497853359E874DEB639B17B04DB264C8828D04FA964A8F.zip)
  and unzip in `data/sh_matching_data`. (Currently only
  `sanger_refs_sh.fasta` and `shs_out.txt` are used, so you can delete
  everything else.)

  ``` sh
  wget https://files.plutof.ut.ee/public/orig/9C/FD/9CFD7C58956E5331F1497853359E874DEB639B17B04DB264C8828D04FA964A8F.zip
  unzip -j 9CFD7C58956E5331F1497853359E874DEB639B17B04DB264C8828D04FA964A8F.zip data/shs_out.txt data/sanger_refs_sh.fasta -d data/sh_matching_data
  rm 9CFD7C58956E5331F1497853359E874DEB639B17B04DB264C8828D04FA964A8F.zip
  ```

- [ ] download or link your own *demultiplexed, paired-end* Illumina
  reads into `sequences/01_raw`. Each sequencing run should be in its
  own directory.

  ``` sh
  mkdir -p sequences/01_raw
  cd sequences/01_raw
  # do your download here
  cd ../..
  ```

  The actual download might be done with `cp` or `ln -s` if it is
  already somewhere else on your system; `wget` for a general HTTP or
  FTP address; `bs download project --name` from Illumina BaseSpace;
  `rclone copy` from AWS or another object storage platform
  (e.g. Allas).

- [ ] (optional) if there are particular species that you are interested
  in, write their full species identifiers as used in ProtaxFungi
  (`Genus_species_mycobank#`, e.g. `Rhodofomes_roseus_127496`), one per
  line, in a file `data/target_taxa.txt`.

- [ ] (optional) if you have additional reference sequences for species
  of interest (which may or may not also be in `target_taxa.txt`), then
  you can create two files: `data/culture_refs.fasta` with the actual
  sequences, and `data/culture_sequences.xlsx` with medadata. In
  particular, `data/culture_sequences.clsx` must include, in the first
  worksheet, a column “Culture ID” giving the name of each reference
  sequence exactly as it appears in `data/culture_refs.fasta`, and a
  column
  “Protax_synonym`giving the full, comma-separated taxonomy for each reference in Protax format, e.g.`Fungi,Basidiomycota,Agaricomycetes,Polyporales,Fomitopsidaceae,Antrodia_17083,Antrodia_piceata_813073`.  Adding reference sequences which are not annotated fully to species (e.g.,`Fungi,Basidiomycota,Agaricomycetes,Polyporales,Fomitopsidaceae\`)
  *likely* works but has not been tested.

- [ ] Install R dependencies. This can be done in several ways:  
  (using renv)

  ``` sh
  R
  ```

  then

  ``` r
  #if renv is not installed
  install.packages("renv")
  renv::init()
  # this will have failed to install github and bioconductor packages
  renv::install("brendanf/optimotu", "bioconductor::dada2")
  renv::snapshot()
  ```

  (using conda)

  ``` sh
  conda env create -f conda/its2_taxonomy_first.yaml
  conda activate its2_taxonomy_first
  ```

  (using tykky on puhti)

      module load tykky
      cd /projapps/{your_project}
      conda-containerize --prefix /projappl/{your_csc_project}/{your_project} conda/its2_taxonomy_first.yaml
      export PATH="/projappl/{your_csc_project}/{your_project}:$PATH"

  (using existing tykky container on puhti)

      export PATH="/projappl/project_2003156/its2_taxonomy_first/bin:$PATH"

### Execution

#### Local exection (command prompt)

Run the full pipeline:

``` sh
R -e 'targets::tar_make(callr_function=NULL)'
```

Test that samples are correctly detected:

``` sh
Rscript _targets.R
```

Check what targets need to be run:

``` sh
R -e 'targets::tar_outdated(callr_function=NULL)'
```

Run to generate a specific target:

``` sh
R -e 'targets::tar_make({name_of_target}, callr_function=NULL)'
```

Check if a specific target (or its dependencies) need to be run:

``` sh
R -e 'targets::tar_outdated({name_of_target}, callr_function=NULL)'
```

#### Local execution (from R)

Run the full pipeline:

``` r
targets::tar_make()
```

Test that samples are correctly detected:

``` r
source("_targets.R")
# you can also interactively investigate the sample_table object
sample_table
sample_table$sample
sample_table$seq_run
# etc.
```

Check what targets need to be run:

``` r
targets::tar_outdated()
```

Run to generate a specific target:

``` r
targets::tar_make({name_of_target})
```

Check if a specific target (or its dependencies) need to be run:

``` r
targets::tar_outdated({name_of_target})
```

#### Cluster execution (Puhti, using a single node)

You may need to modify `run_node.sh`, for instance to change the project
or tykky container.

Run the full pipeline:

``` sh
sbatch run_node.sh
```

Test that samples are correctly detected (on login node):

``` sh
# first line only needed once per session
export PATH="/projappl/{your_csc_project}/{your_project}:$PATH"

Rscript _targets.R
```

Check what targets need to be run (on login node):

``` sh
bash run_node.sh test
```

Run to generate a specific target:

``` sh
sbatch run_node.sh {name_of_target}
```

Check if a specific target (or its dependencies) need to be run:

``` sh
bash run_node.sh test {name_of_target}
```

#### Cluster execution (parallel on multiple nodes)

You may need to modify `run_clustermq.sh`, `run_clustermq.R` and
`slurm/puhti_clustermq.tmpl`, for instance to change the project or
tykky container (`.sh` and `.tmpl`) or the number of workers (`.R`).

Run the full pipeline:

``` sh
sbatch run_clustermq.sh
```

Test that samples are correctly detected (on login node):

``` sh
# first line only needed once per session
export PATH="/projappl/{your_csc_project}/{your_project}:$PATH"

Rscript _targets.R
```

Check what targets need to be run (on login node):

``` sh
bash run_node.sh test
```

Run to generate a specific target:

``` sh
sbatch run_clustermq.sh {name_of_target}
```

Check if a specific target (or its dependencies) need to be run (on
login node):

``` sh
bash run_node.sh test {name_of_target}
```

#### Cluster execution (other)

It should be possible to run on other Slurm-based HPC systems by
additional modification of `run_node.sh`, `run_clustermq.sh`,
`run_clustermq.R`, and `slurm/puhti_clustermq.tmpl`. The least portable
element is the tykky containerization. In the future a Singularity
container will be provided to make installation on other systems easier.

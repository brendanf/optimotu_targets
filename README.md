OptimOTU pipeline
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- badges: start -->

<!-- badges: end -->

### Installation

- [ ] If you do not already have them, you need to [install
  ‘git’](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git),
  create a [GitHub account](https://github.com), and set up [SSH
  keys](https://docs.github.com/en/authentication/connecting-to-github-with-ssh)
  to allow you to clone repositories using SSH.

- [ ] clone the repository from github;
  `git clone git@github.com:brendanf/optimotu_targets.git`\
  This will by default create a directory called `optimotu_targets` for
  your project, but you can put it in a different directory using
  `git clone git@github.com:brendanf/optimotu_targets.git {name_of_directory}`

- [ ] (recommended) start a new branch; if you will do several projects,
  it is recommended to put each branch in a separate directory using
  `git worktree`. This allows them to share one copy of the `.git`
  directory, as well as `protaxFungi`.

  ``` sh
  cd optimotu_targets #or another directory name if you chose one
  git worktree add ../{name_of_project}
  ```

  **Or**, if you only plan on doing one project, you can just create the
  branch in the main directory:

  ``` sh
  cd optimotu_targets #or another directory name if you chose one
  git checkout -b {name_of_project}
  ```

- [ ] (Optional: for ITS fungal metabarcoding) download
  [protaxFungi](https://github.com/psomervuo/protaxfungi) and unzip into
  a sister directory (i.e., `../protaxFungi`). **This currently does not
  work because the github version of protaxfungi is not the lastest
  version.**

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
  sure it is executable. Note that Usearch 12 is not supported.

  ``` sh
  cd {name_of_project}/bin
  wget https://www.drive5.com/downloads/usearch11.0.667_i86linux32.gz
  gunzip usearch11.0.667_i86linux32.gz
  ln -s usearch11.0.667_i86linux32 usearch
  chmod +x usearch
  cd ..
  ```

- [ ] download [reference data for the Unite sh_matching
  pipeline](https://s3.hpc.ut.ee/plutof-public/original/5dcbe93a-b50c-4ff5-8898-528160c4e593.zip)
  and unzip in `data/sh_matching_data`. (Currently only
  `sanger_refs_sh.fasta` and `shs_out.txt` are used, so you can delete
  everything else.)

  ``` sh
  wget https://s3.hpc.ut.ee/plutof-public/original/5dcbe93a-b50c-4ff5-8898-528160c4e593.zip
  unzip -j 5dcbe93a-b50c-4ff5-8898-528160c4e593.zip data/shs_out.txt data/sanger_refs_sh.fasta -d data/sh_matching_data
  rm 5dcbe93a-b50c-4ff5-8898-528160c4e593.zip
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
  column “Protax_synonym” giving the full, comma-separated taxonomy for
  each reference in Protax format,
  e.g. `Fungi,Basidiomycota,Agaricomycetes,Polyporales,Fomitopsidaceae,Antrodia_17083,Antrodia_piceata_813073`.
  Adding reference sequences which are not annotated fully to species
  (e.g.,
  `Fungi,Basidiomycota,Agaricomycetes,Polyporales,Fomitopsidaceae`)
  *likely* works but has not been tested.

- [ ] Install R dependencies. This can be done in several ways:\
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
  renv::install("brendanf/optimotu", "brendanf/optimotu.pipeline", "bioconductor::dada2")
  renv::snapshot()
  ```

  (using conda)

  ``` sh
  conda env create -f conda/OptimOTU_v6.yaml
  conda activate OptimOTU_v6
  ```

  The conda environment file listed above includes only direct
  dependencies of OptimOTU. Occasionally a new version of an indirect
  dependency may cause something to break. In that case, you can instead
  use `conda env create -f conda/OptimOTU_v6_full.yaml`.

  (using tykky on puhti)

      module load tykky
      mkdir /projappl/{your_csc_project}/OptimOTU_v6
      conda-containerize new --prefix /projappl/{your_csc_project}/OptimOTU_v6 --post-install conda/OptimOTU_v6_postinstall.sh conda/OptimOTU_v6.yaml
      export PATH="/projappl/{your_csc_project}/{your_project}:$PATH"

  (using existing tykky container on puhti, if you have access to it)

      export PATH="/projappl/project_2005718/OptimOTU_v6/bin:$PATH"

### Configuration

Configuration of the pipeline is done through a (more or less)
human-readable file called `pipeline_options.yaml`. Comments in the file
(Beginning with `#`) explain the function of most of the options. Listed
here are a few types of options you may want/need to change depending on
different characteristics of your data. Also see the example files
`pipeline_options_example_fungi.yaml` and
`pipeline_options_example_metazoa.yaml`.

#### Basic option: project name

The first option is `project_name`. The default option is
“metabarcoding_project”, but the pipeline will issue a warning if you
don’t change it. The only thing this option actually controls is the
name of the zip file containing all outputs for a run. As such it should
be a good file name (ideally underscores instead of spaces, not too
long). It could be (but does not have to be) the same name as your
project directory/git branch.

#### File names and orientation.

By default, the OptimOTU pipeline auto-detects sample and sequencing-run
names based on file names. It strips some common Illumina elements from
the file names and assumes that what is left is the sample name; and it
assumes the samples are saved in directories which correspond to
sampling runs.

If your situation is more complex than this, you may need to use a
*custom sample table*. This is a TSV file which includes, at a minimum,
columns named “`sample`”, “`seqrun`”, “`fastq_R1`”, and “`fastq_R2`”.
These give, for each sample, the sequencing run it came from, and the
file paths for the R1 and R2 fastq files relative to the directory
`sequences/01_raw`.

One possible reason to use a custom sample table is to define the
orientation of reads. Some lab workflows result in read pairs which all
have the same orientation, while others result in read pairs which are
in mixed orientation; finally others (in order to increase sequence
diversity) are all in the same orientation within each sample, but this
varies between samples. All of these cases are supported by the `orient`
option; the case where orientations vary between samples, however,
requires you to designate the orientation of each sample, which can be
done using a custom sample table with an “`orient`” column.

#### Different primer pairs

The default `pipeline_options.yaml` file (and also
`pipeline_options_example_fungi.yaml`)is set up for ITS2 amplicons
generated using the primer pair ITS3-ITS4.
`pipeline_options_example_metazoa.yaml` is set up for COI amplicons
using the primer pair BF3-BR2. The following options need to be changed
when switching to a different primer pair:

- `forward_primer` and `reverse_primer` should have the sequences of the
  primers. Ambiguous bases (N, C, Y, M, R, …) are OK.

- The settings under `amplicon_model` will need to be changed. Full
  instructions for generating a model (HMM or CM) for a new amplicon is
  beyond the scope of this README, although the process will be
  automated in a future version of OptimOTU. The easiest solution is to
  disable model alignment entirely, by setting `model_type` to “`none`”.

  - `model_type`: “`CM`” is specifically for non-coding RNA genes (like
    ITS, 16S, LSU, …). They can be created using the program `cmbuild`
    from [Infernal](http://eddylab.org/infernal/). “`HMM`” is a more
    general type of model which can be used for non-coding RNA, but also
    for protein-coding genes or sequences in general. They can be
    created using the program `hmmbuild` from
    [Hmmer](http://hmmer.org/).
  - `model_file`: location of the HMM or CM file, relative to the root
    of the project
  - `model_filter`: these parameters should be adjusted when using a new
    CM or HMM. Especially the `min_model_end` (i.e., how far in the
    model the amplicon must go) and the `min_model_score` (i.e., how
    well must the amplicon match the model) are typically higher for
    longer models.

#### Different versions of Protax

OptimOTU currently uses two different versions of Protax for taxonomic
identification: ProtaxFungi and ProtaxAnimal. Although these are named
for the taxonomic group they were designed for, they are also completely
different programs which use a different file structure for their
trained models. It is (in theory) possible to train a model for the
ProtaxFungi software for any amplicon and any group of organisms, given
a defined taxonomy and a set of labeled reference sequences.
ProtaxAnimal is more restrictive but also much faster. It requires that
all input sequences, both references and queries, are globally aligned.
This alignment can be performed in the OptimOTU pipeline by setting
`model_align` (under `amplicon_model`) and `aligned` (under `protax`)
both to “`yes`”. The location of the Protax installation directory
should also be given as `location` under `protax`.

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

Before running, modify `run_node.sh` to change the lines:

``` sh
#SBATCH --account project_2003104
```

and

``` sh
export PATH="/projappl/project_2005718/OptimOTU_v6/bin:$PATH"
```

to use your own CSC project name.

Run the full pipeline:

``` sh
sbatch run_node.sh
```

Test that samples are correctly detected (on login node):

``` sh
# first line only needed once per session
export PATH="/projappl/{your_csc_project}/OptimOTU_v6/bin:$PATH"

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

Before the first time you run, modify `run_crew.sh` and
`slurm/puhti_crew.tmpl` to change the lines

``` sh
#SBATCH --account project_2003104
```

and

``` sh
export PATH="/projappl/project_2005718/OptimOTU_v6/bin:$PATH"
```

to refer to your own CSC project.

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
additional modification of `run_node.sh`, `run_crew.sh`, `run_crew.R`,
and `slurm/puhti_crew.tmpl`. The least portable element is the tykky
containerization. In the future a Singularity container will be provided
to make installation on other systems easier.

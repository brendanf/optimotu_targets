# optimotu_targets 7.0

## Major User-facing Changes

- Add optional UNOISE (vsearch) denoising as an alternative to DADA2,
  configured via `denoising:` and merged-read `filtering:` keys in
  `pipeline_options.yaml`.
- Promote pair-merging options to top-level `merging:` in
  `pipeline_options.yaml` (shared by DADA2 and UNOISE; DADA2 merge parameters
  are configurable for the first time). Add top-level `dist_config:` and
  `executables:`; move `added_reference` under `taxonomy.protax`.
- Add option `supplemental_asv` to merge external final-ASV sets (with optional
  abundances and taxonomy) into clustering and downstream outputs. Supplemental
  ASVs are never rarefied.
- Clustering thresholds can now be optimized from `pipeline_options.yaml`
  (reference data, this dataset, or a custom FASTA), not only loaded from a
  pre-computed file.
- `guilds` in `pipeline_options.yaml` can list annotation databases
  (`funguild`, `carlos`, and/or `{name:, file:}`); `guilds: yes` keeps the
  FUNGuild + Carlos defaults.
- Selections in `pipeline_options.yaml` `output.formats` are now documented
  and implemented.
- Breaking: rename per-read fate target `dada_map` → `read_map`, and phase-1
  plan object `dada_plan` → `phase1_plan`. A shared `denoise_map` now matches
  sequences once per chunk for both community tables and read maps.
- The `asv_tax_prob` table is now in long format (`seq_id`, `rank`, `taxon`,
  `prob`) rather than one column per rank; update any code that reads this
  output.

## Improvements

- Krona plots and guild assignment work with arbitrary taxonomic rank sets.
- LULU post-clustering curation uses a low-memory implementation that scales
  better on large sequencing runs.
- Per-read fate maps (`read_map`) are remapped through LULU before UNCROSS so
  sequence indices match post-LULU / post-UNCROSS abundance tables.
- Hamming clustering defaults to SLINK with merge-based parallelization;
  other distance methods keep tree clustering with concurrent parallelization.
- Document `ksw2` as a `dist_config` distance method in YAML comments.
- OTU reference FASTA outputs now include aligned sequences; ASV FASTAs are
  indexed for faster sequence extraction.
- Output directory and zipped archive now include reproducibility metadata: a
  copy of `pipeline_options.yaml`, git commit hash, uncommitted changes diff,
  `sessionInfo.txt`, and the custom sample table when one is used.
- Route `crew` Slurm worker logs to per-job directories with explicit stdout and
  stderr files to simplify troubleshooting.
- Slurm runs use a second “thin” worker pool for lighter parallel tasks.
- Add versioned Apptainer definition files and a helper build script with
  optional host renv-cache reuse during image builds.

## Bug Fixes

- Fix Infernal covariance-model alignment handling that could fail on current
  Infernal object shapes.
- Fix de novo singleton handling so an ASV that is the only remaining unknown
  in a parent taxon after closed-reference clustering still receives a
  pseudotaxon (and therefore an OTU) instead of being dropped.
- Fix model-based ASV filtering so sequences retained after filtering are
  tracked consistently through clustering, abundance tables, and outputs.
- Fix Protax and BayesANT taxonomy assignment (including indexed Protax input
  and BayesANT threshold optimization).
- Fix missing singleton de novo clusters when using `force_denovo`.
- Fix mapping from candidate to final ASVs, which could mis-assign abundances
  and taxonomy in some cases.
- Forced de novo pseudotaxa that combine known ASVs from different taxa are
  now classified as uncertain rather than known.
- Skip reference-sequence model generation when no amplicon model is
  configured, avoiding errors in that case.
- Fix startup check that compared the `optimotu` version against the wrong
  minimum version.

Requires `optimotu` 0.10.0+ and `optimotu.pipeline` 0.7.0+.

# optimotu_targets 6.0.1

- Add option `force_denovo` to the `clustering` section, to force de-novo
  clustering for certain taxonomic ranks.
- Fix LULU for model-aligned amplicons with Hamming distance.
- Use `optimotu.pipeline` version 0.6.2, which has an important fix for LULU
  implementation.

# optimotu_targets 6.0.0

- Use `FUNGuildR` v0.3.0, which has an important bug fix.
- Parsing of `pipeline_options.yaml` and storage of most global options have
  been ported to the `optimotu.pipeline` package.
- Add new taxonomic classifier options BayesANT, SINTAX, and EPA-ng.
- All clustering is now performed through the `optimotu` package.
- Unconventional situations where some sequencing runs have all reads oriented
  the same, while other sequencing runs are mixed, now work.
- Read-based rarefaction is now supported via the `rarefy` option in
  `pipeline_options.yaml`
- The LULU algorithm for post-clustering curation is now implemented.


# optimotu_targets 5.1.0

- Quality filtering parameters `maxEE_R1` and `maxEE_R2` can now be given
  sample-specific values in a custom sample table using columns with those
  names.
- Fix several errors occurring with empty or almost empty samples/sequencing
  runs.
- Fix an error when using `dense_table: yes` in combination with
  `orient: mixed`.

# optimotu_targets 5.0.0

- Update conda environment to v5.
- Improved handling of very large outgroup reference files.
- Add option `local_threads` to specify the maximum number of threads to use
  in local execution.

# optimotu_targets 4.1.1

- Fixed a bug which caused an error when using positive control sequences.

# optimotu_targets 4.1.0

- Detection of spike sequences is now optional, and an option is also included
  for positive control sequences.  Both of these options are under `controls:`
  in the configuration file.  If either type of sequences should be detected,
  the file containing sequences to detect should be given as `controls: spikes:`
  or `controls: positive:` respectively. Note this is a breaking change for
  projects that did in fact use the default synmock spikes. These are included
  in protaxFungi as `protaxFungi/amptk_synmock.fasta`, but must now be
  explicitly specified.

# optimotu_targets 4.0.1

- Improved parsing of BOLD release datasets when used as outgroup references.
  It is now possible to directly download the `BOLD_Public.DD-Mmm-YYYY.fasta.gz`
  to supply as the file for `outgroup_reference:sequences:`. No separate
  taxonomy file is required. Note that execution will be a bit faster if the
  FASTA file is pre-filtered to include only sequences with "COI-5P" in the
  header, but this is not crucial, since this is more than 80% of the
  sequences in recent BOLD snapshots.

# optimotu_targets 4.0.0

Initial public release, as described in arXiv preprint.

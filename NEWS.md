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

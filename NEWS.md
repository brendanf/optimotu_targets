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

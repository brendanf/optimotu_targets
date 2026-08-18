

#### samplewise_plan ####
# The sample-wise plan consists of targets where individual samples are
# processed separately, but there may be interactions between different reads
# in the same sample. Thus these targets are not independent of rarefaction.
samplewise_plan <- c(
  list(
    ##### samplewise_meta_{.orient?}_{.seqrun}_{.rarefaction?}_{.replicate?} #####
    # grouped tibble:
    #  `seqrun` character; name of sequencing run (directory in sequences/01_raw)
    #  `sample` character; name of sample, based on parsing file name
    #  `readwise_key`character; common prefix of trim_R1, trim_R2, filt_R1 and
    #      filt_R2
    #  `sample_key` character; as `sample_key` but also with rarefy_text
    #  `fastq_R1` character; file name with path for raw R1 file
    #  `trim_R1` character; file name with path for trimmed R1 file
    #  `filt_R1` character; file name with path for filtered R1 file
    #  `filt_R2` character; file name with path for filtered R2 file
    #  `to_denoise_R1` character; file name with path for R1 to denoise
    #  `to_denoise_R2` character; file name with path for R2 to denoise
    #  `rarefy_text` (optional) character; specification of the rarefaction
    #  `numerator` (optional) integer; numerator for fractional rarefacation
    #  `denominator` (optional) integer; denominator for fractional rarefacation
    #  `number` (optional) integer; count for count.based rarefaction
    #  `tar_seed` integer; random seed to use for dereplication
    samplewise_meta = tar_target(
      samplewise_meta,
      dplyr::select(readwise_meta, orient, readwise_key) |>
        dplyr::left_join(sample_table, by = c("orient", "readwise_key")) |>
        dplyr::filter(
          !!(if (optimotu.pipeline::do_rarefy()) {
            quote(rarefy_text == .rarefy_text)
          } else {
            TRUE
          })
        ) |>
        dplyr::semi_join(filt_read_counts, by = "filt_R1") |>
        dplyr::select(
          seqrun,
          sample,
          readwise_key,
          sample_key,
          fastq_R1,
          trim_R1,
          filt_R1,
          filt_R2,
          to_denoise_R1,
          to_denoise_R2,
          any_of(c(
            "rarefy_text",
            "numerator",
            "denominator",
            "number",
            "tar_seed"
          ))
        ),
      pattern = map(readwise_meta)
    )
  ),

  ##### map over R1 and R2 #####
  # inside the tar_map, every occurrence of `read` is replaced by "R1" or "R2"
  # the read name is also appended to all target names
  # so all of this gets done separately for forward and reverse reads.
  # pattern=map() means we are also keeping the different sequencing runs
  # separate
  tar_map(
    values = list(read = c("R1", "R2")),

    ###### predenoise_{read}_{.orient?}_{.seqrun}_{.rarefaction?}_{.replicate?} ######
    # character: path and file name of filtered reads; fastq.gz
    #
    # select only the files corresponding to the read we are working on
    predenoise = if (optimotu.pipeline::do_rarefy()) {
      tar_file(
        predenoise,
        mapply(
          optimotu.pipeline::fastq_sample,
          infile = samplewise_meta[[paste0("filt_", read)]],
          outfile = samplewise_meta[[paste0("to_denoise_", read)]],
          n = !!(if (is.null(optimotu.pipeline::rarefy_number())) {
            quote(round(
              .numerator * filt_read_counts$filt_nread / .denominator
            ))
          } else {
            quote(.number)
          }),
          sample = mapply(
            \(n, seed) {
              tar_seed_set(seed)
              sample(n)
            },
            filt_read_counts$filt_nread,
            samplewise_meta$tar_seed,
            SIMPLIFY = FALSE
          ),
          SIMPLIFY = TRUE
        ),
        pattern = map(samplewise_meta, filter_pairs, filt_read_counts),
        resources = tar_resources(
          crew = tar_resources_crew(controller = "thin")
        )
      )
    } else {
      tar_file(
        predenoise,
        purrr::keep(filter_pairs, endsWith, paste0(read, "_filt.fastq.gz")),
        pattern = map(filter_pairs),
        deployment = "main"
      )
    },

    ###### derep_{read}_{.orient?}_{.seqrun}_{.rarefaction?}_{.replicate?} ######
    # list of dada2 `derep` objects
    #
    # dereplicate
    derep = tar_target(
      derep,
      optimotu.pipeline::derepFastq(
        predenoise,
        verbose = TRUE,
        names = optimotu.pipeline::file_to_sample_key(predenoise)
      ),
      pattern = map(predenoise),
      resources = tar_resources(crew = tar_resources_crew(controller = "wide")) # memory
    ),

    ###### err_{read}_{.orient?}_{.seqrun}_{.rarefaction?}_{.replicate?} ######
    # list: see dada2::LearnErrors
    #
    # fit error profile
    err = tar_target(
      err,
      optimotu.pipeline::learnErrors(
        # TODO: which samples to train the model on?
        purrr::discard(predenoise, grepl, pattern = "BLANK|NEG"),
        errorEstimationFunction = errfun,
        multithread = optimotu.pipeline::local_cpus(),
        verbose = TRUE
      ),
      resources = tar_resources(crew = tar_resources_crew(controller = "wide"))
    ),

    ###### denoise_{read}_{.orient?}_{.seqrun}_{.rarefaction?}_{.replicate?} ######
    # list of dada2 `dada` objects
    denoise = tar_target(
      denoise,
      optimotu.pipeline::dada(
        derep,
        err = err,
        errorEstimationFunction = errfun,
        multithread = optimotu.pipeline::local_cpus(),
        verbose = TRUE
      ),
      pattern = map(derep),
      resources = tar_resources(crew = tar_resources_crew(controller = "wide"))
    )
  ),

  ##### merged_{.orient?}_{.seqrun}_{.rarefaction?}_{.replicate?} #####
  # list of data.frame; see dada2::mergePairs
  #
  # Merge paired reads and make a sequence table for each sequencing run
  merged = tar_target(
    merged,
    optimotu.pipeline::mergePairs(
      denoise_R1,
      derep_R1,
      denoise_R2,
      derep_R2,
      minOverlap = 10,
      maxMismatch = 1,
      verbose = TRUE
    ),
    pattern = map(denoise_R1, derep_R1, denoise_R2, derep_R2),
    resources = tar_resources(crew = tar_resources_crew(controller = "wide"))
  ),

  ##### seqtable_raw_{.orient?}_{.seqrun}_{.rarefaction?}_{.replicate?} #####
  # `tibble` with columns:
  #   `sample` (character) sample name as given in sample_table$sample_key
  #   `seq_idx` (integer) index of a sequence in seq_all
  #   `nread` (integer) number of reads
  seqtable_raw = tar_fst_tbl(
    seqtable_raw,
    optimotu.pipeline::make_mapped_sequence_table(
      merged,
      seq_all,
      rc = .orient == "rev"
    ),
    pattern = map(merged),
    resources = tar_resources(crew = tar_resources_crew(controller = "thin"))
  ),

  ##### dada_map_{.orient?}_{.seqrun}_{.rarefaction?}_{.replicate?} #####
  # map the raw reads to the nochim ASVs
  # indexes are per-sample
  #
  # a tibble:
  #  sample: character identifying the sample, as in sample_table
  #  raw_idx: integer index of read in the un-rarified fastq file
  #  seq_idx: integer index of ASV in seq_all
  #  flags: raw, bits give presence/absence of the read after different stages:
  #    0x01: trim
  #    0x02: filter
  #    0x04: denoise & merge
  #    0x08: tag-jump removal (if performed)
  dada_map = tar_target(
    dada_map,
    mapply(
      FUN = optimotu.pipeline::seq_map,
      sample = samplewise_meta$sample_key,
      fq_raw = samplewise_meta$fastq_R1,
      fq_trim = samplewise_meta$trim_R1,
      fq_filt = samplewise_meta$filt_R1,
      dadaF = denoise_R1,
      derepF = derep_R1,
      dadaR = denoise_R2,
      derepR = derep_R2,
      merged = merged,
      MoreArgs = list(
        seq_all = seq_all,
        rc = .orient == "rev"
      ),
      SIMPLIFY = FALSE
    ) |>
      purrr::list_rbind(
        ptype = tibble::tibble(
          sample = character(),
          raw_idx = integer(),
          seq_idx = integer(),
          flags = raw()
        )
      ),
    pattern = map(
      samplewise_meta,
      denoise_R1,
      derep_R1,
      denoise_R2,
      derep_R2,
      merged
    ),
    resources = tar_resources(crew = tar_resources_crew(controller = "wide")) # for memory
  )
)

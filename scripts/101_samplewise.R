#### samplewise_plan ####
# The sample-wise plan consists of targets where individual samples are
# processed separately, but there may be interactions between different reads
# in the same sample. Thus these targets are not independent of rarefaction.
samplewise_plan <- c(
  if (optimotu.pipeline::do_dada2()) {
    c(
      list(
        ##### samplewise_meta_{.orient?}_{.seqrun}_{.rarefaction?}_{.replicate?} #####
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
      tar_map(
        values = list(read = c("R1", "R2")),

        ###### predenoise_{read} ######
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

        ###### derep_{read} ######
        derep = tar_target(
          derep,
          optimotu.pipeline::derepFastq(
            predenoise,
            verbose = TRUE,
            names = optimotu.pipeline::file_to_sample_key(predenoise)
          ),
          pattern = map(predenoise),
          resources = tar_resources(
            crew = tar_resources_crew(controller = "wide")
          )
        ),

        ###### err_{read} ######
        err = tar_target(
          err,
          optimotu.pipeline::learnErrors(
            purrr::discard(predenoise, grepl, pattern = "BLANK|NEG"),
            errorEstimationFunction = errfun,
            multithread = optimotu.pipeline::local_cpus(),
            verbose = TRUE
          ),
          resources = tar_resources(
            crew = tar_resources_crew(controller = "wide")
          )
        ),

        ###### denoise_{read} ######
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
          resources = tar_resources(
            crew = tar_resources_crew(controller = "wide")
          )
        )
      ),

      ##### merged #####
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
        resources = tar_resources(
          crew = tar_resources_crew(controller = "wide")
        )
      ),

      ##### denoise_map #####
      denoise_map = tar_fst_tbl(
        denoise_map,
        optimotu.pipeline::make_denoise_map(
          merged,
          seq_all,
          rc = .orient == "rev"
        ),
        pattern = map(merged),
        resources = tar_resources(
          crew = tar_resources_crew(controller = "thin")
        )
      ),

      ##### seqtable_raw #####
      seqtable_raw = tar_fst_tbl(
        seqtable_raw,
        optimotu.pipeline::denoise_map_to_seqtable(denoise_map),
        pattern = map(denoise_map),
        resources = tar_resources(
          crew = tar_resources_crew(controller = "thin")
        )
      ),

      ##### read_map #####
      read_map = tar_target(
        read_map,
        optimotu.pipeline::dada2_read_map(
          sample = samplewise_meta$sample_key,
          fq_raw = samplewise_meta$fastq_R1,
          fq_trim = samplewise_meta$trim_R1,
          fq_filt = samplewise_meta$filt_R1,
          dadaF = denoise_R1,
          derepF = derep_R1,
          dadaR = denoise_R2,
          derepR = derep_R2,
          merged = merged,
          denoise_map = denoise_map
        ),
        pattern = map(
          samplewise_meta,
          denoise_R1,
          derep_R1,
          denoise_R2,
          derep_R2,
          merged,
          denoise_map
        ),
        resources = tar_resources(
          crew = tar_resources_crew(controller = "wide")
        )
      )
    )
  },

  if (optimotu.pipeline::do_unoise()) {
    list(
      ##### samplewise_meta #####
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
          dplyr::semi_join(merge_read_counts, by = "merged") |>
          dplyr::select(
            seqrun,
            sample,
            readwise_key,
            sample_key,
            fastq_R1,
            trim_R1,
            tidyselect::all_of("merged"),
            to_denoise_merged,
            any_of(c(
              "rarefy_text",
              "numerator",
              "denominator",
              "number",
              "tar_seed"
            ))
          ),
        pattern = map(readwise_meta)
      ),

      ##### predenoise_merged #####
      predenoise_merged = if (optimotu.pipeline::do_rarefy()) {
        tar_file(
          predenoise_merged,
          mapply(
            optimotu.pipeline::fastq_sample,
            infile = samplewise_meta$merged,
            outfile = samplewise_meta$to_denoise_merged,
            n = !!(if (is.null(optimotu.pipeline::rarefy_number())) {
              quote(round(
                .numerator * merge_read_counts$merge_nread / .denominator
              ))
            } else {
              quote(.number)
            }),
            sample = mapply(
              \(n, seed) {
                tar_seed_set(seed)
                sample(n)
              },
              merge_read_counts$merge_nread,
              samplewise_meta$tar_seed,
              SIMPLIFY = FALSE
            ),
            SIMPLIFY = TRUE
          ),
          pattern = map(samplewise_meta, premerge_seqs, merge_read_counts),
          resources = tar_resources(
            crew = tar_resources_crew(controller = "thin")
          )
        )
      } else {
        tar_file(
          predenoise_merged,
          premerge_seqs,
          pattern = map(premerge_seqs),
          deployment = "main"
        )
      },

      ##### unoise #####
      unoise = tar_target(
        unoise,
        {
          files <- predenoise_merged
          names(files) <- samplewise_meta$sample_key
          optimotu.pipeline::vsearch_cluster_unoise2(
            files,
            min_size = !!optimotu.pipeline::unoise_minsize(),
            alpha = !!optimotu.pipeline::unoise_alpha(),
            threads = 1L,
            shards = optimotu.pipeline::local_cpus()
          ) |>
            stats::setNames(samplewise_meta$sample_key)
        },
        pattern = map(samplewise_meta, predenoise_merged),
        resources = tar_resources(
          crew = tar_resources_crew(controller = "wide")
        )
      ),

      ##### denoise_map #####
      denoise_map = tar_fst_tbl(
        denoise_map,
        optimotu.pipeline::make_denoise_map(
          unoise,
          seq_all,
          rc = .orient == "rev"
        ),
        pattern = map(unoise),
        resources = tar_resources(
          crew = tar_resources_crew(controller = "thin")
        )
      ),

      ##### seqtable_raw #####
      seqtable_raw = tar_fst_tbl(
        seqtable_raw,
        optimotu.pipeline::denoise_map_to_seqtable(denoise_map),
        pattern = map(denoise_map),
        resources = tar_resources(
          crew = tar_resources_crew(controller = "thin")
        )
      ),

      ##### read_map #####
      read_map = tar_target(
        read_map,
        optimotu.pipeline::unoise_read_map(
          sample = samplewise_meta$sample_key,
          fq_raw = samplewise_meta$fastq_R1,
          fq_trim = samplewise_meta$trim_R1,
          fq_merged = predenoise_merged,
          uc = unoise,
          denoise_map = denoise_map
        ),
        pattern = map(samplewise_meta, predenoise_merged, unoise, denoise_map),
        resources = tar_resources(
          crew = tar_resources_crew(controller = "wide")
        )
      )
    )
  }
)

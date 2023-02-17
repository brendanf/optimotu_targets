addedmodel_dir <- "protaxFungi/addedmodel"

if (file.exists("data/culture_refs.fasta") &&
    file.exists("data/culture_sequences.xlsx")) {
  restorationmodel_dir <- "restorationmodel"
  if (!dir.exists(restorationmodel_dir)) dir.create(restorationmodel_dir)
  
  refseq_plan <- list(
    tar_file(
      taxonomy_addedmodel_file,
      file.path(addedmodel_dir, "taxonomy")
    ),
    tar_fst_tbl(
      taxonomy_addedmodel,
      readr::read_tsv(
        taxonomy_addedmodel_file,
        col_names = c("taxon_id", "parent_id", "rank", "classification", "prior"),
        col_types = "iiicd"
      )
    ),
    tar_file(
      taxonomy_ascii7_addedmodel_file,
      file.path(addedmodel_dir, "taxonomy.ascii7")
    ),
    tar_fst_tbl(
      taxonomy_ascii7_addedmodel,
      readr::read_tsv(
        taxonomy_ascii7_addedmodel_file,
        col_names = c("taxon_id", "parent_id", "rank", "classification", "prior"),
        col_types = "iiicd"
      )
    ),
    tar_file(
      new_refseq_file,
      "data/culture_refs.fasta"
    ),
    tar_target(
      new_refseq,
      Biostrings::readDNAStringSet(new_refseq_file)
    ),
    tar_file(
      new_refseq_metadata_file,
      "data/culture_sequences.xlsx"
    ),
    tar_fst_tbl(
      new_refseq_metadata,
      readxl::read_excel(new_refseq_metadata_file) %>%
        dplyr::filter(Culture_ID %in% names(new_refseq))
    ),
    tar_fst_tbl(
      taxonomy_new,
      build_taxonomy(
        taxonomy_addedmodel$classification,
        new_refseq_metadata$Protax_synonym
      )
    ),
    tar_file(
      write_protax_taxonomy_new,
      write_and_return_file(
        taxonomy_new,
        file.path(restorationmodel_dir, "taxonomy"),
        "tsv",
        col_names = FALSE
      )
    ),
    tar_file(
      write_its2_new,
      {
        outfile <- file.path(restorationmodel_dir, "its2.fa")
        file.copy("protaxFungi/addedmodel/its2.fa", outfile, overwrite = TRUE)
        file.append(outfile, new_refseq_file)
        outfile
      }
    ),
    tar_file(
      write_sintaxits2_new,
      {
        outfile <- file.path(restorationmodel_dir, "sintaxits2train.fa")
        file.copy("protaxFungi/addedmodel/sintaxits2train.fa", outfile, overwrite = TRUE)
        tibble::enframe(as.character(new_refseq), name = "Culture_ID") %>%
          dplyr::left_join(
            dplyr::select(new_refseq_metadata, Culture_ID, Protax_synonym),
            by = "Culture_ID"
          ) %>%
          dplyr::transmute(
            name = paste(Culture_ID, sintax_format(Protax_synonym), sep = ";"),
            value = value
          ) %>%
          tibble::deframe() %>%
          Biostrings::DNAStringSet() %>%
          Biostrings::writeXStringSet(outfile, append = TRUE)
        outfile
      }
    ),
    tar_file(
      write_its2udb_new,
      build_udb(
        write_its2_new,
        file.path(restorationmodel_dir, "its2.udb"),
        type = "usearch",
        usearch = "protaxFungi/scripts/usearch10.0.240_i86linux32"
      )
    ),
    tar_file(
      write_sintaxits2udb_new,
      build_udb(
        write_sintaxits2_new,
        file.path(restorationmodel_dir, "sintaxits2.udb"),
        type = "sintax",
        usearch = "protaxFungi/scripts/usearch10.0.240_i86linux32"
      )
    ),
    tar_file(
      write_amptksynmockudb,
      {
        outfile <- file.path(restorationmodel_dir, "amptk_synmock.udb")
        file.symlink("../protaxFungi/addedmodel/amptk_synmock.udb", outfile)
        outfile
      }
    ),
    tar_file(
      write_protax_taxonomy.ascii7_new,
      write_and_return_file(
        dplyr::mutate(
          taxonomy_new,
          classification = ascii_clean(classification)
        ),
        file.path(restorationmodel_dir, "taxonomy.ascii7"),
        "tsv",
        col_names = FALSE
      )
    ),
    tar_map(
      values = list(.rank = 2:7),
      tar_file(
        write_protax_tax,
        write_and_return_file(
          dplyr::filter(taxonomy_new, rank <= .rank),
          file.path(restorationmodel_dir, paste0("tax", .rank)),
          "tsv",
          col_names = FALSE
        )
      ),
      tar_file(
        write_protax_ref.tax,
        {
          outfile <- paste0("ref.tax", .rank)
          file.copy(
            from = file.path(addedmodel_dir, outfile),
            to = file.path(restorationmodel_dir, outfile),
            overwrite = TRUE
          )
          dplyr::transmute(
            new_refseq_metadata,
            Culture_ID = Culture_ID,
            Protax_synonym = truncate_taxonomy(Protax_synonym, .rank)
          ) %>%
            write_and_return_file(
              file.path(restorationmodel_dir, outfile),
              "tsv",
              append = TRUE
            )
        }
      ),
      tar_file(
        write_protax_rseqs,
        dplyr::bind_rows(
          readr::read_tsv(
            file.path(addedmodel_dir, sprintf("rseqs%d", .rank)),
            col_names = c("taxon_id", "accno"),
            col_types = "ic"
          ) %>%
            tidyr::separate_rows(accno, sep = ","),
          dplyr::transmute(
            new_refseq_metadata,
            accno = Culture_ID,
            classification = truncate_taxonomy(Protax_synonym, .rank)
          ) %>%
            dplyr::left_join(
              dplyr::filter(taxonomy_new, rank == .rank),
              by = "classification"
            ) %>%
            dplyr::select(taxon_id, accno)
        ) %>%
          dplyr::group_by(taxon_id) %>%
          dplyr::summarize(accno = paste(accno, collapse = ",")) %>%
          write_and_return_file(
            file.path(restorationmodel_dir, sprintf("rseqs%d", .rank)),
            type = "tsv",
            col_names = FALSE
          )
      )
    )
  )
  
  restorationmodel_files <- purrr::keep(get_target_names(refseq_plan), startsWith, "write_")
  
  paste("{", paste(restorationmodel_files, collapse = "\n"), shQuote(restorationmodel_dir), "}", sep = "\n")
  
  refseq_plan <- c(
    refseq_plan,
    tar_target_raw(
      name = "restorationmodel",
      command = parse(
        text = paste(
          "{",
          paste(restorationmodel_files, collapse = "\n"),
          shQuote(restorationmodel_dir),
          "}",
          sep = "\n"
        )
      ),
      format = "file"
    ),
    tar_file(
      protax_model,
      restorationmodel
    )
  )
} else {
  refseq_plan <- list(
    tar_file(
      protax_model,
      addedmodel_dir
    )
  )
}

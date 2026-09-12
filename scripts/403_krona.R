krona_plan <- c(
  list(
    #### krona_script ####
    # character: KronaTools script for embedding in Krona plots
    krona_script = tar_target(
      krona_script,
      readLines(
        withr::local_connection(
          url("http://marbl.github.io/Krona/src/krona-2.0.js")
        )
      ),
      deployment = "main"
    ),

    #### krona_shortcut_icon ####
    # character: base64-encoded KronaTools shortcut icon for embedding in Krona
    #   plots
    krona_shortcut_icon = tar_target(
      krona_shortcut_icon,
      base64enc::dataURI(
        withr::local_connection(
          url("http://marbl.github.io/Krona//img/favicon.ico", open = "rb")
        ),
        mime = "image/x-icon"
      ),
      deployment = "main"
    ),

    #### krona_hiddenimage ####
    # character: base64-encoded KronaTools "hidden" icon for embedding in Krona
    #   plots
    krona_hiddenimage = tar_target(
      krona_hiddenimage,
      base64enc::dataURI(
        withr::local_connection(
          url("http://marbl.github.io/Krona//img/hidden.png", open = "rb")
        ),
        mime = "image/png"
      ),
      deployment = "main"
    ),

    #### krona_loadingimage ####
    # character: base64-encoded KronaTools "loading" icon for embedding in Krona
    #   plots
    krona_loadingimage = tar_target(
      krona_loadingimage,
      base64enc::dataURI(
        withr::local_connection(
          url("http://marbl.github.io/Krona//img/loading.gif", open = "rb")
        ),
        mime = "image/gif"
      ),
      deployment = "main"
    ),

    #### krona_logo ####
    # character: base64-encoded KronaTools logo for embedding in Krona plots
    krona_logo = tar_target(
      krona_logo,
      base64enc::dataURI(
        withr::local_connection(
          url("http://marbl.github.io/Krona//img/logo-small.png", open = "rb")
        ),
        mime = "image/png"
      ),
      deployment = "main"
    )
  ),

  tar_map(
    values = tibble::tibble(
      .conf_level = c("plausible", "reliable"),
      otu_taxonomy = paste0("otu_taxonomy_", .conf_level) |>
        rlang::syms()
    ),
    names = .conf_level,
    #### otu_krona_data_{.conf_level} ####
    # tibble from generate_krona_data(): rank, taxon, parent_taxonomy,
    # per-rank unknown fractions, nread/nocc/notu, child_unknown_*,
    # fread/focc/fotu
    tar_fst_tbl(
      otu_krona_data,
      optimotu.pipeline::generate_krona_data(
        otu_taxonomy,
        ranks = !!optimotu.pipeline::tax_ranks()
      ),
      deployment = "main"
    ),

    #### write_otu_krona_{.conf_level} ####
    # character (output filename)
    #
    # write a stand-alone HTML file containing the Krona plot
    tar_file(
      write_otu_krona,
      file.path(
        !!optimotu.pipeline::output_path(),
        !!(if (optimotu.pipeline::do_rarefy()) {
          quote(sprintf("otu_krona_%s_%s.html", .conf_level, .rarefy_text))
        } else {
          quote(sprintf("otu_krona_%s.html", .conf_level))
        })
      ) |>
        optimotu.pipeline::krona_xml_nodes(
          data = dplyr::filter(
            otu_krona_data,
            (nocc >= 5) | (notu >= 5) | (nread > 1000)
          ),
          .rank = !!optimotu.pipeline::root_rank(),
          maxrank = !!optimotu.pipeline::tip_rank(),
          outfile = _,
          node_data_format = !!optimotu.pipeline::krona_node_data_format(
            optimotu.pipeline::tax_ranks()
          ),
          taxonomy = NULL,
          pre = c(
            '<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">',
            '<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">',
            ' <head>',
            '  <meta charset="utf-8"/>',
            paste0(
              '  <link rel="shortcut icon" href="',
              krona_shortcut_icon,
              '"/>'
            ),
            '  <script id="notfound" type="text/javascript">window.onload=function(){document.body.innerHTML=""}</script>',
            '  <script language="javascript" type="text/javascript">',
            krona_script,
            '  </script>',
            ' </head>',
            ' <body>',
            paste0(
              '  <img id="hiddenImage" src="',
              krona_hiddenimage,
              '" style="display:none" alt="Hidden Image"/>'
            ),
            paste0(
              '  <img id="loadingImage" src="',
              krona_loadingimage,
              '" style="display:none" alt="Loading Indicator"/>'
            ),
            paste0(
              '  <img id="logo" src="',
              krona_logo,
              '" style="display:none" alt="Logo of Krona"/>'
            ),
            '  <noscript>Javascript must be enabled to view this page.</noscript>',
            '  <div style="display:none">',
            '<krona>',
            !!optimotu.pipeline::krona_html_attributes(),
            '<datasets>',
            '<dataset>Occurence weighting</dataset>',
            '<dataset>Read abundance weighting</dataset>',
            '<dataset>OTU richness weighting</dataset>',
            '</datasets>'
          ),
          post = "</krona>"
        ),
      deployment = "main"
    )
  )
)

optimotu_plan <- c(optimotu_plan, krona_plan)

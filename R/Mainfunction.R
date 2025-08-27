#' Data information for reference samples
#'
#' \code{DataInfobeta3Dmeta} provides basic data information for each combination of site and treatment for
#' (1) the gamma reference sample in the pooled assemblage, and
#' (2) the alpha reference sample in the joint assemblage.
#'
#' @param data (a) For \code{datatype = "abundance"}, a data.frame whose first two columns are study/site and treatment,
#' followed by species columns (an assemblage = site × treatment). \cr
#' (b) For \code{datatype = "incidence_raw"}, a data.frame whose first three columns are study/site, treatment, and patch,
#' followed by species columns with entries in \{0,1\}.
#' @param diversity Diversity type: 'TD' (Taxonomic), 'PD' (Phylogenetic), or 'FD' (Functional).
#' @param datatype "abundance" for individual-based data, or "incidence_raw" for incidence-by-sampling-units data (0/1).
#' @param PDtree (required for \code{diversity = "PD"}) a phylogenetic tree in Newick format for all observed species in the pooled data.
#' @param PDreftime (required only for \code{diversity = "PD"}) a numeric reference time; default \code{NULL} uses the root age of \code{PDtree}.
#' @param FDdistM (required for \code{diversity = "FD"}) a species pairwise distance matrix for the pooled data.
#' @param FDtype (required for \code{diversity = "FD"}) either \code{"tau_value"} for a specified threshold or \code{"AUC"} (default) for
#' the area under the tau-profile.
#' @param FDtau (required for \code{diversity = "FD"} and \code{FDtype = "tau_value"}) a value in [0,1] specifying the threshold.
#' If \code{FDtau = NULL} (default), the threshold is set to the mean pairwise distance (quadratic entropy).
#'
#' @importFrom dplyr group_by summarise filter select across pick everything ungroup reframe
#' @importFrom tidyr unnest
#' @importFrom iNEXT.beta3D DataInfobeta3D
#' @importFrom magrittr %>%
#' @importFrom stringr str_to_title
#' @importFrom utils globalVariables
#'
#' @return A data.frame with basic information restricted to the pooled and joint assemblages. Columns follow
#'         \code{iNEXT.beta3D::DataInfobeta3D} output.
#'
#' @examples
#' ## (Data Information) Taxonomic diversity for abundance data
#' data("Beetle_abundance_data")
#' info_TD_abu <- DataInfobeta3Dmeta(data = Beetle_abundance_data,
#'                                   diversity = 'TD',
#'                                   datatype = 'abundance')
#' info_TD_abu
#'
#' ## (Data Information) Taxonomic diversity for incidence data (Bird)
#' data("Bird_incidence_data")
#' info_TD_inc <- DataInfobeta3Dmeta(data = Bird_incidence_data,
#'                                   diversity = 'TD',
#'                                   datatype = 'incidence_raw')
#' info_TD_inc
#'
#' ## (Data Information) Mean phylogenetic diversity for abundance data
#' data("Beetle_abundance_data")
#' data("Beetle_tree")
#' info_PD_abu <- DataInfobeta3Dmeta(data = Beetle_abundance_data,
#'                                   diversity = 'PD',
#'                                   datatype = 'abundance',
#'                                   PDtree = Beetle_tree, PDreftime = NULL)
#' info_PD_abu
#'
#' ## (Data Information) Mean phylogenetic diversity for incidence data (Bird)
#' data("Bird_incidence_data")
#' data("Bird_tree")
#' info_PD_inc <- DataInfobeta3Dmeta(data = Bird_incidence_data,
#'                                   diversity = 'PD',
#'                                   datatype = 'incidence_raw',
#'                                   PDtree = Bird_tree,
#'                                   PDreftime = NULL)
#' info_PD_inc
#'
#' ## (Data Information) Functional diversity for abundance data (tau-value)
#' data("Beetle_abundance_data")
#' data("Beetle_distM")
#' info_FDtau_abu <- DataInfobeta3Dmeta(data = Beetle_abundance_data,
#'                                      diversity = 'FD',
#'                                      datatype = 'abundance',
#'                                      FDdistM = Beetle_distM,
#'                                      FDtype = "tau_value", FDtau = NULL)
#' info_FDtau_abu
#'
#' ## (Data Information) Functional diversity for abundance data (AUC)
#' data("Beetle_abundance_data")
#' data("Beetle_distM")
#' info_FDAUC_abu <- DataInfobeta3Dmeta(data = Beetle_abundance_data,
#'                                      diversity = 'FD',
#'                                      datatype = 'abundance',
#'                                      FDdistM = Beetle_distM,
#'                                      FDtype = "AUC", FDtau = NULL)
#' info_FDAUC_abu
#'
#' ## (Data Information) Functional diversity for incidence data (Bird)
#' data("Bird_incidence_data")
#' data("Bird_distM")
#' info_FDtau_inc <- DataInfobeta3Dmeta(data = Bird_incidence_data,
#'                                      diversity = 'FD',
#'                                      datatype = 'incidence_raw',
#'                                      FDdistM = Bird_distM,
#'                                      FDtype = "tau_value", FDtau = NULL)
#' info_FDtau_inc
#'
#' @export
DataInfobeta3Dmeta <- function(data, diversity = "TD", datatype = "abundance",
                               PDtree = NULL, PDreftime = NULL,
                               FDdistM = NULL, FDtype = "AUC", FDtau = NULL){
  if (datatype == "abundance"){
    data |>
      group_by(across(1:2)) |>
      summarise(mat = list(t(pick(everything())) %>%
                             DataInfobeta3D(data = ., diversity, datatype = "abundance",
                                            PDtree, PDreftime, FDdistM, FDtype, FDtau) |>
                             filter(Assemblage %in% c("Pooled assemblage", "Joint assemblage"))),
                .groups = "drop") |>
      unnest(mat) |>
      select(-Dataset)
  } else if (datatype == "incidence_raw") {
    data |>
      group_by(across(1:3)) |>
      summarise(mat = list(t(pick(everything()))), .groups = "drop") %>%
      filter(!sapply(mat, function(m) all(as.matrix(m) == 0))) |>
      group_by(across(1:2)) |>
      summarise(mat = list(mat), .groups = "keep") %>%
      reframe(mat %>%
                DataInfobeta3D(data = ., diversity, datatype = "incidence_raw",
                               PDtree, PDreftime, FDdistM, FDtype, FDtau) |>
                filter(Assemblage %in% c("Pooled assemblage", "Joint assemblage"))) |>
      select(-Dataset)
  }
}

#' Difference of standardized 3D diversity (fixed-effect meta-analysis)
#'
#' \code{iNEXTbeta3Dmeta} estimates the difference of standardized 3D (taxonomic, phylogenetic, functional)
#' diversity between two treatments (e.g., Enhanced vs. Control) and synthesizes results across sites using a
#' fixed-effect model.
#'
#' @param data Same structure as \code{iNEXT.beta3D}; for \code{datatype = "abundance"}, the first two columns are site and treatment;
#' for \code{datatype = "incidence_raw"}, the first three columns are site, treatment, and patch.
#' @param diversity 'TD', 'PD', or 'FD'.
#' @param order.q Numeric vector of diversity orders (e.g., \code{c(0,1,2)}).
#' @param datatype "abundance" or "incidence_raw".
#' @param level Sample-coverage level in [0,1]. If \code{NULL}, the function uses the minimum coverage achieved at
#' double the reference sample size across all site × treatment combinations.
#' @param nboot Number of bootstrap replications for uncertainty.
#' @param treatment_order Character vector of length 2 specifying the treatment order; differences are computed as
#' \code{treatment_order[1] - treatment_order[2]}.
#' @param conf Confidence level (default 0.95).
#' @param PDtree,PDreftime,PDtype,FDdistM,FDtype,FDtau,FDcut_number As in \code{iNEXT.beta3D}.
#'
#' @import iNEXT.3D
#' @importFrom dplyr group_by group_split summarise filter select mutate ungroup rename rename_with bind_rows across everything
#' @importFrom purrr map map2_dfr map_dfr map_df imap
#' @importFrom stringr str_to_title
#' @importFrom tibble tibble
#' @importFrom stats qnorm
#' @importFrom metafor rma
#' @importFrom iNEXT.beta3D iNEXTbeta3D
#'
#' @return A list with seven components: \code{Summary_Gamma}, \code{Summary_Alpha}, \code{Summary_Beta},
#' \code{Summary_1-C}, \code{Summary_1-U}, \code{Summary_1-V}, and \code{Summary_1-S}. Each component contains:
#' (1) a data frame of site-level differences (with SE, CI, weights) plus a fixed-effect combined row, and
#' (2) a one-row meta summary with Cochran’s Q (\code{Q_val}), degrees of freedom (\code{df_val}),
#' p-value (\code{p_val}), and heterogeneity percentage (\code{I2_val}).
#'
#' @examples
#' ## TD / abundance (Beetle), coverage-standardized, fixed effect
#' data("Beetle_abundance_data")
#' out_TD_abu <- iNEXTbeta3Dmeta(data = Beetle_abundance_data, diversity = "TD",
#'                               order.q = c(0,1,2), datatype = "abundance", level = NULL,
#'                               nboot = 10, treatment_order = c("Enhanced","Control"), conf = 0.95)
#' out_TD_abu
#'
#' ## PD / abundance (Beetle) — requires Beetle_tree
#' data("Beetle_abundance_data"); data("Beetle_tree")
#' out_PD_abu <- iNEXTbeta3Dmeta(data = Beetle_abundance_data, diversity = "PD",
#'                               order.q = c(0,1,2), datatype = "abundance", level = NULL,
#'                               nboot = 10, treatment_order = c("Enhanced","Control"), conf = 0.95,
#'                               PDtree = Beetle_tree, PDreftime = NULL, PDtype = "meanPD")
#' out_PD_abu
#'
#' ## FD / abundance (Beetle) — requires Beetle_distM
#' data("Beetle_abundance_data"); data("Beetle_distM")
#' out_FD_abu <- iNEXTbeta3Dmeta(data = Beetle_abundance_data, diversity = "FD",
#'                               order.q = c(0,1,2), datatype = "abundance", level = NULL,
#'                               nboot = 10, treatment_order = c("Enhanced","Control"), conf = 0.95,
#'                               FDdistM = Beetle_distM, FDtype = "AUC", FDcut_number = 30)
#' out_FD_abu
#'
#' ## TD / incidence_raw (Bird)
#' data("Bird_incidence_data")
#' out_TD_inc <- iNEXTbeta3Dmeta(data = Bird_incidence_data, diversity = "TD",
#'                               order.q = c(0,1,2), datatype = "incidence_raw", level = NULL,
#'                               nboot = 10, treatment_order = c("Enhanced","Control"), conf = 0.95)
#' out_TD_inc
#'
#' @export
iNEXTbeta3Dmeta <- function(data, diversity = "TD", order.q = 0, datatype = "abundance",
                            level = NULL, nboot = 10, treatment_order, conf = 0.95,
                            PDtree, PDreftime = NULL, PDtype = "meanPD",
                            FDdistM, FDtype = "AUC", FDtau = NULL, FDcut_number = 30){
  data[, 1] <- factor(data[, 1])

  # ---- decide Cmin (auto) ----
  if (is.null(level)) {
    Cmin <- if (datatype == "abundance") {
      data |>
        group_by(across(c(1, 2))) |>
        group_split() |>
        lapply(function(x){
          x |> select(-c(1, 2)) |> unlist() |>
            iNEXT.3D:::Coverage(., "abundance", 2 * sum(.))
        }) |>
        unlist() |> min()
    } else {
      data |>
        group_by(across(1:3)) |>
        summarise(mat = list(t(pick(everything()))), .groups = "drop") |>
        group_by(across(1:2)) |>
        summarise(mat = list(do.call(rbind, mat)), .groups = "drop") |>
        # .$mat |>
        pull(mat) |>
        lapply(function(x) iNEXT.3D:::Coverage(x, "incidence_raw", 2 * ncol(x))) |>
        unlist() |> min()
    }
  } else {
    Cmin <- level
  }

  # ---- split data by treatment & site ----
  if (datatype == "abundance") {
    data_T1 <- data %>% filter(.[[2]] == treatment_order[1]) %>% split(.[[1]]) %>%
      lapply(function(x) x |> select(-c(1, 2)) |> t())
    data_T2 <- data %>% filter(.[[2]] == treatment_order[2]) %>% split(.[[1]]) %>%
      lapply(function(x) x |> select(-c(1, 2)) |> t())
  } else { # incidence_raw
    data_T1 <- data %>% filter(.[[2]] == treatment_order[1]) %>% split(.[[1]]) %>%
      map(~ split(.x, .[[3]]) %>% map(~ .x |> select(-c(1,2,3)) |> t()))
    data_T2 <- data %>% filter(.[[2]] == treatment_order[2]) %>% split(.[[1]]) %>%
      map(~ split(.x, .[[3]]) %>% map(~ .x |> select(-c(1,2,3)) |> t()))
  }

  z <- qnorm((1 - conf) / 2, lower.tail = FALSE)

  # ---- compute standardized diversity for two treatments ----
  ibeta_div_T1 <- iNEXTbeta3D(data_T1, diversity = diversity, q = order.q, datatype = datatype, level = Cmin, nboot = nboot,
                              PDtree = PDtree, PDreftime = PDreftime, PDtype = PDtype,
                              FDdistM = FDdistM, FDtype = FDtype, FDtau = FDtau, FDcut_number = FDcut_number)
  ibeta_div_T2 <- iNEXTbeta3D(data_T2, diversity = diversity, q = order.q, datatype = datatype, level = Cmin, nboot = nboot,
                              PDtree = PDtree, PDreftime = PDreftime, PDtype = PDtype,
                              FDdistM = FDdistM, FDtype = FDtype, FDtau = FDtau, FDcut_number = FDcut_number)

  tidy_fun <- function(x){
    x |>
      map(~ imap(.x, ~ mutate(.x, Type = str_to_title(.y)) |>
                   rename_with(~ ifelse(. %in% c("Gamma","Alpha","Beta","Dissimilarity"), "qD", .)))) |>
      map_dfr(~ bind_rows(.x))
  }

  gab_T1 <- tidy_fun(ibeta_div_T1)
  gab_T2 <- tidy_fun(ibeta_div_T2)

  # ---- combine, compute FE weights & meta summary ----
  one_component <- function(kind) {
    d1 <- gab_T1 |> filter(Type == kind)
    d2 <- gab_T2 |> filter(Type == kind)

    merged <- merge(d1, d2, by = c("Dataset","Order.q"), suffixes = c("_T1","_T2")) |>
      mutate(yi = qD_T1 - qD_T2, vi = s.e._T1^2 + s.e._T2^2,
             SE = sqrt(vi),
             LCL = yi - z * SE, UCL = yi + z * SE) |>
      rename(Difference = yi) |>
      group_split(Order.q)

    # fixed-effect meta via metafor::rma(..., method = "FE")
    models <- merged |> map(~ rma(yi = Difference, sei = SE, method = "FE", data = .x))

    meta_note <- map2_dfr(models, order.q, \(m, qv){
      tibble(Order.q = qv,
             Q_val = round(m$QE, 2),
             df_val = m$k - 1,
             p_val = ifelse(m$QEp < 0.001, "< 0.001", format(round(m$QEp, 3), nsmall = 3)),
             I2_val = round(m$I2, 1))
    })

    meta_rows <- lapply(seq_along(models), function(i){
      m <- models[[i]]
      m$data |>
        mutate(Weight = ((1/vi) / sum(1/vi)) * 100) |>
        select(Dataset, Difference, SE, LCL, UCL, Order.q, Diversity_T1, qD_T1, qD_T2, Weight) |>
        bind_rows(
          tibble(Dataset   = "Fixed-effect model",
                 Difference = as.numeric(m$beta),
                 SE         = m$se,
                 LCL        = m$ci.lb,
                 UCL        = m$ci.ub,
                 Order.q    = unique(m$data$Order.q),
                 Diversity_T1 = diversity,
                 qD_T1 = NA_real_, qD_T2 = NA_real_,
                 Weight = 100)
        )
    }) |> bind_rows() |>
      rename(Site = Dataset, Diversity = Diversity_T1) |>
      rename_with(~ treatment_order, .cols = c(qD_T1, qD_T2)) |>
      as.data.frame()

    list(meta_rows, meta_note)
  }

  out_list <- lapply(c("Gamma","Alpha","Beta","1-C","1-U","1-V","1-S"), one_component)
  names(out_list) <- paste("Summary", c("Gamma","Alpha","Beta","1-C","1-U","1-V","1-S"), sep = "_")
  out_list
}

#' Forest plot of the standardized 3D diversity difference (fixed effect)
#'
#' \code{ggiNEXTmeta} draws a forest plot for the difference of standardized 3D diversity between two treatments,
#' summarized using a fixed-effect model.
#'
#' @param output The output from \code{iNEXTbeta3Dmeta}.
#' @param q A value of \code{Order.q} present in \code{output}.
#' @param num_round Number of decimal places shown.
#' @param range Optional lower/upper bounds to clip CIs (pass to \code{forestplot}'s \code{clip}; commented out by default).
#' @param type One of "Gamma", "Alpha", "Beta", "1-C", "1-U", "1-V", or "1-S".
#' @param level Optional sample coverage (in percent) to annotate in the title; if \code{NULL}, no annotation is added.
#'
#' @importFrom dplyr filter
#' @importFrom forestplot forestplot fp_set_style fp_add_header fp_append_row fp_decorate_graph fp_set_zebra_style fp_add_lines fpTxtGp
#' @importFrom tibble tibble
#' @importFrom magrittr %>%
#' @importFrom grid grid.text gpar
#' @importFrom gridExtra grid.arrange
#'
#' @return A forest plot showing site-level differences and the fixed-effect combined estimate.
#'
#' @examples
#' data("Beetle_abundance_data")
#' out <- iNEXTbeta3Dmeta(Beetle_abundance_data, diversity = "TD",
#'                        order.q = c(0,1,2), datatype = "abundance",
#'                        level = NULL, nboot = 10,
#'                        treatment_order = c("Enhanced","Control"))
#' ggiNEXTmeta(out, q = 0, num_round = 3, range = c(-20,15), type = "Gamma", level = NULL)
#' @export
ggiNEXTmeta <- function(output, q, num_round = 3, range, type = NULL, level = NULL){
  if (!is.null(type)) output <- output[paste("Summary", type, sep = "_")]

  meta_output <- output[[1]][[1]]
  meta_note   <- output[[1]][[2]]

  meta_output <- meta_output |> filter(Order.q == q)
  meta_note   <- meta_note   |> filter(Order.q == q)

  meta_note_q <- paste0(
    meta_output[nrow(meta_output), 1], " (Q = ", meta_note$Q_val,
    ", df = ", meta_note$df_val,
    ", p = ", meta_note$p_val,
    "; I2 = ", meta_note$I2_val, ")"
  )

  name        <- colnames(meta_output)[1]
  name_treat  <- colnames(meta_output)[c(8, 9)]
  ll          <- nrow(meta_output)
  name_site   <- meta_output[, 1][-ll]
  diversity   <- meta_output$Diversity[1]
  order       <- meta_output$Order.q[1]

  forestplot_q <- tibble(
    mean   = round(meta_output$Difference[-ll], num_round),
    lower  = round(meta_output$LCL[-ll], num_round),
    upper  = round(meta_output$UCL[-ll], num_round),
    study  = name_site,
    q_T1   = round(meta_output[, 8][-ll], num_round),
    q_T2   = round(meta_output[, 9][-ll], num_round),
    diff   = round(meta_output$Difference[-ll], num_round),
    LCL    = round(meta_output$LCL[-ll], num_round),
    UCL    = round(meta_output$UCL[-ll], num_round),
    w_fixed= paste(round(meta_output$Weight[-ll], 2), "%", sep = "")
  )

  fig <- forestplot_q |>
    forestplot(labeltext = c(study, q_T1, q_T2, diff, LCL, UCL, w_fixed),
               title = if (is.null(level)) type else paste0(type, " (SC = ", level, "%)"),
               # clip = range,
               xlog = FALSE) |>
    fp_set_style(txt_gp = fpTxtGp(ticks = gpar(cex = 1))) |>
    fp_add_header(
      study = c("", name),
      q_T1  = c(paste(diversity, " (q = ", order, ")", sep = ""), name_treat[1]),
      q_T2  = c(paste(diversity, " (q = ", order, ")", sep = ""), name_treat[2]),
      diff  = c("", "Difference"),
      LCL   = c("", "LCL"),
      UCL   = c("", "UCL"),
      w_fixed = c("", "Weight")
    ) |>
    fp_append_row(
      mean  = round(meta_output$Difference[ll], num_round),
      lower = round(meta_output$LCL[ll], num_round),
      upper = round(meta_output$UCL[ll], num_round),
      study = "Fixed-effect model",
      diff  = round(meta_output$Difference[ll], num_round),
      LCL   = round(meta_output$LCL[ll], num_round),
      UCL   = round(meta_output$UCL[ll], num_round),
      is.summary = TRUE
    ) |>
    fp_decorate_graph(graph.pos = 4) |>
    fp_set_zebra_style("#EFEFEF") |>
    fp_add_lines()

  print(fig)
  grid.text(meta_note_q, x = 0.5, y = 0, just = "bottom", check.overlap = TRUE)
}

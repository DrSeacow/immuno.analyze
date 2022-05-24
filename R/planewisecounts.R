#'
#'
#' Summarize Immunohistochemistry Data with Graph and Unpaired-Samples T-Tests (Cell Counts, differentiated by coronal plane)
#'
#' @description Processes .CSV files output by the "cell outline" batch macro in the ImageJ custom macros plugin (Timothy and Forlano, 2019).
#'
#' "Cell counts" data refers to a tally of rows, each denoting an above-threshold cell marked by a particular wavelength.
#'
#' Can only be used to analyze one wavelength per analysis. This R function identifies two groups, produces graphs, and runs an unpaired-samples T-Test.
#'
#' It is generally recommended that images are taken and analyzed from multiple coronal planes of brain tissue. This particular function differentiates between the anterior, middle, and posterior planes and analyzes them separately, using 'patchwork' to display three barplots and vectorizing the P-values from the T-Tests as the anterior comparison P-Value, middle comparison P-Value, and posterior comparison P-Value, respectively.
#'
#' All parameters must be in quotations, besides "n_per_group." The .CSV files to be analyzed must be placed in your working directory.
#'
#' When imaging your tissue initially, it is critical that the output filenames are saved in the following format: GroupName_ExperimentName_SubjectNumber_MagnificationX_CoronalPlane_SideOfTissue_RegionOfTissue
#'
#' To mitigate experimenter bias (throughout all stages of the experimental process, but including this one) it is recommended that the group names are codified in a way agnostic to experimental treatment.
#'
#' @param graph_name_1 The title of the output barplot for the analysis of your "anterior" sections.
#' @param graph_name_2 The title of the output barplot for the analysis of your "middle" sections.
#' @param graph_name_3 The title of the output barplot for the analysis of your "posterior" sections.
#' @param DV The dependent variable used for the analysis. Should denote the cell type, tissue type or region, and/or wavelength analyzed.
#' @param first_group How one of your groups is codified. Note: When imaging tissue for processing, this string must be put at the beginning of the saved filenames.
#' @param second_group How the other group is codified. Note: "first_group" and "second_group" must be alphabetized in the respective order of your codified group names.
#' @param n_per_group How many subjects are included in each group. Used to compute standard error of the mean for error bars on barplot.
#' @param x The name of your .CSV file.
#' @return The analysis table for a two-sample T-Test comparing the groups, as well as a rudimentary barplot.
#' @examples cell.counts.plane("VTA_Counts", "med", "pos", "Cell_Counts", "Ctrl", "Switch", 6, data = "outline_VTA.csv")
#'
#' @references Timothy, M., & Forlano, P. M. (2019). A versatile macro-based neurohistological image analysis suite for ImageJ focused on automated and standardized user interaction and reproducible data output. Journal of neuroscience methods, 324, 108286. https://doi.org/10.1016/j.jneumeth.2019.04.009


cell.counts.plane <- function(graph_name_1, graph_name_2, graph_name_3, DV, first_group, second_group, n_per_group, data = x){

  library(stringr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)

  immunot_data <- as.data.frame(read.csv(data, header = TRUE, sep = ",", dec = ".")) %>%
    select(c("Label", "Mean", "X.Area"))
  immunot_data <- separate(immunot_data, "Label", c("Group", NA, "Subject", NA, "Plane", "Side", NA), "_")

  counts_by_plane <- tibble(aggregate(immunot_data$Plane, FUN = length, by = list(Plane = immunot_data$Plane, Group = immunot_data$Group, immunot_data$Subject)))

  gp1_ant <- counts_by_plane[counts_by_plane$Group == first_group & counts_by_plane$Plane == "Ant", ]
  aa <- as.vector(t(gp1_ant %>% select(x)))
  gp2_ant <- counts_by_plane[counts_by_plane$Group == second_group & counts_by_plane$Plane == "Ant", ]
  ab <- as.vector(t(gp2_ant %>% select(x)))

  gp1_med <- counts_by_plane[counts_by_plane$Group == first_group & counts_by_plane$Plane == "Med", ]
  ba <- as.vector(t(gp1_med %>% select(x)))
  gp2_med <- counts_by_plane[counts_by_plane$Group == second_group & counts_by_plane$Plane == "Med", ]
  bb <- as.vector(t(gp2_med %>% select(x)))

  gp1_pos <- counts_by_plane[counts_by_plane$Group == first_group & counts_by_plane$Plane == "Pos", ]
  ca <- as.vector(t(gp1_pos %>% select(x)))
  gp2_pos <- counts_by_plane[counts_by_plane$Group == second_group & counts_by_plane$Plane == "Pos", ]
  cb <- as.vector(t(gp2_pos %>% select(x)))

  plane_count <- tibble(Group = rep(rep(c("gp1", "gp2"), each = 6), 3),
                            Plane = rep(c("ant", "med", "pos"), each = 12),
                            Counts = c(aa, ab, ba, bb, ca, cb)
  )

  means_by_plane <- aggregate(counts_by_plane$x, FUN = mean, by = list(Plane = counts_by_plane$Plane, Group = counts_by_plane$Group))

  sds_by_plane <- tibble(aggregate(counts_by_plane$x, FUN = sd, by = list(Plane = counts_by_plane$Plane, Group = counts_by_plane$Group)))

  gp1_mean_ant <- means_by_plane[1,3]
  gp1_mean_med <- means_by_plane[2,3]
  gp1_mean_pos <- means_by_plane[3,3]
  gp2_mean_ant <- means_by_plane[4,3]
  gp2_mean_med <- means_by_plane[5,3]
  gp2_mean_pos <- means_by_plane[6,3]

  gp1_sterr_ant <- as.numeric(sds_by_plane[1,3]) / sqrt(length(n_per_group))
  gp1_sterr_med <- as.numeric(sds_by_plane[2,3]) / sqrt(length(n_per_group))
  gp1_sterr_pos <- as.numeric(sds_by_plane[3,3]) / sqrt(length(n_per_group))
  gp2_sterr_ant <- as.numeric(sds_by_plane[4,3]) / sqrt(length(n_per_group))
  gp2_sterr_med <- as.numeric(sds_by_plane[5,3]) / sqrt(length(n_per_group))
  gp2_sterr_pos <- as.numeric(sds_by_plane[6,3]) / sqrt(length(n_per_group))

  ant_compare <- c(gp1_mean_ant, gp2_mean_ant)

  med_compare <- c(gp1_mean_med, gp2_mean_med)

  pos_compare <- c(gp1_mean_pos, gp2_mean_pos)

  ant_sterrs <- as.numeric(c(gp1_sterr_ant, gp2_sterr_ant))

  med_sterrs <- as.numeric(c(gp1_sterr_med, gp2_sterr_med))

  pos_sterrs <- as.numeric(c(gp1_sterr_pos, gp2_sterr_pos))

  graph_df_ant <- tibble(Ant = rep(c(first_group, second_group)),
                         DV = ant_compare)

  graph_df_med <- tibble(Med = rep(c(first_group, second_group)),
                         DV = med_compare)

  graph_df_pos <- tibble(Pos = rep(c(first_group, second_group)),
                         DV = pos_compare)

  a <- ggplot(graph_df_ant, aes(x = Ant, y = ant_compare)) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin = (ant_compare - ant_sterrs), ymax = (ant_compare + ant_sterrs))) +
    ggtitle(graph_name_1) +
    ylab(DV)

  m <- ggplot(graph_df_med, aes(x = Med, y = DV)) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin = (med_compare - med_sterrs), ymax = (med_compare + med_sterrs))) +
    ggtitle(graph_name_2) +
    ylab(DV)

  p <- ggplot(graph_df_pos, aes(x = Pos, y = DV)) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin = (pos_compare - pos_sterrs), ymax = (pos_compare + pos_sterrs))) +
    ggtitle(graph_name_3) +
    ylab(DV)

  ant_t <- t.test(aa, ab, var.equal = TRUE)
  med_t <- t.test(ba, bb, var.equal = TRUE)
  pos_t <- t.test(ca, cb, var.equal = TRUE)

  analysis_list <- c(ant_t$p.value, med_t$p.value, pos_t$p.value)

  plot <- a+m+p

  return(list(analysis_list, plot))

}

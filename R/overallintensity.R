#'
#'
#' Summarize Immunohistochemistry Data with Graph and Unpaired-Samples T-Tests (Cell fluorescence intensity, not differentiated by coronal plane)
#'
#' @description Processes .CSV files output by the "cell outline" batch macro in the ImageJ custom macros plugin (Timothy and Forlano, 2019).
#'
#' "Intensity" data refers to the average of the ImageJ "mean grey value" for all above-threshold cells on a particular image, marked by a particular wavelength.
#'
#' Can only be used to analyze one wavelength per analysis. This R function identifies two groups, produces graphs, and runs an unpaired-samples T-Test.
#'
#' It is generally recommended that images are taken and analyzed from multiple coronal planes of brain tissue. This particular function averages across the anterior, middle, and posterior planes.
#'
#' All parameters must be in quotations, besides "n_per_group." The .CSV files to be analyzed must be placed in your working directory.
#'
#' When imaging your tissue initially, it is critical that the output filenames are saved in the following format: GroupName_ExperimentName_SubjectNumber_MagnificationX_CoronalPlane_LateralSideOfTissue_RegionOfTissue
#'
#' To mitigate experimenter bias (throughout all stages of the experimental process, but including this one) it is recommended that the group names are codified in a way agnostic to experimental treatment.
#'
#' @param graph_name The title of the output barplot.
#' @param DV The dependent variable used for the analysis. Should denote the cell type, tissue type or region, and/or wavelength analyzed.
#' @param first_group How one of your groups is codified. Note: When imaging tissue for processing, this string must be put at the beginning of the saved filenames.
#' @param second_group How the other group is codified. Note: "first_group" and "second_group" must be alphabetized in the respective order of your codified group names.
#' @param n_per_group How many subjects are included in each group. Used to compute standard error of the mean for error bars on barplot.
#' @param x The name of your .CSV file.
#' @return The analysis table for a two-sample T-Test comparing the groups, as well as a rudimentary barplot.
#' @examples overall.intensity("VTA_Intensity", "Intensity", "Ctrl", "Switch", 6, data = "outline_VTA.csv")
#'
#' @references Timothy, M., & Forlano, P. M. (2019). A versatile macro-based neurohistological image analysis suite for ImageJ focused on automated and standardized user interaction and reproducible data output. Journal of neuroscience methods, 324, 108286. https://doi.org/10.1016/j.jneumeth.2019.04.009


overall.intensity <- function(graph_name, DV, first_group, second_group, n_per_group, data = x){

  immunot_data <- as.data.frame(read.csv(data, header = TRUE, sep = ",", dec = ".")) %>%
    select(c("Label", "Mean", "X.Area"))
  immunot_data <- separate(immunot_data, "Label", c("Group", NA, "Subject", NA, "Plane", "Side", NA), "_")

  gp1_grand_mean <- mean(immunot_data[immunot_data$Group == first_group, 'Mean'])
  gp2_grand_mean <- mean(immunot_data[immunot_data$Group == second_group, 'Mean'])

  gp1_grand_sterr <- sd(immunot_data[immunot_data$Group == first_group, 'Mean']) / sqrt(n_per_group)
  gp2_grand_sterr <- sd(immunot_data[immunot_data$Group == second_group, 'Mean']) / sqrt(n_per_group)

  cross_plane_sterr <- c(gp1_grand_sterr, gp2_grand_sterr)

  cross_plane_compare <- as.numeric(c(gp1_grand_mean, gp2_grand_mean))

  graph_df_without_plane <- tibble(Group = c(first_group, second_group),
                                   DV = cross_plane_compare)

  plot <- ggplot(graph_df_without_plane, aes(x = Group, y = cross_plane_compare)) +
          geom_bar(stat = "identity") +
          geom_errorbar(aes(ymin = (cross_plane_compare - cross_plane_sterr), ymax = (cross_plane_compare + cross_plane_sterr))) +
          ggtitle(graph_name) +
          ylab(DV)

  values_by_group <- tibble(aggregate(immunot_data$Mean, FUN = mean, by = list(Group = immunot_data$Group, Subject = immunot_data$Subject)))

  values_by_group_avg_by_plane <- tibble(aggregate(values_by_group$x, FUN = mean, by = list(Group = values_by_group$Group, Subject = values_by_group$Subject)))

  gp1_values <- as.vector(values_by_group$x[values_by_group$Group == first_group])
  gp2_values <- as.vector(values_by_group$x[values_by_group$Group == second_group])

  overall_t <- t.test(gp1_values, gp2_values, var.equal = TRUE)
  return(list(overall_t, plot))

}

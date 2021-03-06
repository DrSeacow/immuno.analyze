% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/overallpercentarea.R
\name{overall.percent.area}
\alias{overall.percent.area}
\title{Summarize Immunohistochemistry Data with Graph and Unpaired-Samples T-Tests (Percent area, not differentiated by coronal plane)}
\usage{
overall.percent.area(
  graph_name,
  DV,
  first_group,
  second_group,
  n_per_group,
  data = x
)
}
\arguments{
\item{graph_name}{The title of the output barplot.}

\item{DV}{The dependent variable used for the analysis. Should denote the cell type, tissue type or region, and/or wavelength analyzed.}

\item{first_group}{How one of your groups is codified. Note: When imaging tissue for processing, this string must be put at the beginning of the saved filenames.}

\item{second_group}{How the other group is codified. Note: "first_group" and "second_group" must be alphabetized in the respective order of your codified group names.}

\item{n_per_group}{How many subjects are included in each group. Used to compute standard error of the mean for error bars on barplot.}

\item{x}{The name of your .CSV file.}
}
\value{
The analysis table for a two-sample T-Test comparing the groups, as well as a rudimentary barplot.
}
\description{
Processes .CSV files output by the "percent-area" batch macro in the ImageJ custom macros plugin (Timothy and Forlano, 2019).

"Percent area" data refers to the proportion of a selected region on an image that contains above-threshold cell activity marked by a particular wavelength.

Can be used for single-channel or fluorophor colocalization analyses. This R function identifies two groups, produces graphs, and runs an unpaired-samples T-Test.

It is generally recommended that images are taken and analyzed from multiple coronal planes of brain tissue. This particular function averages across the anterior, middle, and posterior planes.

All parameters must be in quotations, besides "n_per_group." The .CSV files to be analyzed must be placed in your working directory.

When imaging your tissue initially, it is critical that the output filenames are saved in the following format: GroupName_ExperimentName_SubjectNumber_MagnificationX_CoronalPlane_SideOfTissue_RegionOfTissue

To mitigate experimenter bias (throughout all stages of the experimental process, but including this one) it is recommended that the group names are codified in a way agnostic to experimental treatment.
}
\examples{
immuno.analyze::overall.percent.area("percent_ps6_BLA", "percent", "Ctrl", "Switch", 6, data = "percent_ps6_BLA.csv")

}
\references{
Timothy, M., & Forlano, P. M. (2019). A versatile macro-based neurohistological image analysis suite for ImageJ focused on automated and standardized user interaction and reproducible data output. Journal of neuroscience methods, 324, 108286. https://doi.org/10.1016/j.jneumeth.2019.04.009
}

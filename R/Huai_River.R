#' Huai River Dataset
#'
#' This dataset contains aggregated PM10 values and related variables for 161 Disease Surveillance Points (DSPs)
#' in China, used to study the impact of sustained exposure to air pollution on life expectancy based on China's
#' Huai River Policy.
#'
#' @format A data frame with 161 rows and 4 variables:
#' \describe{
#'   \item{county_code}{Numeric, the unique code identifying each county or DSP.}
#'   \item{north_huai}{Numeric, indicator variable: 1 if the county is north of the Huai River, 0 otherwise.}
#'   \item{pm10}{Numeric, the concentration of particulate matter smaller than 10 μm (PM10) in μg/m³.}
#'   \item{dist_huai}{Numeric, the distance in degrees from the Huai River.}
#' }
#'
#' @details
#' This dataset is part of the replication files for the study "New Evidence on the Impact of Sustained
#' Exposure to Air Pollution on Life Expectancy from China’s Huai River Policy." The study leverages
#' quasi-experimental variation in PM10 levels caused by China’s Huai River Policy, which provides free
#' or subsidized coal for winter heating north of the Huai River but not south of it. The analysis,
#' based on a regression discontinuity design, finds that a 10-μg/m³ increase in PM10 reduces life
#' expectancy by 0.64 years (95% CI = 0.21–1.07), primarily due to cardiorespiratory mortality.
#' The dataset includes PM10 measurements for 161 DSPs and is used in tables and figures of the original paper.
#'
#' @source
#' Ebenstein, Avraham, Michael Greenstone, Maoyong Fan, Guojun He, and Maigeng Zhou. 2017.
#' "New evidence on the impact of sustained exposure to air pollution on life expectancy from China’s
#' Huai River Policy." *Proceedings of the National Academy of Sciences* 114(39): 10384-10389.
#' DOI: 10.1073/pnas.1616784114.
#'
#' @examples
#' data(Huai_River)
#' head(Huai_River)
#'
#' @seealso
#' The original replication files are available in the "Supporting_Datasets" folder of the study’s archive (https://www.pnas.org/doi/10.1073/pnas.1616784114#data-availability),
#' with "DSP_PM10.dta" being the source file for this dataset.
#'
"Huai_River"

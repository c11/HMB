#' Generate Report for Estimation Method Comparison
#'
#' @param n
#' @param SampFr
#' @param SampFrlim
#' @param EnLargeFactor
#' @param TileSide
#' @param PixelSize
#' @param meanY_S
#' @param sigmaY_S
#' @param meanX_S
#' @param sigmaX_S
#' @param meanZ_S
#' @param sigmaZ_S
#' @param meanX_Sa
#' @param sigmaX_Sa
#' @param meanZ_Sa
#' @param sigmaZ_Sa
#' @param meanZ_U
#' @param CorrYX
#' @param CorrYZ
#' @param CorrXZ
#' @param AutoCorrF
#' @param AutoCorrG
#' @param AutoCorrG_star
#' @param AutoCorr
#' @param module
#' @param ylim
#' @param nameForestAtt
#' @param nameRSdata
#' @param nameStudyArea
#' @param output_file
#'
#' @return
#' @export
#'
#' @examples
MethComp <- function(n, SampFr, SampFrlim, EnLargeFactor, TileSide, PixelSize,
                     meanY_S,  sigmaY_S,
                     meanX_S,  sigmaX_S,
                     meanZ_S,  sigmaZ_S,
                     meanX_Sa, sigmaX_Sa,
                     meanZ_Sa, sigmaZ_Sa,
                     meanZ_U,
                     CorrYX, CorrYZ, CorrXZ,
                     AutoCorrF = 0.8, AutoCorrG = 0.95, AutoCorrG_star = 0.8, AutoCorr = TRUE,
                     module = "All", ylim = c(15, 15, 50, 20, 100),
                     nameForestAtt = "AGB",
                     nameRSdata = c("GEDI", "TanDEM-X"),
                     nameStudyArea = " ",
                     output_file = NULL) {

  listData = list(n = n,
                  SampFr = SampFr,
                  SampFrlim = SampFrlim,
                  EnLargeFactor = EnLargeFactor,
                  TileSide = TileSide,
                  PixelSize = PixelSize,
                  nameStudyArea = nameStudyArea,
                  nameForestAtt = nameForestAtt,
                  nameRSdata = nameRSdata,
                  meanY_S  = meanY_S,
                  sigmaY_S  = sigmaY_S,
                  meanX_S  = meanX_S,
                  sigmaX_S  = sigmaX_S,
                  meanZ_S  = meanZ_S,
                  sigmaZ_S  = sigmaZ_S,
                  meanX_Sa = meanX_Sa,
                  sigmaX_Sa = sigmaX_Sa,
                  meanZ_Sa = meanZ_Sa,
                  sigmaZ_Sa = sigmaZ_Sa,
                  meanX_U = meanX_Sa,
                  meanZ_U = meanZ_U,
                  CorrYX = CorrYX,
                  CorrYZ = CorrYZ,
                  CorrXZ = CorrXZ,
                  AutoCorrX = AutoCorrF,
                  AutoCorrZ = AutoCorrG,
                  AutoCorrZ_Sa = AutoCorrG_star,
                  AutoCorr = AutoCorr,
                  #module = module,
                  ylim = ylim);

  library("rmarkdown")
  RMDname = system.file("rmd", "method_report.Rmd", package=getPackageName())
  render(RMDname, params=listData, output_file = output_file, output_format="html_document")
}

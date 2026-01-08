library(reticulate)
source_python("process.py")

#' Process data for use in keju.
#'
#' @param df a data.frame object containing TODO
#' @param control_treatment name of control treatment or phenotype in data
#' @param G number of enhancers per dispersion parameter
#' @returns a list TODO
#' @export
process <- function(df, 
                    control_treatment,
                    G=50
) {
    return(py_process(df, control_treatment, G))
}

#' Further process data to use keju with motif-level shrinkage. 
#'
#' @param df a data.frame object containing TODO
#' @param control_treatment name of control treatment or phenotype in data
#' @param G number of enhancers per dispersion parameter
#' @param motif column name of motif TODO
#' @returns a list TODO
#' @export
use_motif_shrinkage <- function(df,
                                control_treatment,
                                G=50,
                                motif='motif'

) {
    keju = py_process(df, control_treatment, G)
    keju = py_use_motif_shrinkage(keju, motif)
    return(keju)
}

#' Process data to use keju to infer covariate-level effects on transcription rate with motif-level shrinkage.
#'
#' @param df a data.frame object containing TODO
#' @param control_treatment name of control treatment or phenotype in data
#' @param G number of enhancers per dispersion parameter
#' @param motif column name of motif TODO
#' @param covariate column name of covariate TODO
#' @returns a list TODO
#' @export
use_covariate_slope_intercept <- function(df,
                                          control_treatment,
                                          G=50,
                                          motif='motif',
                                          covariate='covariate'
) {
    keju = py_process(df, control_treatment, G)
    keju = py_use_motif_shrinkage(keju, motif)
    keju = py_use_covariate_slope_intercept(keju, motif, covariate)
    return(keju)
}

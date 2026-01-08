#' Fit keju model to data.
#'
#' @param df TODO
#' @returns a list TODO
#' @export
fit <- function(
    keju,
    output_dir,
    model= 'no_motif',
    seed = 1,
    chains = 4,
    parallel_chains = 4,
    refresh = 100
) {
    # check model
    if (model %in% c('no_motif', 'motif_shrinkage', 'covariate_motif_slope_intercept')) {
        stan_path <- system.file('stan', paste0('keju_', model, '.stan'), package='keju')
    } else {
        warning(paste0('Model ', model, ' not recognized. Must be one of no_motif, motif_shrinkage, or covariate_motif_slope_intercept.'))
        # assert error
    }

    # create output directory and folder with outputs
    if (!dir.exists(output_dir)) {
        warning(paste0("Output directory ", output_dir, " does not exist. Creating...\n"))
        dir.create(output_dir, recursive=T)
    }
    output_folder <- gsub("//", "/", file.path(output_dir, "keju"))
    if (!dir.exists(output_folder)) {
        dir.create(output_folder)
    }

    mod <- cmdstan_model(stan_path)
    fit <- mod$sample(data = keju$data,
                      seed = seed,
                      chains = chains,
                      parallel_chains = parallel_chains,
                      refresh = refresh
    )

    keju$diagnostics <- fit$diagnostic_summary()
    print(keju$diagnostics)

    fit_summary = fit$summary()
    rownames(fit_summary) = fit_summary$variable
    keju$fit = fit_summary
    write.csv(fit_summary, TODO)

    pretty_summarize(keju, model, output_folder)
}

pretty_summarize <- function(
                             keju,
                             model,
                             output_folder
) {
    fit <- keju$fit
    tres <- keju$tres
    effects <- keju$effects
    correction_effects <- keju$correction_effects

    alphas <- fit[startsWith(rownames(fit), 'alpha['),]
    betas <- fit[startsWith(rownames(fit), 'beta['),]
    covariates <- fit[startsWith(rownames(fit), 'beta_correction['),]

    betas$is_significant = (betas$q5 > 0) | (betas$q95 < 0)

    rownames(alphas) <- keju$tres
    rownames(betas) <- keju$effects
    rownames(covariates) <- keju$covariates

    keju$alphas <- alphas
    keju$betas <- betas
    keju$covariates <- covariates

    write.csv(alphas, gsub("//", "/", file.path(output_folder, 'alphas.csv')))
    write.csv(betas, gsub("//", "/", file.path(output_folder, 'betas.csv')))
    write.csv(covariates, gsub("//", "/", file.path(output_folder, 'covariate_corrections.csv')))

    if (model != 'no_motif') {
        alpha_motifs <- fit[startsWith(rownames(fit), 'motif_transcription_rate['),]
        beta_motifs <- fit[startsWith(rownames(fit), 'motif_effect['),]

        betas_motifs$is_significant = (betas_motifs$q5 > 0) | (betas_motifs$q95 < 0)

        rownames(alpha_motifs) <- keju$alpha_motifs
        rownames(beta_motifs) <- keju$beta_motifs

        keju$alpha_motifs <- alpha_motifs
        keju$beta_motifs <- beta_motifs

        write.csv(alpha_motifs, gsub("//", "/", file.path(output_folder, 'motifs_alphas.csv')))
        write.csv(beta_motifs, gsub("//", "/", file.path(output_folder, 'motif_betas.csv')))

        if (model == 'covariate_motif_slope_intercept') {
            slope <- fit[startsWith(rownames(fit), 'slope['),]
            intercept <- fit[startsWith(rownames(fit), 'intercept['),]

            rownames(slope) <- keju$covariates
            rownames(intercept) <- keju$covariates

            keju$covariate_slope <- slope
            keju$covariate_intercept <- intercept

            write.csv(slope, gsub("//", "/", file.path(output_folder, 'covariate_slopes.csv')))
            write.csv(intercept, gsub("//", "/", file.path(output_folder, 'covariate_intercepts.csv')))
        }
    }

}
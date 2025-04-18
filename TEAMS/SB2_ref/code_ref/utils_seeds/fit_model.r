fit_biospear <- function(data,
                         biomarkers,
                         surv,
                         cvrts,
                         inter,
                         treatment,
                         methods) {
    set.seed(1993)
    fit <- BMsel(data = data,
                    x = biomarkers,
                    y = surv,
                    z = cvrts,
                    tt = treatment,
                    inter = inter,
                    std.x = TRUE,
                    std.i = FALSE,
                    std.tt = FALSE,
                    method = methods,
                    folds = 10,
                    uni.fdr = 0.05,
                    uni.test = 1,
                    ss.rando = FALSE,
                    ss.nsub = 100,
                    ss.fsub = 0.5,
                    ss.fwer = 1,
                    ss.thr = 0.6,
                    dfmax = 70,
                    pct.rep = 1,
                    pct.qtl = 0.95,
                    showWarn = TRUE,
                    trace = TRUE)
    # Get the summary of the Cox model
    #cox_summary <- summary(fit)
    #return(list(model = fit, summary = cox_summary))
    return(fit)
}
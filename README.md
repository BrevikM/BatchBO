# BatchBO

![version](https://img.shields.io/badge/version-0.1.0-blue)
![R](https://img.shields.io/badge/R-%E2%89%A5%204.2-276DC3)
![license](https://img.shields.io/badge/license-MIT-green)

`BatchBO` is an R package for Bayesian optimization of laboratory experiments
that are run in batches. In real experiments, batches can often be prone to shifts up or down, affecting all experiments inside.  
`BatchBO` separates these batch offsets from the effect of the settings, so the optimization follows the
underlying target response. 

### BatchBO enables

**Batch-aware modelling**
- A Gaussian process with a random batch offset: response = surface + batch offset + noise.
- A reference condition, repeated in every batch, anchors the batches and measures run-to-run noise.

**Experiments planned in batches**
- Space-filling starting designs, and adaptive batches by noisy expected improvement, Kriging believer or Thompson sampling.
- Exportable run sheets.

**Interpretable results**
- A recommendation for final operating conditions.
- Sobol sensitivity indices, main effects, response surfaces, and the best result as the campaign progressed.

**Resumable campaigns**
- One object holds the data, every fitted model, proposal and analysis; save it with `saveRDS()` and continue later.
- Add each batch when its results come back from the lab.

# Installation

```r
# install.packages("remotes")
remotes::install_github("BrevikM/BatchBO", build_vignettes = TRUE)
```

# Getting started

The vignette walks through a complete simulated campaign, from the first
batches to the final recommendation and analyses:

```r
vignette("hartmann_ESI", package = "BatchBO")
```

Every function has a help page with examples:

| Step | Help page |
|---|---|
| Declare factors and ranges | `?design_space` |
| Create the campaign | `?campaign` |
| Propose batches | `?proposeBatch` |
| Add completed batches | `?addBatch` |
| Fit the surrogate | `?runSurrogate` |
| Recommend the best condition | `?recommend` |
| Analyse the surface | `?runSobol`, `?runMainEffects`, `?runResponseSurfaces` |
| Follow progress | `?runCumulativeBest` |

`?BatchBO` gives an overview, and `help(package = "BatchBO")` lists every page.

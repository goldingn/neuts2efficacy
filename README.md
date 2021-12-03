# Analyses to predict the efficacy and waning of vaccines and previous infection against transmission and clinical outcomes of SARS-CoV-2 variants.

This repo contains analyses to:
 - fit a Bayesian model of the relationship between neutralising antibody titres and protective efficacies of vaccines and prior infection
 - use this model to predict waning efficacy, and efficacy angainst new variants
 - use a next generation matrix approach to predict the reduction in transmission achievable by different levels of vaccine efficacy against acquisition and onward transmission.

This model is a Bayesian implementation of the model described by [Khoury et al. 2020 Nature Medicine](https://doi.org/10.1038/s41591-021-01377-8) and used in [Cromer et al. 2021 Lancet Microbe](https://doi.org/10.1016/S2666-5247(21)00267-6). The code for that work can be found [here](https://github.com/InfectionAnalytics/COVID19-ProtectiveThreshold) and [here](https://github.com/InfectionAnalytics/SARS-CoV-2-Variants-and-Boosting---Lancet-Microbe). If you want to refer to the analyses in this repo, you should probably be citing them.

To run the code in this repo, you will need to install the packages listed in `packages.R`, including the most recent Github version of [greta](https://github.com/greta-dev/greta). greta links to python and depends on TensorFlow so can be tricky to install. It's best to follow the instructions greta gives you, and look for help on the [greta forum](https://forum.greta-stats.org/) if you have any trouble.

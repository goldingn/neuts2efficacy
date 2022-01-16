# Analyses to predict the efficacy and waning of vaccines and previous infection against transmission and clinical outcomes of SARS-CoV-2 variants.

This repo contains analyses to:
 - fit a Bayesian model of the relationship between neutralising antibody titres and protective efficacies of vaccines and prior infection
 - use this model to predict VEs for other outcomes, vaccines, and waning
 - infer the plausible range of immune escape and intrinsic transmissibility for Omicron
 
I'm planning to extend it to:
 - use a next generation matrix approach to predict the reduction in transmission achievable by different levels of vaccine efficacy against acquisition and onward transmission.

See the [Methods Description](methods.md) for a technical explanation of the analyses.

This model is a Bayesian implementation of the model described by [Khoury et al. 2020 Nature Medicine](https://doi.org/10.1038/s41591-021-01377-8) and used in [Cromer et al. 2021 Lancet Microbe](https://doi.org/10.1016/S2666-5247(21)00267-6). The code for that work can be found [here](https://github.com/InfectionAnalytics/COVID19-ProtectiveThreshold) and [here](https://github.com/InfectionAnalytics/SARS-CoV-2-Variants-and-Boosting---Lancet-Microbe). If you want to refer to the analyses in this repo, you should probably be citing them.

The omicron analysis relies on detailed, rapid analyses of the ratio of reproduction numbers for Omicron and Delta ([here](https://twitter.com/cap1024/status/1466840869852651529)), and the increasing risk of reinfection ([here](https://www.medrxiv.org/content/10.1101/2021.11.11.21266068v2)) in South Africa, by South African scientists as part of [SACMC](https://sacovid19mc.github.io/). 

To run the code in this repo, you will need to install the packages listed in `packages.R`, including the most recent Github version of [greta](https://github.com/greta-dev/greta). greta links to python and depends on TensorFlow so can be tricky to install. It's best to follow the instructions greta gives you, and look for help on the [greta forum](https://forum.greta-stats.org/) if you have any trouble.

Then you run the whole analysis using the targets package function `tar_make()` to run the whole analysis; fittng the model, producing validation statistics, producing plots, and outputting predictions.

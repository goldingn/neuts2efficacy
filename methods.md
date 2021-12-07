methods
================
Nick Golding
07/12/2021

# Overview

This document describes the methods used to delineate plausible values
of the intrinsic transmissibility and immune evasion of the Omicron
variant of SARS-CoV-2. It is structured in three parts: - A summary of
Bayesian implementation of the model of vaccine efficacy - The extension
of this Bayesian model to available data on Omicron - Details on model
fitting by MCMC using the greta R package

# Vaccine efficacy model

The core of this analysis is a Bayesian implementation of a predictive
model of vaccine efficacy documented and validated in [Khoury et
al. (2020)](https://doi.org/10.1038/s41591-021-01377-8) and [Cromer et
al. (2021)](https://doi.org/10.1016/S2666-5247(21)00267-6). The original
publications estimate parameters by maximum likelihood, and focus on
predicting vaccine efficacies for different combinations of vaccine dose
and product, the outcome against which efficacy is measured, and the
SARS-CoV-2 variant, based on neutralising antibody titres. See those
publications for detailed analysis and validation.

This model assumes that each immune individual *i* in a population has
some neutralisation level*n*<sub>*i*, *v*</sub> to variant *v*, given by
the common (ie. base 10) logarithm of the titre of neutralising
antibody, relative to the mean neutralisation level produced by
infection with wild-type SARS-CoV-2. This indexing on wild-type
convalescent neutralisation levels enables comparability between studies
and across variants. The individual’s neutralisation level is assumed to
be drawn from a normal distribution with mean
*μ*<sub>*s*, *d*, *v*</sub> differing based on the source of their
immunity *s* (e.g. the vaccine dose or product), the number of days *d*
post-peak immunity (ie. the degree of waning) adnd the variant, and with
variance *σ*<sup>2</sup> giving the inter-individual variation in
neutralisation levels - assumed to be constant across variants, sources
of immunity, and levels of waning:

*n*<sub>*i*, *v*</sub> ∼ *N*(*μ*<sub>*s*, *d*, *v*</sub>, *σ*<sup>2</sup>)

For each individual and for each type of outcome *o* (e.g. death, severe
disease, infection, onward transmission) the probability that the
outcome is averted *E*<sub>*i*, *o*</sub> (the vaccine is effective for
that outcome) is given by a sigmoid function, parameterised by a
threshold neutralisation level *n*<sub>*o*, 50</sub> at which 50% of
outcome events are prevented, and a slope parameter *k* determining the
steepness of this relationship. Crucially, these parameters are assumed
to independent of the variant and source of immunity, enabling
prediction of vaccine efficacy to new situations.

*E*<sub>*o*</sub>(*n*<sub>*i*, *v*</sub>) = 1/(1 + *e*<sup> − *k*(*n*<sub>*i*, *v*</sub> − *n*<sub>*o*, 50</sub>)</sup>)

Note that this model can equivalently be interpreted as there being a
different deterministic threshold for each different *event* whereby
that outcome would occur in the absence of vaccination, and would not
occur if the neutralisation level is above the random threshold, where
the random thresholds are drawn from a logistic distribution with
location parameter *n*<sub>*o*, 50</sub> and scale parameter 1/*k*.

At the population level, the efficacy of immunity against a given
outcome from a given variant in a cohort with mean neutralisation level
*μ*<sub>*s*, *d*, *v*</sub> is the average probability over the whole
population of the outcome being averted. This is computed by integrating
the sigmoid function with respect to the normal distribution of
neutralisation levels. This integral has no closed form, so must be
computed by numerical approximation.

*P*<sub>*s*, *d*, *o*, *v*</sub> = ∫<sub> − ∞</sub><sup>∞</sup>*E*<sub>*o*</sub>(*n*<sub>*v*</sub>)*f*(*n*<sub>*v*</sub>\|*μ*<sub>*s*, *d*, *v*</sub>, *σ*<sup>2</sup>)*d**n*<sub>*v*</sub>

The mean neutralisation level *μ*<sub>*s*, *d*, *v*</sub> for a cohort
with immunity source *s* and number of days *d* since peak immunity from
that source is assumed to follow exponential decay with half-life of *h*
days, from a peak mean neutralisation level against that variant of
*μ*<sub>*s*, *v*</sub><sup>\*</sup> for each source:

*μ*<sub>*s*, *d*</sub> = *l**o**g*10(10<sup>*μ*<sub>*s*, *v*</sub><sup>\*</sup></sup>*e*<sup> − *d*/*h*</sup>)

The peak mean neutralisation level against a given variant is in turn
modelled as a log10 fold increase or decrease in neutralising antibody
titres that variant, relative to an index variant:

*μ*<sub>*s*, *v*</sub><sup>\*</sup> = *μ*<sub>*s*, 0</sub><sup>\*</sup> + *F*<sub>*v*</sub>

Where for the index variant *F*<sub>*v*</sub> = 0, for a variant with
immune escape relative to the index the *F*<sub>*v*</sub> &lt; 1, and
for a variant more susceptible to neutralisation than the index,
*F*<sub>*v*</sub> &gt; 1.

Among these model parameters, *μ*<sub>*s*, 0</sub><sup>\*</sup>,
*F*<sub>*v*</sub>, *h*, *σ*<sup>2</sup> can be estimated from
neutralisation assay experiments. The parameters *n*<sub>*o*, 50</sub>
(one per outcome) and *k* (one in total) must be learned by fitting the
model to data on population-level vaccine efficacy. In this Bayesian
implementation, the parameter *σ*<sup>2</sup> is kept fixed at the value
estimated by [Cromer et
al. (2021)](https://doi.org/10.1016/S2666-5247(21)00267-6),
*μ*<sub>*s*, 0</sub><sup>\*</sup>, *F*<sub>*v*</sub>, and *h* are given
informative priors based on estimates from [Khoury et
al. (2020)](https://doi.org/10.1038/s41591-021-01377-8) and [Cromer et
al. (2021)](https://doi.org/10.1016/S2666-5247(21)00267-6), and the
remaining parameters are given less informative priors. This enables the
model to update these parameters slightly based on the data, and fully
incorporate uncertainty in these parameters into predictions.

This model is fitted to estimates from [Andrews et
al. (2021)](https://doi.org/10.1101/2021.09.15.21263583) of the
population-level efficacy of the Pfizer and AstraZeneca vaccines (two
doses) against clinical outcomes (death, severe disease, symptomatic
infection) from the Delta variant over different periods of time
post-administration, and estimates from [Pouwels et
al. (2021)](ttps://doi.org/10.1101/2021.09.28.21264260) and from [Eyre
et al. (2021)](https://doi.org/10.1101/2021.09.28.21264260) against
acquisition (symptomatic or asymptomatic) and onward transmission of
breakthrough infections, respectively, of the Delta variant.

The model likelihood for vaccine efficacy estimate *j* is defined as a
normal distribution over the logit-transformed estimate
VE<sub>*j*</sub>, with mean given by the logit-transformed predicted
efficacy for that combination of source, days post-administration,
outcome, and variant
*P*<sub>*s*<sub>*j*</sub>, *d*<sub>*j*</sub>, *o*<sub>*j*</sub>, *v*<sub>*j*</sub></sub>
and with variance given by the sum of the square of the standard error
of the estimate on the logit scale logit-SE<sub>*j*</sub><sup>2</sup>
(approximated from provided uncertainty intervals), and an additional
variance term *σ*<sub>ve</sub><sup>2</sup>, to represent any additional
errors in these estimates that may arise from inference on observational
data:

logit(VE<sub>*j*</sub>) ∼ *N*(logit(*P*<sub>*s*<sub>*j*</sub>, *d*<sub>*j*</sub>, *o*<sub>*j*</sub>, *v*<sub>*j*</sub></sub>), logit-SE<sub>*j*</sub><sup>2</sup> + *σ*<sub>ve</sub><sup>2</sup>)

All VE estimates covered events over a period of time, and
*d*<sub>*j*</sub> was taken as the midpoint of that period.

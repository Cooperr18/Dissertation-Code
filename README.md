# Beyond the Signals: An evaluation of the Signal Selection Test using simulated data

## R Markdown

This repo includes all scripts, figures and tables created to carry out
the analysis of my dissertation, within the framework of the MPhil in
Human Evolutionary Studies at the University of Cambridge. The aim of
this project is to critically assess the performance of the core
statistical component of the
[`signatselect`](https://github.com/benmarwick/signatselect), the
**Frequency Increase Test (FIT)**.

The FIT is a statistical test Created by [Feder et
al. (2014)](https://doi.org/10.1534/genetics.113.158220), which aims to
detect signals of selection in time series data. It employs a *t*-test
to compare two successive frequency changes (i.e. three time steps)
observed in our data against null expectations, in this case, neutral
evolution driven by stochastic processes alone. If the changes
significantly deviate from the null, they are flagged by the FIT and
interpreted as **selection** events.

This investigation is driven by the lack of systematic research of the
`signatselect` in archaeology. There are few sources which explore the
scope of this technique: the original publication, [Feder et
al. (2014)](https://doi.org/10.1534/genetics.113.158220); an application
to language time series data [Newberry et
al. (2017)](https://doi.org/10.1038/nature24455); a brief overview of
softwares to detect selection [Vlachos et
al. (2020)](https://doi.org/10.1186/s13059-019-1770-8), and an
evaluation of Newberry et al.’s results plus the effect of time binning
on the test [Karjus et al. (2019)](https://doi.org/10.5334/gjgl.909).
However, the only source which explicitly applies and evaluates the
`signatselect` with archaeological data is a public repo posted by Ben
Marwick, Hezekiah A. Bacovcin and Sergey Kryazhimskiy, at
<https://github.com/benmarwick/signatselect>. This source has served of
great inspiration to this study.

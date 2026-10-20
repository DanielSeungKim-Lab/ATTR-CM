Black representation across benchmarks in ATTR-CM trials
Data and code to reproduce the reported results.

RUN
  python CODE-01-reproduce-results.py     prints every number reported in the paper
  python CODE-02-make-figure.py           writes FIGURE1.png
Requires Python 3.13 with pandas, numpy, matplotlib. Both scripts run from this directory
and take no arguments. CODE-01 writes nothing; the figure is the only generated file.

FILES
  CODE-01-reproduce-results.py   pooled proportion and bootstrap interval, the 8 benchmarks
                                 and their participation ratios, per-trial catchment
                                 comparisons, and the ENDEAVOUR comparisons
  CODE-02-make-figure.py         the two-panel figure, recomputed from the same inputs
  DATA-01-trial-race-counts.csv  race counts for the 6 landmark ATTR-CM trials reporting race
                                 in a cardiomyopathy population
  DATA-02-benchmark-inputs.csv   the published counts and percentages the benchmarks are
                                 built from, one row per quantity
  DATA-03-site-counties.csv      the US counties hosting each trial's sites, with their
                                 Black-alone and White-alone populations aged 75 and over
  DATA-04-landmark-disposition.csv  all 23 trials meeting the prespecified landmark
                                 definition, with tier and how each was used
  DOC-02-data-dictionary.txt     column definitions

SCOPE
  Five trials are pooled; ENDEAVOUR is reported separately because enrollment was restricted
  to pathogenic TTR variant carriers. Proportions throughout are Black divided by Black plus
  White. Participation ratios divide the pooled proportion by a reference proportion, so a
  ratio below 1 means a smaller Black share than the reference population.

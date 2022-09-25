---
title: 'edecob: Detection of sustained change in high-frequency and high-variance time series'
tags:
  - R
  - Digital Biomarker
authors:
  - name: Zheng Chen Man
    orcid: 0000-0001-5025-3960
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
  - name: Fabian Model
    corresponding: true # (This is how you can denote equal contributions between multiple authors)
    affiliation: 3
  - name: Frank Dondelinger
    corresponding: true # (This is how to denote the corresponding author)
    affiliation: 3
  - name: Stanislas Hubeaux
    corresponding: true # (This is how to denote the corresponding author)
    affiliation: 2
affiliations:
 - name: Name, Country
   index: 1
 - name: F. Hoffmann-La Roche AG, Switzerland
   index: 2
 - name: Independent Researcher, Country
   index: 3
date: 15 July 2022
bibliography: paper.bib

---

# Summary

`edecob` (**E**vent **De**tection using **Co**fidence **B**ounds) is an R package that detects sustained change (denoted “event” in this paper and in the package documentation) in high-frequency, longitudinally-collected data. We approximate the data using a smoother (e.g. moving median). The residuals of this smoother are then modelled using an autoregressive model and the corresponding error is bootstrapped to obtain pointwise confidence intervals. These confidence intervals are then used to construct simultaneous confidence bands for the smoothed data [@buhlmann1998]. We define an event as the occasion where confidence bands stay within a prespecified range compared to a reference (e.g. the baseline) for a predefined amount of time. Our approach takes into account the increased variability and variations in data density that can occur when digital assessments are performed outside of a controlled setting [@Roussos_2022]. Formally, the methodology assumes that the measurements are equally spaced and the time series remains stationary. In practice, we observed that the method provides meaningful results under a variety of conditions. The parameters of the model are adjustable to enable customization of the methodology to specific use cases. 


# Statement of need

Our methodology is suitable for event detection in continuous time-dependent high-frequency datasets, leading to a wide range of applications, including finance, meteorology, and digital biomarkers. We will focus here on the latter, as this was the setting that prompted the initial development of the package.

A common question occurring in the diagnosis and monitoring of chronic diseases is the detection of sustained change on a disease severity measurement that needs to be distinguished from short-term improvement or worsening of symptoms. Examples of sustained worsening are the monitoring of disability progression in Multiple Sclerosis [@MSBase] or the diagnosis of chronic obstructive pulmonary disease (COPD) [@COPD]. Outcomes of sustained worsening are of particular interest for clinical trials where a new intervention needs to show that it can reduce the risk of long-term disability accumulation. On the other hand, therapeutic interventions may lead to sustained improvement of symptoms [@improvement].

Digital sensor-based technologies such as wearables allow for frequent functional assessments of impairment outside of the clinic, i.e. by patients themselves at home. Digital biomarkers collected with these technologies are an active field of research, in particular for neurological conditions such as multiple sclerosis [@biomarkersinMS]. Examples of digital biomarker analysis include distinguishing patients from healthy subjects [@currentsoftware] and using machine-learning based approach to quantify impairment [@machinelearningapproach]. For digital biomarkers collected at weekly or even daily frequency, traditional definitions of sustained change designed for in-clinic supervised assessments performed every 3 or 6 months cannot be applied. To close this gap, our package provides a method for reliably detecting sustained change events in time series of digital biomarkers. 



# Overview

The package can be downloaded from [CRAN](https://CRAN.R-project.org/package=edecob). The user can provide the data to the main function edecob in the long format. The first column specifies the identifier (e.g. a subject identifier in case of a digital biomarker), the second column the time point at which the measure is collected, and the third column the values for this time point. Additionally, the user can supply the detection range in the fourth and fifth columns. Note that this enables detection of either increase or decrease by choosing appropriate bounds.

```
head(example_data, 3)
#>     subject study_day jump_height detect_lower detect_upper
#> 1 Subject 1         1    55.60844         -Inf     54.41227
#> 2 Subject 1         4    57.77688         -Inf     54.41227
#> 3 Subject 1         7    57.59584         -Inf     54.41227
```

The data frame constitutes the first argument of the main edecob function. The user can then specify the following:

* `min_change_dur` specifies the minimal number of time units the confidence bounds need to stay within the detection range to detect an event.

* `detect` specifies how the detection range is determined. When using `above` or `below`, one of the detection bounds will be the median over the data points in the first `bline_period` time units, then multiplied with `detect_factor`. The other detection bound will be `Inf` when choosing `above` or `-Inf` when choosing `below`. `detect` can also be chosen to be `custom`, in which case the detection range will be provided by the user for every subject.

* `resample_method` specifies whether bootstrapping for simultaneous confidence should be done over all residuals `all` or only over a window `window`. If a window is chosen, the window from which residuals should be bootstrapped from can be specified with `resample_win`.

The use of other variables is described in the documentation. If the user wants the confidence band to reflect the local frequency of the data, resampling from all data may be more appropriate. However, if the user wants the confidence band to reflect the local variation of the data rather than the global variation, restricting the resampling over a window may be the right choice. This is especially relevant when analyzing time-dependent data.
```
# We apply the main function of the package onto our example_data
example_event <- edecob(example_data, 
			 smoother = "mov_med",
			 resample_method = "window",
			 min_change_dur = 50,
			 conf_band_lvl = 0.95,
			 bt_tot_rep = 50,
			 time_unit = "day",
			 detect = "custom",
			 detect_factor = 1,
             bline_period = 14,
			 resample_win = c(-5,5))
```
The `example_event` object contains 3 objects, one corresponding to each subject and a data frame summarizing the event information.
```
names(example_event)
#> [1] "Subject 1"  "Subject 2"  "event_info"
```
We can then choose a subject and plot the corresponding confidence bands.

```
plot(example_event$`Subject 1`)
```
![A plot generated using simulated data. The confidence band of the smoothed data is used to detect events. If the confidence bands stay within the detection bounds for a prespecified amount of time, an event will be detected. \label{fig:plot}](plot.png){width=80%}

The event information is summarized and accessible in a table called `event_info`. It contains the information about whether an event was detected, the time point of the onset of the event, how long the change is sustained, and whether the event is sustained until the end of available observations for the subject.

```
example_event$event_info
#>           event_detected event_onset event_duration event_stop
#> Subject 1           TRUE         169             87       TRUE
#> Subject 2           TRUE         205             51       TRUE
#> Subject 3          FALSE         306             38      FALSE
```

Using this table, we can generate a survival plot using the survival package [@survival-package].

```
library("survival")
plot(survfit(Surv(time = event_onset, event = event_detected) ~ 1,
             data = example_event$event_info),
     conf.int = FALSE, xlim = c(0,350), ylim = c(0,1), mark.time = TRUE,
     xlab = "Study Day", ylab = "Survival Probability", main = "Survival plot")
```

![The survival plot generated using the code above.\label{fig:survplot}](survplot.png){width=80%}

In case of a large number of subjects, the main function edecob is parallelizable by separately calling the function for every subject. 

# References

# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements


# References
---
title: 'edecob: Detection of sustained change in high-frequency and high-variance time series'
tags:
  - R
  - Digital Biomarker
authors:
  - name: Zheng Chen Man
    orcid: 0000-0001-5025-3960
    affiliation: "1" # (Multiple affiliations must be quoted)
  - name: Fabian Model
    affiliation: 2
  - name: Frank Dondelinger
    orcid: 0000-0003-1816-6300
    affiliation: 2
  - name: Stanislas Hubeaux
    affiliation: 2
  - name: Peter Bühlmann
    orcid: 0000-0002-1782-6015
    affiliation: "1"
affiliations:
 - name: ETH Zurich, Switzerland
   index: 1
 - name: F. Hoffmann-La Roche AG, Switzerland
   index: 2

date: 26 October 2022
bibliography: paper.bib

---

# Summary

`edecob` (**E**vent **De**tection using **Co**nfidence **B**ounds) is an R package that detects sustained change (denoted as “event” in this paper and in the package documentation) in high-frequency, longitudinally-collected data. We model the data as a smooth trend plus dependent noise. The first part is approximated with a smoother (e.g. running median). The residuals of this smoother are then modelled using an autoregressive model and the corresponding error is bootstrapped to obtain pointwise confidence intervals. These confidence intervals are then used to construct simultaneous confidence bands for the underlying smooth trend [@buhlmann1998]. We define an event as the occasion where confidence bands stay within a prespecified range compared to a reference (e.g. the baseline) for a predefined amount of time. Our approach takes into account the varying density of data ocurrence which may happen when convenience measurements are performed outside of a controlled setting, e.g. patients performing self-assessments using digital health technology at home instead of in a controlled clinical environment [@Roussos_2022]. In practice, we observed that the method provides meaningful results under a variety of settings. The parameters of the model are adjustable to enable customization of the methodology to specific use cases. 


# Statement of need

Our methodology is suitable for event detection in continuous time-dependent high-frequency datasets, leading to a wide range of applications, including finance, meteorology, and digital biomarkers. We will focus here on the latter, as this was the setting that prompted the initial development of the package.

A common question occurring in the diagnosis and monitoring of chronic diseases is the detection of sustained change on a disease severity measurement that needs to be distinguished from short-term improvement or worsening of symptoms. Examples of sustained worsening are the monitoring of disability progression in Multiple Sclerosis [@MSBase] or the diagnosis of chronic obstructive pulmonary disease (COPD) [@COPD]. Outcomes of sustained worsening are of particular interest for clinical trials where a new intervention needs to show that it can reduce the risk of long-term disability accumulation. On the other hand, therapeutic interventions may lead to sustained improvement of symptoms [@improvement].

Digital sensor-based technologies such as wearables allow for frequent functional assessments of impairment outside of the clinic, i.e. by patients themselves at home. Digital biomarkers collected with these technologies are an active field of research, in particular for neurological conditions such as multiple sclerosis [@biomarkersinMS]. Examples of digital biomarker analysis include distinguishing patients from healthy subjects [@currentsoftware] or using machine-learning based approaches to quantify impairment [@machinelearningapproach]. For digital biomarkers collected at weekly or even daily frequency, traditional definitions of sustained change designed for in-clinic supervised assessments performed every 3 or 6 months cannot be applied. To close this gap, our package provides a method for reliably detecting sustained change, i.e. events, in time series of digital biomarkers. 



# Overview

The R package can be downloaded from [CRAN](https://CRAN.R-project.org/package=edecob). The user can provide the data to the main function edecob. The first column specifies the identifier (e.g. a subject identifier in case of a digital biomarker), the second column the time point at which the measure is collected, and the third column the values for this time point. Additionally, the user can supply the detection range in the fourth and fifth columns. Note that this enables detection of either increase or decrease by choosing appropriate bounds.

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

The use of other variables is described in the documentation. If the user wants the confidence band to reflect the local frequency of the data, resampling from all data may be more appropriate. However, if the user wants the confidence band to reflect the local variation of the data rather than the global variation, restricting the resampling over a window may be the right choice. This is especially relevant when analyzing data exhibiting varying density of occurences.
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
We can then choose a subject and plot the corresponding confidence bands \autoref{fig:example}.

```
plot(example_event$`Subject 1`)
```
![A plot generated using simulated data. The confidence band of the smooth trend (in the data) is used to detect events. If the confidence bands stay within the detection bounds for a prespecified amount of time, an event will be detected. \label{fig:plot}](plot.png){width=80%}

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

![The results obtained using the package can be easily further analyzed using the survival package. Here, we generate a survival plot using the example data provided with the package. \label{fig:survplot}](survplot.png){width=80%}

In case of a large number of subjects, the main function edecob can trivially be parallelized by separately calling the function for every subject. In fact, the method of edecob uses one time series from a single subject to detect events for this particular subject.

# References

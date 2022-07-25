---
title: 'edecob: Event Detection using Confidence Bounds'
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

`edecob` (Event Detection using Confidence Bounds) is an R package that detects sustained change (which we call “event”) in high-frequency, longitudinally-collected-digital biomarker data. We approximate a subject’s performance using a smoother and detect sustained change by constructing confidence bands for the smoother as described in @buhlmann:1998. We define an event as the occasion where confidence bands stay within a prespecified range compared to the baseline for a predefined amount of time. Our approach is robust to noise and the increased variability that can result when the digital assessments are performed outside of a controlled setting (citation needed). The parameters of the model are adjustable to enable customization of the methodology to specific user’s cases. 

# Statement of need

A common question occurring in the diagnosis and monitoring of chronic diseases is the detection of sustained worsening on a disease severity measurement that needs to be distinguished from more acute short-term worsening of symptoms. Examples are the monitoring of disability progression in Multiple Sclerosis [@MSBase:2006] or the diagnosis of chronic obstructive pulmonary disease (COPD) [@COPD:2006]. Outcomes of sustained worsening are of particular interest for clinical trials where a new intervention needs to show that it can reduce the risk of long-term disability accumulation. Interventions that can improve symptoms or have the potential to cure the reverse outcome of sustained improvement can be of interest [@improvement:2021].
Digital technologies allow for frequent functional assessments of impairment by patients themselves at home. For measurements at weekly or even daily frequency, traditional definitions of sustained change designed for professionally monitored assessments performed every 3 or 6 months cannot be applied. Digital biomarkers are an active field of research [@biomarkersinMS:2021]. Current software development has mainly focused on using digital biomarker data to distinguishing patients from healthy subjects [@currentsoftware:2021] or on using machine-learning based approach to quantify impairment [@machinelearningapproach:2018].


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
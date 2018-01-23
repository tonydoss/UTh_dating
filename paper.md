---
title: "UTh_dating"
tags:
- uranium-thorium dating
author: 
- name: Anthony Dosseto
  orcid: 0000-0002-3575-0106
  affiliation: "1""
affiliations: 
- name: Wollongong Isotope Geochronology Laboratory, School of Earth and Environmental Sciences. University of Wollongong. Wollongong, NSW, Australia
  index: 1
date: 05 December 2017
bibliography: SURFACE.bib
---

# Summary

Uranium-thorium dating is a major geochronological technique widely used in Earth sciences. It can be used to determine the age of formation of corals (and thus of important for palaeo-climate and sea-level studies) [@RN4495], speleothems (also of significance for palaeo-climate studies) [@RN3182]. If the sample verifies a set of assumptions (see [@RN4495, @RN3182], for instance), the measurement of (<sup>230</sup>Th/<sup>238</sup>U) and (<sup>234</sup>U/<sup>238</sup>U) activity ratios can be used to calculate a closed-system <sup>230</sup>Th-U age and the initial (<sup>234</sup>U/<sup>238</sup>U) ratio, i.e. the (<sup>234</sup>U/<sup>238</sup>U) ratio at the time of formation (parentheses denote activities throughout this article).


The R script provided allows you to calculate closed-system <sup>230</sup>Th-U ages. 

Simply, set the working directory “path”, name of the sample to solve “sample_name” (can solve several samples if they share the same character string as defined by the user), and run the code.
Data need to be in a tab-separated .txt file “IoliteExport_All_Integrations.txt” (although the file name can be changed), with sample names, measured (<sup>230</sup>Th/<sup>238</sup>U), (<sup>234</sup>U/<sup>238</sup>U) activity ratios  and their 2$\sigma$ errors in columns named respectively “X”, “Th230_U238_CORR”, “U234_U238_CORR”, “Th230_U238_CORR_Int2SE” and “U234_U238_CORR_Int2SE”.
You can also set the number of optimisation for each sample, “nbit”, and the lower and upper bounds, respectively “lowerbound” and “upperbound”, for log10 of the age (in log10(yr)) and initial (<sup>234</sup>U/<sup>238</sup>U).

The code returns ages (in kyr), initial (<sup>234</sup>U/<sup>238</sup>U) and their 2$\sigma$ errors in a comma-separated .csv file with “sample_name” as file name.

# References

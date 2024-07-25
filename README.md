# EfficiencySimulation
This repository contains R codes to reproduce the numerical experiments in the paper "Pitman Efficiency Lower Bounds for Multivariate Distribution-Free Tests Based on Optimal Transport" ([https://arxiv.org/pdf/2104.01986](https://arxiv.org/pdf/2104.01986)).
The experiments were replicated with the terremoto cluster to speed up computation. We have added an R file named Testingtwosamexec.R and a generic scipt file named Testingtwosam.sh which contains details of the implementation. The other files have self-explanatory names and inline comments.
As an example, a filename "GaussianLocation_(A1).R gives the code for setting (A1) from the aforementioned paper, which is based on Gaussian location alternatives.
The analysis of the SONAR dataset is contained in the file SONARanalysis.R

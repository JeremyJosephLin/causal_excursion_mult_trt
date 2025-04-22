# paper_MRTSampleSizeBinary

Reproducible code for paper "Causal excursion effect for continuous outcome with multilevel treatment" by Jeremy Lin and Tianchen Qian

## File structure

- application code: code to reproduce results in Section 5 "Data Example".
- functions: code containing relevant user-defined R functions
- simulation code: code to reproduce simulation results in Section 4 "Sample Size Formula with Categorical Treatment" (Sec. 4.3) and Appendix C and D.

## How to reproduce the results

- To reproduce results in Section 5 "Application", run all R scripts in folder "application code". The R scripts do not depend on one another so no particular ordering in running.
- To reproduce simulation results in Section 6 "Detailed Simulation Results" (Sec. 6.2-6.8) and Supplementary Material D, do the following for each subfolder inside "simulation code":
    - First, run the R script(s) inside all subfolder(s) "simulationX(.X)". This conducts Monte Carlo simulations and saves result file. **Caution: each R script may take a long time (days) to finish. Also, X.X in simulationX.X does not correspond to Sections in the paper; these indices were created during the development of the paper.**
    - Second, run the R script(s) named "make figure X.R". This makes plots using the simulation result files. **Caution: X in the "make figure X.R" does not correspond to figure index in the paper. See table below for the figure index correspondence.**
    - For example, to reproduce everything in Section 6.3, go inside folder "2. WA-a violated (Sec 6.3)", then run the R scripts in subfolders "simulation2.2", "simulation2.3", "simulation2.4", "simulation2.5", "simulation4.1". Then, run the three R scripts "make figure 1.R", "make figure 2.R", "make figure 3.R".
- To reproduce Figure 4 in in Sec. 6.1 of the paper, run the three R scripts in folder "misc code".


| Figure in paper | R script to make the figure                                                            |
|-----------------|----------------------------------------------------------------------------------------|
| 1               | application code/size Drink Less 1 - determine ATE, ASPN.R                             |
| 2               | application code/size Drink Less 2 - assess sensitivity to theta.R                     |
| 3               | application code/size Drink Less 3 - assess impact of AA.R                             |
| 4               | misc code/make figure tau.R, make figure theta-f.R, make figure theta-g.R              |
| 5               | simulation code/1. all working assumptions hold (Sec 6.2)/make figure 0.R              |
| 6               | simulation code/2. WA-a violated (Sec 6.3)/make figure 1.R                             |
| 7a              | simulation code/2. WA-a violated (Sec 6.3)/make figure 2.R                             |
| 7b              | simulation code/2. WA-a violated (Sec 6.3)/make figure 3.R                             |
| 8               | simulation code/3. WA-b violated (Sec 6.4)/make figure 4.R                             |
| 9a              | simulation code/3. WA-b violated (Sec 6.4)/make figure 5.R                             |
| 9b              | simulation code/3. WA-b violated (Sec 6.4)/make figure 6.R                             |
| 10              | simulation code/3. WA-b violated (Sec 6.4)/make figure 7.R                             |
| 11              | simulation code/4. WA-c violated (Sec 6.5)/make figure 8.R                             |
| 12              | simulation code/4. WA-c violated (Sec 6.5)/make figure 9.R                             |
| 13              | simulation code/5. WA-d violated (Sec 6.6)/make figure 10.R                            |
| 14              | simulation code/6. WA-e violated (Sec 6.7)/make figure 11.R                            |
| 15              | simulation code/7. multiple WAs violated (Sec 6.8)/make figure 12.R                    |
| D.1             | simulation code/8. WA-a violated (additional results in Appendix D.1)/make figure A1.R |
| D.2             | simulation code/9. WA-b violated (additional results in Appendix D.2)/make figure A2.R |

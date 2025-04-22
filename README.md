# paper_MRTSampleSizeBinary

Reproducible code for paper "Causal excursion effect for continuous outcome with multilevel treatment" by Jeremy Lin and Tianchen Qian

## File structure

- application code: code to reproduce results in Section 5 "Data Example".
- functions: code containing relevant user-defined R functions that includes sample size calculator and CEE estimator
- simulation code: code to reproduce simulation results in Section 4 "Sample Size Formula with Categorical Treatment" (Sec. 4.3) and Appendix C and D.

## How to reproduce the results

To reproduce the results in Section 5, "Data Example," run all R scripts in the folder "application code." The R scripts do not depend on one another, so there is no particular order for running.
- To reproduce simulation results in Section  "Detailed Simulation Results" (Appendix D), do the following for each subfolder inside "simulation code":
    - First, run the R script(s) inside all subfolder(s) "simulationX(.X)". This conducts Monte Carlo simulations and saves result file. **Caution: each R script may take a long time (days) to finish. Also, X.X in simulationX.X does not correspond to Sections in the paper; these indices were created during the development of the paper.**
    - Second, run the Rmd file(s) named "sim_X.X_report.RMD". This makes plots and results using the simulation result files. **Caution: X in the "make figure X.R" does not correspond to figure index in the paper. See table below for the figure index correspondence.**
    - For example, to reproduce everything in Section 6.3, go inside folder "2. WA-a violated (Sec 6.3)", then run the R scripts in subfolders "simulation2.2", "simulation2.3", "simulation2.4", "simulation2.5", "simulation4.1". Then, run the three R scripts "sim_X.X_report.RMD", "make figure 2.R", "make figure 3.R".
- To reproduce Figure 4 in in Sec. 6.1 of the paper, run the three R scripts in folder "misc code".


| Figure in paper | R script to make the figure                                                            |
|-----------------|----------------------------------------------------------------------------------------|
| 1.a             | application code/graph_constant.Rmd                                                    |
| 1.b             | application code/theta_sensitivity_graph.Rmd                                           |
| 1.c,d           | application code/AA_sensitivity_graph.Rmd                                              |
| S1, S2          | simulation/power simulation/1. all working assumptions hold/sim_1a_report.Rmd          |
| S3              | simulation/power simulation/2.WA-a violated/sim 3.1a/sim_3_1a_report.Rmd               |
| S4a             | simulation/power simulation/2.WA-a violated/sim 2.1a/sim_2_1a_report.Rmd               |
| S4b             | simulation/power simulation/2.WA-a violated/sim 2.2a/sim_2_2a_report.Rmd               |
| S4c             | simulation/power simulation/2.WA-a violated/sim 2.5a/sim_2_5a_report.Rmd               |
| S4d             | simulation/power simulation/2.WA-a violated/sim 2.3a/sim_2_3a_report.Rmd               |
| S4e             | simulation/power simulation/2.WA-a violated/sim 2.4a/sim_2_4a_report.Rmd               |
| S4f             | simulation/power simulation/2.WA-a violated/sim 2.6a/sim_2_6a_report.Rmd               |
| S5a             | simulation/power simulation/3.WA-b violated/sim 4.1a/result/sim_4_1_report.Rmd         |
| S5b             | simulation/power simulation/3.WA-b violated/sim 4.1a/result/sim_const_true_eo.Rmd      |
| S5c,d           | power simulation/3.WA-b violated/sim 4.1a/result/sim_const_working_eo.Rmd              |
| S6              | simulation code/5. WA-d violated (Sec 6.6)/make figure 10.R                            |
| 14              | simulation code/6. WA-e violated (Sec 6.7)/make figure 11.R                            |
| 15              | simulation code/7. multiple WAs violated (Sec 6.8)/make figure 12.R                    |
| D.1             | simulation code/8. WA-a violated (additional results in Appendix D.1)/make figure A1.R |
| D.2             | simulation code/9. WA-b violated (additional results in Appendix D.2)/make figure A2.R |



# simstudy_incubationtime

This folder contains the code for two simulation studies, described in the corresponding paper.

There are files related to three sections: simulation study I, simulation study II and application.

**
Simulation study I
**
- fun_simulate_data_new.R to generate interval censored observations of incubation time.
- fun_estimate_single_interval_censored.R to estimate incubation time using different approaches.
- simulate.R to perform the simulation study. We performed the study in the ALICE server (https://www.universiteitleiden.nl/en/research/research-facilities/alice-leiden-computer-cluster) but the code can also be run locally.
- Graphs_simstudy.Rmd to create the graphics for the paper.
- exposure_window_lengths_based_on_real_data.RData exposure window widths in the opensource data used in Application.

**
Simulation study II
**
- fun_simulate_data_renewal.R to generate forward time observations of incubation time.
- fun_estimate_renewal_new.R to estimate incubation time using an approach inspired by renewal process theory.
- simulate_renewal.R to perform the simulation study.
- Graphs_simstudy_renewal.Rmd to create the graphics for the paper.

**
Application
**
- opensource_data_combined_Dec2021.RData collected opensource data.
- Application_opensource_data.Rmd the estimation on opensource data and corresponding graphs. Note that here, the function from Simulation study I is used.

If you have questions or comments, feel free to reach out to me by email (v.h.arntzen@math.leidenuniv.nl).

Kind regards,
Vera Arntzen

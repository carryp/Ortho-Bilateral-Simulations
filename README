# Introduction

This file describes the R programs that support the paper entitled “When 1 + 1 Does Not Equal 2: Re-Evaluating the Impact of Including of Patients with Bilateral Conditions in Orthopedic Clinical Research Studies”. The purpose of this study was to evaluate the impact of including bilateral patients in orthopedic clinical studies under different scenarios defined based on the prevalence of bilateral patients in the study population as well as, among bilateral patients, the strength of correlation between both observations. We compared the performance of three common modelling strategies for managing bilateral observations: (1) analyze all limbs as independent patients (naïve), (2) randomly select 1 limb per patient (random), and (3) account for correlation between limbs using linear mixed models (LMM).

# Study Design

We simulated hypothetical randomized controlled trials in which Western Ontario and McMaster Universities Arthritis Index (WOMAC) scores were collected from patients in two groups (treatment vs control) at baseline and two-years post-operatively. The primary outcome represented change in WOMAC scores from baseline. We simulated two scenarios, there was truly no difference between groups (null hypothesis was true) and there was truly a difference between groups (alternative hypothesis was true). Under the alternative hypothesis, the difference between groups (treatment effect) was simulated to be a 10-point difference in the change in WOMAC scores. This represents >0.5 standard deviation difference between groups (standard deviation = 19). Assuming a common standard deviation of 19 and a mean difference of 10 units, the sample size for the two groups (n=60) in our hypothetical randomized controlled trials was selected as the sample size expected to achieve >80% power to reject the null hypothesis of no difference between groups based on a two-sample independent t-test.

There are two simulation scripts provided in this repository. The `Full Sim Script_High Correlation_Final` file provides code for simulating a high level of correlation, Pearson correlation coefficient ~0.82, between bilateral observations. The `Full Sim Script_Low Correlation_Final` file provides code for simulating a modest level of correlation, Pearson correlation of ~0.33, between bilateral observations.

Each condition was evaluated 10,000 times, this can be modified by changing the `replicates` input parameter within the code.

# Parameters within the Simulation Scripts

The simulation code provides the flexibility to update the parameters under new study conditions that were not evaluated in our study. These parameters are defined below:

- `N.Bilateral.List`: List of bilateral conditions (number of patients simulated to be bilateral patients in each simulated randomized controlled trial)
- `u1`: Random intercept representing subject-specific intercept
- `E`: Random error (altering E and u1 will change the variation and correlation between limbs)
- `B0`: Intercept for linear model (represents mean of the outcome in the reference/control group)
- `B1`: Group effect (represents mean difference between treatment vs control groups, simulated to be 10 under the alternative hypothesis and 0 under the null hypothesis)
- `B2`: Indicator variable for bilateral patients (1 = bilateral patient, 0 = unilateral, simulated to be 0 no difference between bilateral vs unilateral subjects)

# Output Files

The code will output two files (both n x 30 matrices, where n represents the number of bilateral conditions – see `N.Bilteral.List` input variable) that summarize the results of the simulations. For each model, the output file includes the mean difference between groups, standard error, p-value, and 95% confidence intervals. One output file represents model performance under the null hypothesis, the other output file represents performance under the alternative hypothesis.

# Load necessary libraries
install.packages('tidyverse')
install.packages('data.table')
install.packages('nlme') # For linear mixed models
library(tidyverse)
library(data.table)
library(nlme)

# Define the number of bilateral observations for each iteration
N.Bilateral.List = as.matrix(
  x = c(60, 57, 54, 51, 48, 45, 42, 39, 36, 33, 30, 27, 24, 21, 18, 15, 12, 9, 6)
)

# Define subject-specific intercept and error term
u1 = 16.75
E = 8.25

# Coefficients for the simulation
B0 = 5
B1 = 0
B2 = 0

# Number of replicates for the simulation (5 for testing, 10000 for actual simulation)
replicates = 10

# Create empty matrices to store simulation and error results
Simulation.Output = matrix(NA, nrow = replicates, ncol = 33)
Error.Output = matrix(NA, nrow = replicates, ncol = 3)

# Loop through each N.Bilateral value
for(j in 1:nrow(N.Bilateral.List)){
  
  # Number of bilateral subjects for this iteration
  N.Bilateral = N.Bilateral.List[j]
  
  # Inner loop for each replicate
  for (n in 1:replicates){
    
    # Reset seed to ensure reproducibility
    set.seed(NULL)
    new_seed <- sample(1:2147483647, 1)
    set.seed(new_seed)
    
    # Initialize variables to avoid carryover between iterations
    Limb = NA
    Limb.Order = NA
    ID.Table = NA
    Data.Out = NA
    Data.Out.Include = NA
    LM.1 = NA
    LM.2 = NA
    LMM = NA
    
    # Create data for intervention group (Group 1)
    Group.1.L = data.frame(1:60 + 100, rep(1, 60), c(rep(1, N.Bilateral), rep(0, 60 - N.Bilateral)))
    colnames(Group.1.L) = c("ID", "Intervention", "Bilateral")
    Group.1.L$Limb = "L"
    
    Group.1.R = data.frame(1:60 + 100, rep(1, 60), c(rep(1, N.Bilateral), rep(0, 60 - N.Bilateral)))
    colnames(Group.1.R) = c("ID", "Intervention", "Bilateral")
    Group.1.R$Limb = "R"                    
    
    # Create data for control group (Group 2)
    Group.2.L = data.frame(1:60 + 200, rep(0, 60), c(rep(1, N.Bilateral), rep(0, 60 - N.Bilateral)))
    colnames(Group.2.L) = c("ID", "Intervention", "Bilateral")
    Group.2.L$Limb = "L"
    
    Group.2.R = data.frame(1:60 + 200, rep(0, 60), c(rep(1, N.Bilateral), rep(0, 60 - N.Bilateral)))
    colnames(Group.2.R) = c("ID", "Intervention", "Bilateral")
    Group.2.R$Limb = "R"
    
    # Combine all groups into a single data frame
    Data = rbind(Group.1.L, Group.1.R, Group.2.L, Group.2.R)
    
    # Create a limb table with random order
    Limb = data.frame(c(rep("L", 60), rep("R", 60)), rnorm(120))
    colnames(Limb) = c("Limb", "Random")
    Limb.Order = Limb[order(Limb[, "Random"]), ]
    
    # Create an ID table with random intercepts
    ID.Table = data.frame(unique(Data$ID), Limb.Order$Limb)
    colnames(ID.Table) = c("ID", "Limb.Include")
    ID.Table$Random.Int = rnorm(nrow(ID.Table), 0, u1)
    
    # Merge data with ID table
    Data.Out = merge(Data, ID.Table, by = "ID")
    
    # Determine which observations to include based on bilateral status and limb
    Data.Out$Include = ifelse(Data.Out$Bilateral == 1, 1,
                              ifelse(Data.Out$Limb == Data.Out$Limb.Include, 1, 0))
    
    Data.Out.Include = Data.Out[Data.Out$Include == 1, ]
    
    # Generate the true outcome variable with subject-specific and random effects
    Data.Out.Include$Y.True = B0 + B1 * Data.Out.Include$Intervention + B2 * Data.Out.Include$Bilateral +
      rnorm(nrow(Data.Out.Include), 0, E) + Data.Out.Include$Random.Int
    
    # Fit naive linear model
    LM.1 = lm(Y.True ~ Intervention, data = Data.Out.Include)
    
    # Fit random limb linear model
    Data.Out.Include.Random = Data.Out.Include[Data.Out.Include$Limb.Include == Data.Out.Include$Limb, ]
    LM.2 = lm(Y.True ~ Intervention, data = Data.Out.Include.Random)
    
    # Fit linear mixed model
    LMM.Try = try(lme(Y.True ~ Intervention, data = Data.Out.Include, random = ~1 | ID))
    
    # Handle potential errors in fitting the linear mixed model
    ifelse("try-error" %in% class(LMM.Try) | "try-error" %in% class(try(intervals(lme(Y.True ~ Intervention, data = Data.Out.Include, random = ~1 | ID))$fixed)),
           {LMM.Coef = matrix(NA, nrow = 2, ncol = 5)
           LMM.CI = matrix(NA, nrow = 2, ncol = 3)},
           {
             LMM = lme(Y.True ~ Intervention, data = Data.Out.Include, random = ~1 | ID)
             LMM.Coef = summary(lme(Y.True ~ Intervention, data = Data.Out.Include, random = ~1 | ID))$tTable
             LMM.CI = intervals(lme(Y.True ~ Intervention, data = Data.Out.Include, random = ~1 | ID))$fixed
           })
    
    # Safely compute intervals for LMM and capture errors
    safelme <- safely(.f = intervals, otherwise = NA, quiet = F)
    lmm.err <- lme(Y.True ~ Intervention, data = Data.Out.Include, random = ~1 | ID)
    lmm.int.err <- safelme(lmm.err)
    
    Error.Output[n, 1] <- n
    Error.Output[n, 2] <- ifelse(test = is.null(lmm.int.err$error[[1]]), yes = NA, no = lmm.int.err$error[[1]])
    
    # Compute within-subject correlation for bilateral subjects
    R = Data.Out.Include[Data.Out.Include$Bilateral == 1 & Data.Out.Include$Limb == "R", ]
    L = Data.Out.Include[Data.Out.Include$Bilateral == 1 & Data.Out.Include$Limb == "L", ]
    
    P.Corr = cor.test(x = R$Y.True, y = L$Y.True, method = 'pearson')
    B.2 = P.Corr$estimate
    
    # Store results for this iteration
    Simulation.Output[n, 1] = n
    # Naive model results
    Simulation.Output[n, 2] = summary(LM.1)$coefficients[2, 1]
    Simulation.Output[n, 3] = summary(LM.1)$coefficients[2, 2]
    Simulation.Output[n, 4] = summary(LM.1)$coefficients[2, 3]
    Simulation.Output[n, 5] = summary(LM.1)$coefficients[2, 4]
    # Random limb model results
    Simulation.Output[n, 6] = summary(LM.2)$coefficients[2, 1]
    Simulation.Output[n, 7] = summary(LM.2)$coefficients[2, 2]
    Simulation.Output[n, 8] = summary(LM.2)$coefficients[2, 3]
    Simulation.Output[n, 9] = summary(LM.2)$coefficients[2, 4]
    # Linear mixed model results
    Simulation.Output[n, 10] = LMM.Coef[2, 1]
    Simulation.Output[n, 11] = LMM.Coef[2, 2]
    Simulation.Output[n, 12] = LMM.Coef[2, 4]
    Simulation.Output[n, 13] = LMM.Coef[2, 5]
    # Within-subject correlation
    Simulation.Output[n, 14] = P.Corr$estimate
    # Standard deviations
    Simulation.Output[n, 15] = sd(Data.Out.Include.Random$Y.True)
    Simulation.Output[n, 16] = sd(Data.Out.Include.Random[Data.Out.Include.Random$Intervention == 1, "Y.True"])
    Simulation.Output[n, 17] = sd(Data.Out.Include.Random[Data.Out.Include.Random$Intervention == 0, "Y.True"])
    # Confidence intervals
    Simulation.Output[n, 18] = confint(LM.1)[2, 1]
    Simulation.Output[n, 19] = confint(LM.1)[2, 2]
    Simulation.Output[n, 20] = confint(LM.2)[2, 1]
    Simulation.Output[n, 21] = confint(LM.2)[2, 2]
    Simulation.Output[n, 22] = LMM.CI[2, 1]
    Simulation.Output[n, 23] = LMM.CI[2, 3]
    # Seed, N.Bilateral, and hypothesis
    Simulation.Output[n, 24] = new_seed
    Simulation.Output[n, 25] = N.Bilateral.List[j]
    Simulation.Output[n, 26] = "Null Hyp."
    
    # Convert simulation output to data frame
    Simulation.Output.Data = data.frame(Simulation.Output)
    
    # Define column names for the simulation output
    colnames(Simulation.Output.Data) = 
      c(
        'Iteration',
        'Naive.Est',
        'Naive.SE',
        'Naive.TValue',
        'Naive.PValue',
        'Random.Est', 
        'Random.SE', 
        'Random.TValue', 
        'Random.Pvalue',
        'LMM.Est', 
        'LMM.SE', 
        'LMM.TValue', 
        'LMM.Pvalue',
        'Correlation',
        'SD.All',
        'SD.Tx',
        'SD.NoTx',
        'Naive.LCL',
        'Naive.UCL',
        'Random.LCL',
        'Random.UCL',
        'LMM.LCL',
        'LMM.UCL',
        'Seed',
        'N.Bilateral',
        'Hypothesis'
      )
    
    # Convert error output to data frame
    Error.Output <- as.data.frame(Error.Output)
    
    # Define column names for the error output
    colnames(Error.Output) <- c("Iteration", "LME Error")
    
    # Save the results to files
    saveRDS(
      object = Simulation.Output.Data, 
      file = paste0("/Users/harrysmith/Documents/Carry_lab/misc/twoLimbSimulationStudy/Data_Out/Null/High_Corr_Data_Out_Null_", N.Bilateral.List[j], ".Rds")
    )
    saveRDS(
      object = Error.Output, 
      file = paste0("/Users/harrysmith/Documents/Carry_lab/misc/twoLimbSimulationStudy/Error_Out/Null/High_Corr_Error_Out_Null_", N.Bilateral.List[j], ".Rds")
    )
  }
}


rm(list=ls())   # Removes environment for next simulation

                                ############################ 
                                ## ALTERNATIVE HYPOTHESIS ## 
                                ############################ 

############################ 
## ALTERNATIVE HYPOTHESIS ## 
############################ 

# List of Number of Bilateral Limbs
N.Bilateral.List = as.matrix(
  x = c(60, 57, 54, 51, 48, 45, 42, 39, 36, 33, 30, 27, 24, 21, 18, 15, 12, 9, 6)
)

# Subject-specific intercept and error term
u1 = 16.75
E = 8.25

# Coefficients for the alternative hypothesis simulation
B0 = 5
B1 = 10
B2 = 0

# Number of replicates for the simulation (5 for testing, 10000 for actual simulation)
replicates = 10

# Create empty matrices to store simulation and error results
Simulation.Output = matrix(NA, nrow = replicates, ncol = 33)
Error.Output = matrix(NA, nrow = replicates, ncol = 3)

# Loop through each N.Bilateral value
for(j in 1:nrow(N.Bilateral.List)){
  
  # Number of bilateral subjects for this iteration
  N.Bilateral = N.Bilateral.List[j]
  
  # Inner loop for each replicate
  for (n in 1:replicates){
    
    # Reset seed to ensure reproducibility
    set.seed(NULL)
    new_seed <- sample(1:2147483647, 1)
    set.seed(new_seed)
    
    # Initialize variables to avoid carryover between iterations
    Limb = NA
    Limb.Order = NA
    ID.Table = NA
    Data.Out = NA
    Data.Out.Include = NA
    LM.1 = NA
    LM.2 = NA
    LMM = NA
    
    # Create data for intervention group (Group 1)
    Group.1.L = data.frame(1:60 + 100, rep(1, 60), c(rep(1, N.Bilateral), rep(0, 60 - N.Bilateral)))
    colnames(Group.1.L) = c("ID", "Intervention", "Bilateral")
    Group.1.L$Limb = "L"
    
    Group.1.R = data.frame(1:60 + 100, rep(1, 60), c(rep(1, N.Bilateral), rep(0, 60 - N.Bilateral)))
    colnames(Group.1.R) = c("ID", "Intervention", "Bilateral")
    Group.1.R$Limb = "R"                    
    
    # Create data for control group (Group 2)
    Group.2.L = data.frame(1:60 + 200, rep(0, 60), c(rep(1, N.Bilateral), rep(0, 60 - N.Bilateral)))
    colnames(Group.2.L) = c("ID", "Intervention", "Bilateral")
    Group.2.L$Limb = "L"
    
    Group.2.R = data.frame(1:60 + 200, rep(0, 60), c(rep(1, N.Bilateral), rep(0, 60 - N.Bilateral)))
    colnames(Group.2.R) = c("ID", "Intervention", "Bilateral")
    Group.2.R$Limb = "R"
    
    # Combine all groups into a single data frame
    Data = rbind(Group.1.L, Group.1.R, Group.2.L, Group.2.R)
    
    # Create a limb table with random order
    Limb = data.frame(c(rep("L", 60), rep("R", 60)), rnorm(120))
    colnames(Limb) = c("Limb", "Random")
    Limb.Order = Limb[order(Limb[, "Random"]), ]
    
    # Create an ID table with random intercepts
    ID.Table = data.frame(unique(Data$ID), Limb.Order$Limb)
    colnames(ID.Table) = c("ID", "Limb.Include")
    ID.Table$Random.Int = rnorm(nrow(ID.Table), 0, u1)
    
    # Merge data with ID table
    Data.Out = merge(Data, ID.Table, by = "ID")
    
    # Determine which observations to include based on bilateral status and limb
    Data.Out$Include = ifelse(Data.Out$Bilateral == 1, 1,
                              ifelse(Data.Out$Limb == Data.Out$Limb.Include, 1, 0))
    
    Data.Out.Include = Data.Out[Data.Out$Include == 1, ]
    
    # Generate the true outcome variable with subject-specific and random effects
    Data.Out.Include$Y.True = B0 + B1 * Data.Out.Include$Intervention + B2 * Data.Out.Include$Bilateral +
      rnorm(nrow(Data.Out.Include), 0, E) + Data.Out.Include$Random.Int
    
    # Fit naive linear model
    LM.1 = lm(Y.True ~ Intervention, data = Data.Out.Include)
    
    # Fit random limb linear model
    Data.Out.Include.Random = Data.Out.Include[Data.Out.Include$Limb.Include == Data.Out.Include$Limb, ]
    LM.2 = lm(Y.True ~ Intervention, data = Data.Out.Include.Random)
    
    # Fit linear mixed model
    LMM.Try = try(lme(Y.True ~ Intervention, data = Data.Out.Include, random = ~1 | ID))
    
    # Handle potential errors in fitting the linear mixed model
    ifelse("try-error" %in% class(LMM.Try) | "try-error" %in% class(try(intervals(lme(Y.True ~ Intervention, data = Data.Out.Include, random = ~1 | ID))$fixed)),
           {LMM.Coef = matrix(NA, nrow = 2, ncol = 5)
           LMM.CI = matrix(NA, nrow = 2, ncol = 3)},
           {
             LMM = lme(Y.True ~ Intervention, data = Data.Out.Include, random = ~1 | ID)
             LMM.Coef = summary(lme(Y.True ~ Intervention, data = Data.Out.Include, random = ~1 | ID))$tTable
             LMM.CI = intervals(lme(Y.True ~ Intervention, data = Data.Out.Include, random = ~1 | ID))$fixed
           })
    
    # Safely compute intervals for LMM and capture errors
    safelme <- safely(.f = intervals, otherwise = NA, quiet = F)
    lmm.err <- lme(Y.True ~ Intervention, data = Data.Out.Include, random = ~1 | ID)
    lmm.int.err <- safelme(lmm.err)
    
    Error.Output[n, 1] <- n
    Error.Output[n, 2] <- ifelse(test = is.null(lmm.int.err$error[[1]]), yes = NA, no = lmm.int.err$error[[1]])
    
    # Compute within-subject correlation for bilateral subjects
    R = Data.Out.Include[Data.Out.Include$Bilateral == 1 & Data.Out.Include$Limb == "R", ]
    L = Data.Out.Include[Data.Out.Include$Bilateral == 1 & Data.Out.Include$Limb == "L", ]
    
    P.Corr = cor.test(x = R$Y.True, y = L$Y.True, method = 'pearson')
    B.2 = P.Corr$estimate
    
    # Store results for this iteration
    Simulation.Output[n, 1] = n
    # Naive model results
    Simulation.Output[n, 2] = summary(LM.1)$coefficients[2, 1]
    Simulation.Output[n, 3] = summary(LM.1)$coefficients[2, 2]
    Simulation.Output[n, 4] = summary(LM.1)$coefficients[2, 3]
    Simulation.Output[n, 5] = summary(LM.1)$coefficients[2, 4]
    # Random limb model results
    Simulation.Output[n, 6] = summary(LM.2)$coefficients[2, 1]
    Simulation.Output[n, 7] = summary(LM.2)$coefficients[2, 2]
    Simulation.Output[n, 8] = summary(LM.2)$coefficients[2, 3]
    Simulation.Output[n, 9] = summary(LM.2)$coefficients[2, 4]
    # Linear mixed model results
    Simulation.Output[n, 10] = LMM.Coef[2, 1]
    Simulation.Output[n, 11] = LMM.Coef[2, 2]
    Simulation.Output[n, 12] = LMM.Coef[2, 4]
    Simulation.Output[n, 13] = LMM.Coef[2, 5]
    # Within-subject correlation
    Simulation.Output[n, 18] = P.Corr$estimate
    # Standard deviations
    Simulation.Output[n, 19] = sd(Data.Out.Include.Random$Y.True)
    Simulation.Output[n, 20] = sd(Data.Out.Include.Random[Data.Out.Include.Random$Intervention == 1, "Y.True"])
    Simulation.Output[n, 21] = sd(Data.Out.Include.Random[Data.Out.Include.Random$Intervention == 0, "Y.True"])
    # Confidence intervals
    Simulation.Output[n, 22] = confint(LM.1)[2, 1]
    Simulation.Output[n, 23] = confint(LM.1)[2, 2]
    Simulation.Output[n, 24] = confint(LM.2)[2, 1]
    Simulation.Output[n, 25] = confint(LM.2)[2, 2]
    Simulation.Output[n, 26] = LMM.CI[2, 1]
    Simulation.Output[n, 27] = LMM.CI[2, 3]
    # Seed, N.Bilateral, and hypothesis
    Simulation.Output[n, 28] = new_seed
    Simulation.Output[n, 29] = N.Bilateral.List[j]
    Simulation.Output[n, 30] = "Alt. Hyp."
    
    # Convert simulation output to data frame
    Simulation.Output.Data = data.frame(Simulation.Output)
    
    # Define column names for the simulation output
    colnames(Simulation.Output.Data) = 
      c(
        'Iteration',
        'Naive.Est',
        'Naive.SE',
        'Naive.TValue',
        'Naive.PValue',
        'Random.Est', 
        'Random.SE', 
        'Random.TValue', 
        'Random.Pvalue',
        'LMM.Est', 
        'LMM.SE', 
        'LMM.TValue', 
        'LMM.Pvalue',
        'Correlation',
        'SD.All',
        'SD.Tx',
        'SD.NoTx',
        'Naive.LCL',
        'Naive.UCL',
        'Random.LCL',
        'Random.UCL',
        'LMM.LCL',
        'LMM.UCL',
        'Seed',
        'N.Bilateral',
        'Hypothesis'
      )
    
    # Convert error output to data frame
    Error.Output <- as.data.frame(Error.Output)
    
    # Define column names for the error output
    colnames(Error.Output) <- c(
      "Iteration", 
      "LME Error"
    )
    
    # Save the results to files
    saveRDS(
      object = Simulation.Output.Data, 
      file = paste0("/Users/harrysmith/Documents/Carry_lab/misc/twoLimbSimulationStudy/Data_Out/Alternative/High_Corr_Data_Out_Alternative_", N.Bilateral.List[j], ".Rds")
    )
    
    saveRDS(
      object = Error.Output, 
      file = paste0("/Users/harrysmith/Documents/Carry_lab/misc/twoLimbSimulationStudy/Error_Out/Alternative/High_Corr_Error_Out_Alternative_", N.Bilateral.List[j], ".Rds")
    )
  }
}


rm(list = ls())

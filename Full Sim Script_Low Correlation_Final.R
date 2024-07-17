rm(list = ls())

# Libraries 
library(nlme)
library(purrr)



                                #####################
                                ## NULL HYPOTHESIS ## 
                                #####################



N.Bilateral.List = as.matrix(
  x = c(60, 57, 54, 51, 48, 45, 42, 39, 36, 33, 30, 27, 24, 21, 18, 15, 12, 9, 6)
  )

#Subject specific intercept
u1 = 10
E = 16

B0 = 5
B1 = 0
B2 = 0

#Start of simulation code

replicates = 10000 #5 for testing, 10000 for actual simulation

# Simulation and Error Dataframes (empty)
Simulation.Output = matrix(NA, nrow = replicates, ncol = 33)
Error.Output = matrix(NA, nrow = replicates, ncol = 3)

# Simulation (NULL)

for(j in 1:nrow(N.Bilateral.List)){
 
  # Iterates through different N.Bilateral
  N.Bilateral = N.Bilateral.List[j]
 
   for (n in 1:replicates){
    
    #Resetting seed between iterations
    set.seed( NULL )
    
    #Setting seed so that it is reproducible
    new_seed <- sample(1:2147483647, 1)
    
    set.seed( new_seed )
    
    #Resetting to null so there is not carryover between simulation iterations
    Limb = NA
    Limb.Order = NA
    ID.Table = NA
    Data.Out = NA
    Data.Out.Include = NA
    LM.1 = NA
    LM.2 = NA
    LMM = NA

    
    #Intervention
    Group.1.L = data.frame(1:60 + 100, rep(1, 60), c( rep(1, N.Bilateral), rep(0, 60 - N.Bilateral) ))
    colnames(Group.1.L) = c("ID", "Intervention", "Bilateral")
    Group.1.L$Limb = "L"
    
    Group.1.R = data.frame(1:60+100, rep(1, 60), c( rep(1, N.Bilateral), rep(0, 60-N.Bilateral) ))
    colnames(Group.1.R) = c("ID", "Intervention", "Bilateral")
    Group.1.R$Limb = "R"                    

    #Control
    Group.2.L = data.frame(1:60+200, rep(0, 60), c( rep(1, N.Bilateral), rep(0, 60-N.Bilateral) ))
    colnames(Group.2.L) = c("ID", "Intervention", "Bilateral")
    Group.2.L$Limb = "L"
    
    Group.2.R = data.frame(1:60+200, rep(0, 60), c( rep(1, N.Bilateral), rep(0, 60-N.Bilateral) ))
    colnames(Group.2.R) = c("ID", "Intervention", "Bilateral")
    Group.2.R$Limb = "R"     
    
    Data = rbind(Group.1.L, Group.1.R, Group.2.L, Group.2.R)
  
    #Limb table
    Limb = data.frame( c(rep("L", 60), rep("R", 60)), rnorm( 120 ) )
    colnames(Limb) = c("Limb",  "Random")
    Limb.Order = Limb[order(Limb[, "Random"]), ]
    
    ID.Table = data.frame( unique(Data$ID), Limb.Order$Limb)
    colnames(ID.Table) = c("ID", "Limb.Include")
    ID.Table$Random.Int = rnorm(nrow(ID.Table), 0, u1)
    
    Data.Out = merge(Data, ID.Table, by = "ID")
    
    Data.Out$Include = ifelse(Data.Out$Bilateral==1, 1,
                              ifelse(Data.Out$Limb==Data.Out$Limb.Include, 1, 0))
    
    Data.Out.Include = Data.Out[Data.Out$Include==1, ]
    
    #Creating final dataset
    Data.Out.Include$Y.True = B0 + B1*Data.Out.Include$Intervention + B2*Data.Out.Include$Bilateral + rnorm(nrow(Data.Out.Include), 0, E) + Data.Out.Include$Random.Int
    
    #Naive Model
    LM.1 = lm(Y.True ~ Intervention, data = Data.Out.Include)
    
    #Random limb condition
    Data.Out.Include.Random = Data.Out.Include[Data.Out.Include$Limb.Include==Data.Out.Include$Limb, ]
    
    LM.2 = lm(Y.True ~ Intervention, data = Data.Out.Include.Random)
    
    #Linear mixed model condition
    
    LMM.Try = try(lme(Y.True ~ Intervention, data = Data.Out.Include, random=~1|ID))
    
    ifelse("try-error" %in% class(LMM.Try) | "try-error" %in% class(try(intervals(lme(Y.True ~ Intervention, data = Data.Out.Include, random=~1|ID))$fixed )), 
           {LMM.Coef = matrix(NA, nrow = 2, ncol = 5)
           LMM.CI = matrix(NA, nrow = 2, ncol = 3)},
           {
             LMM = lme(Y.True ~ Intervention, data = Data.Out.Include, random=~1|ID)
             LMM.Coef = summary(lme(Y.True ~ Intervention, data = Data.Out.Include, random=~1|ID))$tTable
             LMM.CI = intervals(lme(Y.True ~ Intervention, data = Data.Out.Include, random=~1|ID))$fixed})
  
    ## Catch Error LME 
    safelme <- safely(
      .f = intervals, 
      otherwise = NA, 
      quiet = F
    )
    
    lmm.err <- lme(Y.True ~ Intervention, data = Data.Out.Include, random = ~1|ID)
    lmm.int.err <- safelme(lmm.err)
    
    Error.Output[n, 1] <- n
    Error.Output[n, 2] <- ifelse(
      test = is.null(lmm.int.err$error[[1]]),
      yes = NA,
      no = lmm.int.err$error[[1]]
    )
    
    
    #Creating dataset to review within subject correlation
    R = Data.Out.Include[Data.Out.Include$Bilateral==1 & Data.Out.Include$Limb == "R", ]
    L = Data.Out.Include[Data.Out.Include$Bilateral==1 & Data.Out.Include$Limb == "L", ]
    
    P.Corr = cor.test(x=R$Y.True, y=L$Y.True, method = 'pearson')
    
    B.2 = P.Corr$estimate
  
    #Output, by simulation iteration
    Simulation.Output[n, 1] = n
    #Naive
    Simulation.Output[n, 2] = summary(LM.1)$coefficients[2, 1]
    Simulation.Output[n, 3] = summary(LM.1)$coefficients[2, 2]
    Simulation.Output[n, 4] = summary(LM.1)$coefficients[2, 3]
    Simulation.Output[n, 5] = summary(LM.1)$coefficients[2, 4]
    #Random
    Simulation.Output[n, 6] = summary(LM.2)$coefficients[2, 1]
    Simulation.Output[n, 7] = summary(LM.2)$coefficients[2, 2]
    Simulation.Output[n, 8] = summary(LM.2)$coefficients[2, 3]
    Simulation.Output[n, 9] = summary(LM.2)$coefficients[2, 4]
    #LMM
    Simulation.Output[n, 10] = LMM.Coef[2, 1]
    Simulation.Output[n, 11] = LMM.Coef[2, 2]
    Simulation.Output[n, 12] = LMM.Coef[2, 4]
    Simulation.Output[n, 13] = LMM.Coef[2, 5]

    #Within subject correlation
    Simulation.Output[n, 14] = P.Corr$estimate
    #SD Overall
    Simulation.Output[n, 15] = sd(Data.Out.Include.Random$Y.True)
    #SD Tx
    Simulation.Output[n, 16] = sd(Data.Out.Include.Random[Data.Out.Include.Random$Intervention==1, "Y.True"])
    #SD Control
    Simulation.Output[n, 17] = sd(Data.Out.Include.Random[Data.Out.Include.Random$Intervention==0, "Y.True"])
    
    #LM confidence intervals
    Simulation.Output[n, 18] = confint(LM.1)[2,1]
    Simulation.Output[n, 19] = confint(LM.1)[2,2]
    Simulation.Output[n, 20] = confint(LM.2)[2,1]
    Simulation.Output[n, 21] = confint(LM.2)[2,2]
    #LMM confidence intervals
    Simulation.Output[n, 22] = LMM.CI[2, 1]
    Simulation.Output[n, 23] = LMM.CI[2, 3]
    
    Simulation.Output[n, 24] = new_seed
    Simulation.Output[n, 25] = N.Bilateral.List[j]
    Simulation.Output[n, 26] = "Null Hyp."
    
    Simulation.Output.Data = data.frame(Simulation.Output)
    
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
    
    Error.Output <- as.data.frame(Error.Output)
    
    colnames(Error.Output) <- c(
      "Iteration", 
      "LME Error"
      )
    
    saveRDS(
      object = Simulation.Output.Data, 
      file = paste0("RDS Files/Data Out/Null/Low_Corr_Data_Out_Null_", N.Bilateral.List[j],".Rds")
      )
    saveRDS(
      object = Error.Output, 
      file = paste0("RDS Files/Error Out/Null/Low_Corr_Error_Out_Null_", N.Bilateral.List[j],".Rds")
      )
  }
}

rm(list=ls())   # Removes environment for next simulation

                                ############################ 
                                ## ALTERNATIVE HYPOTHESIS ## 
                                ############################ 

# List of Number of Bilateral Limbs
N.Bilateral.List = as.matrix(
  x = c(60, 57, 54, 51, 48, 45, 42, 39, 36, 33, 30, 27, 24, 21, 18, 15, 12, 9, 6)
)

#Subject specific intercept
u1 = 10
E = 16

B0 = 5
B1 = 10
B2 = 0

#Start of simulation code

replicates = 10000

Simulation.Output = matrix(NA, nrow = replicates, ncol = 33)
Error.Output = matrix(NA, nrow = replicates, ncol = 3)

for(j in 1:nrow(N.Bilateral.List)){
  
  N.Bilateral = N.Bilateral.List[j]
  
  for (n in 1:replicates){
    
    #Resetting seed between iterations
    set.seed( NULL )
    
    #Setting seed so that it is reproducible
    new_seed <- sample(1:2147483647, 1)
    
    set.seed( new_seed )
    
    #Resetting to null so there is not carryover between simulation iterations
    Limb = NA
    Limb.Order = NA
    ID.Table = NA
    Data.Out = NA
    Data.Out.Include = NA
    LM.1 = NA
    LM.2 = NA
    LMM = NA
    
    #Intervention
    Group.1.L = data.frame(1:60 + 100, rep(1, 60), c( rep(1, N.Bilateral), rep(0, 60 - N.Bilateral) ))
    colnames(Group.1.L) = c("ID", "Intervention", "Bilateral")
    Group.1.L$Limb = "L"
    
    Group.1.R = data.frame(1:60+100, rep(1, 60), c( rep(1, N.Bilateral), rep(0, 60-N.Bilateral) ))
    colnames(Group.1.R) = c("ID", "Intervention", "Bilateral")
    Group.1.R$Limb = "R"                    
    
    #Control
    Group.2.L = data.frame(1:60+200, rep(0, 60), c( rep(1, N.Bilateral), rep(0, 60-N.Bilateral) ))
    colnames(Group.2.L) = c("ID", "Intervention", "Bilateral")
    Group.2.L$Limb = "L"
    
    Group.2.R = data.frame(1:60+200, rep(0, 60), c( rep(1, N.Bilateral), rep(0, 60-N.Bilateral) ))
    colnames(Group.2.R) = c("ID", "Intervention", "Bilateral")
    Group.2.R$Limb = "R"     
    
    Data = rbind(Group.1.L, Group.1.R, Group.2.L, Group.2.R)
    
    #Limb table
    Limb = data.frame( c(rep("L", 60), rep("R", 60)), rnorm( 120 ) )
    colnames(Limb) = c("Limb",  "Random")
    Limb.Order = Limb[order(Limb[, "Random"]), ]
    
    ID.Table = data.frame( unique(Data$ID), Limb.Order$Limb)
    colnames(ID.Table) = c("ID", "Limb.Include")
    ID.Table$Random.Int = rnorm(nrow(ID.Table), 0, u1)
    
    Data.Out = merge(Data, ID.Table, by = "ID")
    
    Data.Out$Include = ifelse(Data.Out$Bilateral==1, 1,
                              ifelse(Data.Out$Limb==Data.Out$Limb.Include, 1, 0))
    
    Data.Out.Include = Data.Out[Data.Out$Include==1, ]
    
    #Creating final dataset
    Data.Out.Include$Y.True = B0 + B1*Data.Out.Include$Intervention + B2*Data.Out.Include$Bilateral + rnorm(nrow(Data.Out.Include), 0, E) + Data.Out.Include$Random.Int
    
    #Naive Model
    LM.1 = lm(Y.True ~ Intervention, data = Data.Out.Include)
    
    #Random limb condition
    Data.Out.Include.Random = Data.Out.Include[Data.Out.Include$Limb.Include==Data.Out.Include$Limb, ]
    
    LM.2 = lm(Y.True ~ Intervention, data = Data.Out.Include.Random)
    
    #Linear mixed model condition
    
    LMM.Try = try(lme(Y.True ~ Intervention, data = Data.Out.Include, random=~1|ID))
    
    ifelse("try-error" %in% class(LMM.Try) | "try-error" %in% class(try(intervals(lme(Y.True ~ Intervention, data = Data.Out.Include, random=~1|ID))$fixed )), 
           {LMM.Coef = matrix(NA, nrow = 2, ncol = 5)
           LMM.CI = matrix(NA, nrow = 2, ncol = 3)},
           {
             LMM = lme(Y.True ~ Intervention, data = Data.Out.Include, random=~1|ID)
             LMM.Coef = summary(lme(Y.True ~ Intervention, data = Data.Out.Include, random=~1|ID))$tTable
             LMM.CI = intervals(lme(Y.True ~ Intervention, data = Data.Out.Include, random=~1|ID))$fixed})
    
    
    ## Catch Error LME 
    safelme <- safely(
      .f = intervals, 
      otherwise = NA, 
      quiet = F
    )
    
    lmm.err <- lme(Y.True ~ Intervention, data = Data.Out.Include, random = ~1|ID)
    lmm.int.err <- safelme(lmm.err)
    
    Error.Output[n, 1] <- n
    Error.Output[n, 2] <- ifelse(
      test = is.null(lmm.int.err$error[[1]]),
      yes = NA,
      no = lmm.int.err$error[[1]]
    )
    
    
    Error.Output[n, 1] <- n
    Error.Output[n, 2] <- ifelse(
      test = is.null(lmm.int.err$error[[1]]),
      yes = NA,
      no = lmm.int.err$error[[1]]
    )
    
    
    
    #Creating dataset to review within subject correlation
    R = Data.Out.Include[Data.Out.Include$Bilateral==1 & Data.Out.Include$Limb == "R", ]
    L = Data.Out.Include[Data.Out.Include$Bilateral==1 & Data.Out.Include$Limb == "L", ]
    
    P.Corr = cor.test(x=R$Y.True, y=L$Y.True, method = 'pearson')
    
    B.2 = P.Corr$estimate
    
    #Output, by simulation iteration
    Simulation.Output[n, 1] = n
    #Naive
    Simulation.Output[n, 2] = summary(LM.1)$coefficients[2, 1]
    Simulation.Output[n, 3] = summary(LM.1)$coefficients[2, 2]
    Simulation.Output[n, 4] = summary(LM.1)$coefficients[2, 3]
    Simulation.Output[n, 5] = summary(LM.1)$coefficients[2, 4]
    #Random
    Simulation.Output[n, 6] = summary(LM.2)$coefficients[2, 1]
    Simulation.Output[n, 7] = summary(LM.2)$coefficients[2, 2]
    Simulation.Output[n, 8] = summary(LM.2)$coefficients[2, 3]
    Simulation.Output[n, 9] = summary(LM.2)$coefficients[2, 4]
    #LMM
    Simulation.Output[n, 10] = LMM.Coef[2, 1]
    Simulation.Output[n, 11] = LMM.Coef[2, 2]
    Simulation.Output[n, 12] = LMM.Coef[2, 4]
    Simulation.Output[n, 13] = LMM.Coef[2, 5]

    #Within subject correlation
    Simulation.Output[n, 18] = P.Corr$estimate
    #SD Overall
    Simulation.Output[n, 19] = sd(Data.Out.Include.Random$Y.True)
    #SD Tx
    Simulation.Output[n, 20] = sd(Data.Out.Include.Random[Data.Out.Include.Random$Intervention==1, "Y.True"])
    #SD Control
    Simulation.Output[n, 21] = sd(Data.Out.Include.Random[Data.Out.Include.Random$Intervention==0, "Y.True"])
    
    #LM confidence intervals
    Simulation.Output[n, 22] = confint(LM.1)[2,1]
    Simulation.Output[n, 23] = confint(LM.1)[2,2]
    Simulation.Output[n, 24] = confint(LM.2)[2,1]
    Simulation.Output[n, 25] = confint(LM.2)[2,2]
    #LMM confidence intervals
    Simulation.Output[n, 26] = LMM.CI[2, 1]
    Simulation.Output[n, 27] = LMM.CI[2, 3]
    
    Simulation.Output[n, 28] = new_seed
    Simulation.Output[n, 29] = N.Bilateral.List[j]
    Simulation.Output[n, 30] = "Alt. Hyp."
    
    Simulation.Output.Data = data.frame(Simulation.Output)
    
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
    
    Error.Output <- as.data.frame(Error.Output)
    
    colnames(Error.Output) <- c(
      "Iteration", 
      "LME Error"
    )
    
    saveRDS(
      object = Simulation.Output.Data, 
      file = paste0("RDS Files/Data Out/Alternative/Low_Corr_Data_Out_Alternative_", N.Bilateral.List[j],".Rds")
    )
    
    saveRDS(
      object = Error.Output, 
      file = paste0("RDS Files/Error Out/Alternative/Low_Corr_Error_Out_Alternative_", N.Bilateral.List[j],".Rds")
    )
  }
}

rm(list = ls())

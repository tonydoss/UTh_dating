# UTh_dating - Closed-system 230Th-U age calculations
#
#    Copyright 2017 Anthony Dosseto
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
# http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


library(deSolve)
library(ggplot2)

l234 <- 2.8262e-6 # 234U decay constant (a-1)
l230 <- 9.1577e-6 # 230Th decay constant (a-1)

# set working directory
#path <- "C:/Users/me/mydatafolder"
#setwd(path)

# name of sample to solve
sample_name <- "mysample"
sample_name <- "MK16"
# import iolite results
iolite_results <- read.table("IoliteExport_All_Integrations.txt", 
                             header = TRUE, sep = "\t", comment.char = "")
# create dataframe with data only for samples to solve
data <- subset(iolite_results, (grepl(sample_name, X)))
# data <- data[-c(8,11),]
# number of samples to solve
number_sampletosolve <- nrow(data)
# nb of times optimisation is repeated (for each sample)
nbit <- 100

lowerbound <- c(2, 1.0) # lower bound values for age (log10(yr)) and initial (234U/238U)
upperbound <- c(6, 5.0) # upper bound values for age (log10(yr)) and initial (234U/238U)

# create vectors
time_results <- vector(mode="numeric", length=number_sampletosolve)
err_time_results <- vector(mode="numeric", length=number_sampletosolve)
R48i_results <- vector(mode="numeric", length=number_sampletosolve)
err_R48i_results <- vector(mode="numeric", length=number_sampletosolve)
time_2se_results <- vector(mode="numeric", length=number_sampletosolve)
R48i_2se_results <- vector(mode="numeric", length=number_sampletosolve)

# repeat loop for each sample
for (count in 1:number_sampletosolve){
  U48meas <- data$U234_U238_CORR[count]
  Th0U8meas <- data$Th230_U238_CORR[count]
  err_R08 <- data$Th230_U238_CORR_Int2SE[count]
  err_R48 <- data$U234_U238_CORR_Int2SE[count]
  
  # create list for optimisation results
  sol=list()
  # create vectors
  U48calc <- vector(mode="numeric", length=nbit)
  Th0U8calc <- vector(mode="numeric", length=nbit)
  time <- vector(mode="numeric", length=nbit)
  R48i <- vector(mode="numeric", length=nbit)
  U48calc_2se <- vector(mode="numeric", length=nbit)
  Th0U8calc_2se <- vector(mode="numeric", length=nbit)
  time_2se <- vector(mode="numeric", length=nbit)
  R48i_2se <- vector(mode="numeric", length=nbit)
  
  
  # repeat optimisation 'nbit' number of times for a given sample
  for (i in 1:nbit){
    # start optimisation with random age and initial (234U/23U) taken from the range of values allowed
    init_time <- runif(1, lowerbound[1], upperbound[1])
    init_R48i <- runif(1, lowerbound[2], upperbound[2])
    paraminit <- c(init_time, init_R48i)
    
    # function to minimise
    funmin <- function(x) {
      t <- 10^x[1] # age in yr
      U48i <- x[2] # intial (234U/238U)
      
      U48calc <- 1 + (U48i - 1)*exp(-l234*t) # (234U/238U)
      Th0U8calc <- 1 - exp(-l230*t) + (U48calc - 1)*(l230/(l230 - l234))*(1 - exp((l234 - l230)*t))
      
      fmin <- sum((U48calc - U48target)^2 + (Th0U8calc - Th0U8target)^2 ) # function to minimise
    }
    
    # optimisations
    # do optimisation with measured (234U/238U) & (230Th/238U) 
    U48target <- runif(1, U48meas - err_R48, U48meas + err_R48)
    Th0U8target <- runif(1, Th0U8meas - err_R08, Th0U8meas + err_R08)
    sol <- optim(paraminit, funmin, method = "L-BFGS-B",
                 lower = lowerbound, upper = upperbound, control = list(factr = 1e-8))
    
    # store calculated age, initial (234U/23U) and calculated activity ratios for each optimisation
    time[i] <- 10^sol$par[1]
    R48i[i] <- sol$par[2]
    U48calc[i] <- 1 + (R48i[i] - 1)*exp(-l234*time[i]) # (234U/238U)
    Th0U8calc[i] <- 1 - exp(-l230*time[i]) - (U48calc[i] - 1)*(l230/(l234 - l230))*(1 - exp((l234 - l230)*time[i]))
  }
  
  # store results from all optimisations for a given sample
  results <- as.data.frame(cbind(time, R48i, 
                                 U48calc, Th0U8calc))
  # take the median of all ages and initial (234U/23U)
  median_time <- median(results$time)
  time_2se <- 2*sd(results$time)
  median_R48i <- median(results$R48i)
  R48i_2se <- 2*sd(results$R48i)
  
  # store age, error on age and initial (234U/23U) for each sample
  time_results[count] <- median_time
  time_2se_results[count] <- time_2se
  R48i_results[count] <- median_R48i
  R48i_2se_results[count] <- R48i_2se
}

# store results
final_results <- as.data.frame(cbind(round(time_results/1000,3), 
                                     round(time_2se_results/1000,3), 
                                     round(R48i_results,3), round(R48i_2se_results,3)))
final_results <- as.data.frame(cbind(as.character(data$X[1:number_sampletosolve]), final_results))

# give column names
colnames(final_results) <- c("ID", "Age (kyr)", "2sd", "(234U/238U)i", "2sd")

# display results in the console
print(final_results)
# display average age (in ka)
mean(time_results/1000)

# save results to a csv file with 'sample_name' as file name
write.table(final_results, file = paste(sample_name,".csv"), sep = ",", row.names = F)

# plot ages
p1 <- ggplot(final_results, aes(ID, `Age (kyr)`)) + # plot ages
  geom_errorbar(aes(ymin = (`Age (kyr)` - `2sd`),ymax = (`Age (kyr)` + `2sd`)), width=0.1) + # plot error bars
  geom_point(size=5) + # plot points
  xlab("Sample ID") + # x axis label
  ylab("Age (ka)") # y axis label

p1

# change column name of initial (234U/238U) error so it can be used to show error bars
colnames(final_results) <- c("ID", "Age (kyr)", "2sd", "(234U/238U)i", "2sd#2")

# plot initial (234U/238U)
p2 <- ggplot(final_results, aes(ID, `(234U/238U)i`)) + # plot ages
  geom_errorbar(aes(ymin = (`(234U/238U)i` - `2sd#2`),ymax = (`(234U/238U)i` + `2sd#2`)), width=0.1) + # plot error bars
  geom_point(size=5) + # plot points
  xlab("Sample ID") + # x axis label
  ylab("Age (ka)") # y axis label

p2

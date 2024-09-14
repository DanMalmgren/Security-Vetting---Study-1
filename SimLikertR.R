#Generate likert scale data

library(LikertMakeR) #https://github.com/WinzarH/LikertMakeR
library(reshape2)
library(tibble)
library(tidyr) 
library(ggplot2)



######################## Start simulation of data ##############################

#ensure repicable results
set.seed(5476)

#For Sweden

#TA = parental
TA <- lfast(
  n = 30,
  mean = 6,
  sd = 1.5,
  lowerbound = 1,
  upperbound = 7,
  items = 3
)

#TB = married
TB <- lfast(
  n = 30,
  mean = 6,
  sd = 1,
  lowerbound = 1,
  upperbound = 7,
  items = 3
)

#TC = born
TC <- lfast(
  n = 30,
  mean = 6,
  sd = 1.5,
  lowerbound = 1,
  upperbound = 7,
  items = 3
)

Sw <- data.frame(TA, TB, TC)
rm(TA, TB, TC)

#For UK

#TA = parental
TA <- lfast(
  n = 30,
  mean = 5.5,
  sd = 1.5,
  lowerbound = 1,
  upperbound = 7,
  items = 3
)

#TB = married
TB <- lfast(
  n = 30,
  mean = 5.5,
  sd = 1,
  lowerbound = 1,
  upperbound = 7,
  items = 3
)

#TC = born
TC <- lfast(
  n = 30,
  mean = 5.5,
  sd = 1.5,
  lowerbound = 1,
  upperbound = 7,
  items = 3
)

Uk <- data.frame(TA, TB, TC)
rm(TA, TB, TC)

#For Germany

#TA = parental
TA <- lfast(
  n = 30,
  mean = 5.0,
  sd = 1.5,
  lowerbound = 1,
  upperbound = 7,
  items = 3
)

#TB = married
TB <- lfast(
  n = 30,
  mean = 5.0,
  sd = 1,
  lowerbound = 1,
  upperbound = 7,
  items = 3
)

#TC = born
TC <- lfast(
  n = 30,
  mean = 5.0,
  sd = 0.5,
  lowerbound = 1,
  upperbound = 7,
  items = 3
)

Ge <- data.frame(TA, TB, TC)
rm(TA, TB, TC)

#For Brazil

#TA = parental
TA <- lfast(
  n = 30,
  mean = 4.0,
  sd = 1.5,
  lowerbound = 1,
  upperbound = 7,
  items = 3
)

#TB = married
TB <- lfast(
  n = 30,
  mean = 4.0,
  sd = 1,
  lowerbound = 1,
  upperbound = 7,
  items = 3
)

#TC = born
TC <- lfast(
  n = 30,
  mean = 4.0,
  sd = 0.5,
  lowerbound = 1,
  upperbound = 7,
  items = 3
)

Br <- data.frame(TA, TB, TC)
rm(TA, TB, TC)

#Thailand

#TA = parental
TA <- lfast(
  n = 30,
  mean = 3.0,
  sd = 1.5,
  lowerbound = 1,
  upperbound = 7,
  items = 3
)

#TB = married
TB <- lfast(
  n = 30,
  mean = 3.0,
  sd = 1,
  lowerbound = 1,
  upperbound = 7,
  items = 3
)

#TC = born
TC <- lfast(
  n = 30,
  mean = 3.0,
  sd = 0.5,
  lowerbound = 1,
  upperbound = 7,
  items = 3
)

Th <- data.frame(TA, TB, TC)
rm(TA, TB, TC)

######################### End simulation of data ##############################



######################### Create the data frame ############################# 

#Convert each frame to long format
#CI = Thailand, CJ = Brazil, CK = Germany, CL = UK, CM = Sweden
CM <- gather(Sw, key = "T", value = "Y", c("TA", "TB", "TC" ))
CL <- gather(Uk, key = "T", value = "Y", c("TA", "TB", "TC" ))
CK <- gather(Ge, key = "T", value = "Y", c("TA", "TB", "TC" ))
CJ <- gather(Br, key = "T", value = "Y", c("TA", "TB", "TC" ))
CI <- gather(Th, key = "T", value = "Y", c("TA", "TB", "TC" ))

#Remove the short frames from the environment
rm(Sw, Uk, Ge, Br, Th)

#Create a variable for factor for country in each data frame

C <- rep(c("CI"), len = 90)
CI <- cbind(C, CI)
rm(C)

C <- rep(c("CJ"), len = 90)
CJ <- cbind(C, CJ)
rm(C)

C <- rep(c("CK"), len = 90)
CK <- cbind(C, CK)
rm(C)

C <- rep(c("CL"), len = 90)
CL <- cbind(C, CL)
rm(C)

C <- rep(c("CM"), len = 90)
CM <- cbind(C, CM)
rm(C)


#Create subj id
id1 <- rep(1:30)
id2 <- rep(31:60)
id3 <- rep(61:90)

id <- data.frame(id1, id2, id3)

idL <- gather(id, key = "id", value = "id", (c(id1, id2, id3)))

idL <- subset(idL, select = -c(id))

#Insert subject id
CI <- cbind(idL, CI)
CJ <- cbind(idL, CJ)
CK <- cbind(idL, CK)
CL <- cbind(idL, CL)
CM <- cbind(idL, CM)

rm(id1, id2, id3, id, idL)

#Creating item
Item <- rep(1:6, each = 30, len = 90)
CI <- cbind(Item, CI)
rm(Item)

Item <- rep(7:12, each = 30, len = 90)
CJ <- cbind(Item, CJ)
rm(Item)

Item <- rep(13:18, each = 30, len = 90)
CK <- cbind(Item, CK)
rm(Item)

Item <- rep(19:24, each = 30, len = 90)
CL <- cbind(Item, CL)
rm(Item)

Item <- rep(25:30, each = 30, len = 90)
CM <- cbind(Item, CM)
rm(Item)

#Create one dataframe
result <- rbind(CI, CJ, CK, CL, CM)

rm(CI, CJ, CK, CL, CM)

#Creating coding for gender on assessed material
G <- rep( c( rep(c ( "M") , each = 15 ), rep( c( "F") , each = 15)) , 
          times = 5 )
result <- cbind(G, result)
rm(G)

#Creating coding for gender of assessor
A <- rep( c( rep(c ( "M") , each = 8 ), rep( c( "F") , each = 7)) , 
          times = 10 )
result <- cbind(A, result)
rm(A)

#Set factor for country
result$C = factor(result$C, 
                  levels = c("CI", "CJ", "CK", "CL", "CM"))

#Set factor for type of relationship
result$T = factor(result$T,
                  levels = c("TA", "TB", "TC"))

#Set factor for gender on assessed material
result$G <- factor(result$G,
                   levels = c("M", "F"))

#Set factor for gender of assessor
result$A <- factor(result$A,
                   levels = c("M", "F"))

contrasts(result$C)
contrasts(result$T)
contrasts(result$G)
contrasts(result$A)

#Set outcome as integer
result$Y <-as.integer(result$Y)

#Get a glimpse
glimpse(result)

table(result$Y, result$C, result$T)

############################ describe the data ################################
#Simple descriptive statistics
with(result, tapply(Y, list(C, T), median))
#     TA TB  TC
# CI 2.5  3 3.0
# CJ 3.5  4 4.0
# CK 5.0  5 4.5
# CL 5.5  5 5.5
# CM 6.5  6 6.5

with(result, tapply(Y, list(C, T), mean))
#          TA       TB       TC
# CI 2.633333 2.733333 2.733333
# CJ 3.666667 3.733333 3.766667
# CK 4.666667 4.833333 4.566667
# CL 5.300000 5.100000 5.133333
# CM 5.800000 5.766667 5.800000

with(result, tapply(Y, list(C, T), sd))
#          TA        TB        TC
# CI 1.496740 0.9444332 0.5832923
# CJ 1.397864 1.0148325 0.5683208
# CK 1.561019 1.0531835 0.6260623
# CL 1.556964 0.9948141 1.5698305
# CM 1.606023 1.0063020 1.6273524


#function to create a bar plot; created by CoPilot
bar_plot <- function(data){
  data %>%
    ggplot(aes(x = Y)) +
    geom_bar() +
    facet_grid(C ~ T)
}
bar_plot(result)




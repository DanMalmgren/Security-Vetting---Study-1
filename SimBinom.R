#Simulation of dichotomous data

library(rethinking)
library(reshape2)
library(tibble)
library(tidyr) 

#ensure repicable results
set.seed(5476)

#For Sweden
#TA = parental
TA <- rbinom(30, size = 1, prob = 0.08)

#TB = married
TB <- rbinom(30, size = 1, prob = 0.08)

#TC = born
TC <- rbinom(30, size = 1, prob = 0.08)

Sw <- data.frame(TA, TB, TC)
rm(TA, TB, TC)

#For UK
#TA = parental
TA <- rbinom(30, size = 1, prob = 0.08)

#TB = married
TB <- rbinom(30, size = 1, prob = 0.1)

#TC = born
TC <- rbinom(30, size = 1, prob = 0.12)

Uk <- data.frame(TA, TB, TC)
rm(TA, TB, TC)

#For Germany
#TA = parental
TA <- rbinom(30, size = 1, prob = 0.08)

#TB = married
TB <- rbinom(30, size = 1, prob = 0.1)

#TC = born
TC <- rbinom(30, size = 1, prob = 0.12)

Ge <- data.frame(TA, TB, TC)
rm(TA, TB, TC)

#For Brazil
#TA = parental
TA <- rbinom(30, size = 1, prob = 0.5)

#TB = married
TB <- rbinom(30, size = 1, prob = 0.6)

#TC = born
TC <- rbinom(30, size = 1, prob = 0.7)

Br <- data.frame(TA, TB, TC)
rm(TA, TB, TC)

#For Thailand
#TA = parental
TA <- rbinom(30, size = 1, prob = 0.6)

#TB = married
TB <- rbinom(30, size = 1, prob = 0.7)

#TC = born
TC <- rbinom(30, size = 1, prob = 0.8)

Th <- data.frame(TA, TB, TC)
rm(TA, TB, TC)

######################### End simulation of data ##############################

######################### Create the data frame ############################# 

#Convert each frame to long format
#CI = Thailand, CJ = Brazil, CK = Germany, CL = UK, CM = Sweden
CI <- gather(Th, key = "T", value = "D", c("TA", "TB", "TC" ))
CJ <- gather(Br, key = "T", value = "D", c("TA", "TB", "TC" ))
CK <- gather(Ge, key = "T", value = "D", c("TA", "TB", "TC" ))
CL <- gather(Uk, key = "T", value = "D", c("TA", "TB", "TC" ))
CM <- gather(Sw, key = "T", value = "D", c("TA", "TB", "TC" ))

#Remove the short frames from the environment
rm(Sw, Uk, Ge, Br, Th)

#Create a factor for country in each dataframe
#CI = Sweden, CJ = UK, CK = Germany, CL = Brazil, CM = Thailand 

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
Dresult <- rbind(CI, CJ, CK, CL, CM)

rm(CI, CJ, CK, CL, CM)

#Creating coding for gender on assessed material
G <- rep( c( rep(c ( "M") , each = 15 ), rep( c( "F") , each = 15)) , 
          times = 5 )
Dresult <- cbind(G, Dresult)
rm(G)

#Creating coding for gender of assessor
A <- rep( c( rep(c ( "M") , each = 8 ), rep( c( "F") , each = 7)) , 
          times = 10 )
Dresult <- cbind(A, Dresult)
rm(A)

#Set factor for country
Dresult$C = factor(Dresult$C, 
                  levels = c("CI", "CJ", "CK", "CL", "CM"))
levels(Dresult$C)

#Set factor for type of relationship
Dresult$T = factor(Dresult$T,
                  levels = c("TA", "TB", "TC"))

levels(Dresult$T)

#Detach rethinking
detach("package:rethinking", unload = TRUE)

#Set factor for country
Dresult$C <- factor(Dresult$C)

#Set factor for type of relationship
Dresult$T <- factor(Dresult$T)

# Apply sum contrasts (effects coding) for the C and T variables
contrasts(Dresult$C) <- contr.sum(n = length(levels(Dresult$C)))
contrasts(Dresult$T) <- contr.sum(n = length(levels(Dresult$T)))

contrasts(Dresult$C)
contrasts(Dresult$T)

#Get a glimpse
glimpse(Dresult)


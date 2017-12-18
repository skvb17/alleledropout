#getting database "genos"
getwd()
setwd("~/Desktop/Shannel-project")
load("~/Desktop/Shannel-project/Rori's Database.Rdata")


#getting allele frequency table
freq = read.table("AFs_AfA.dat.txt") #reads in data
t(freq)



#getting function that counts how many alleles match between two genotype arrays
alleleMatch = function(arr1, arr2) {
  matches = 0
  for(i in 1:13) {
    if((arr1[1,i] == arr2[1,i]) && (arr1[2,i] == arr2[2,i]) || (arr1[1,i] == arr2[2,i] && arr1[2,i] == arr2[1,i]))
      matches = matches + 2
    else if((arr1[1,i] == arr2[1,i]) || (arr1[1,i] == arr2[2,i]) || (arr1[2,i] == arr2[1,i]) || (arr1[2,i] == arr2[2,i]))
      matches = matches + 1
  }
  
  return(matches)
}



#getting function that counts the allelic matches between a profile and the first n profiles in a database
directMatch = function(array, database, nProfiles) {
  directmatches = rep(0,27)
  for (i in 1:nProfiles) {
    directmatches0 = alleleMatch(array,database[i,,])
    directmatches[directmatches0 + 1] = directmatches[directmatches0 + 1] + 1
  }
  
  return(directmatches)
}

#getting random profile for experiment and turning it into an array
dim(genos)  # 200K profiles 2 allelles  13 loci
dim(genos)[1] # 200k
ngenos = dim(genos)[1]
sample(ngenos,1) -> PROFILE  # choosing # of RANDOM profile(s)
PROFILE
genos[PROFILE,,] -> CONTROL # this require a manual insert of the choosen profile # 
CONTROL

#entering the control genotype into function directMatch
directMatch(CONTROL, genos, ngenos) 


# Creating a empty DataFrame name DROP13
DROP13 <- matrix(0, ncol = 27, nrow = 13)
DROP13 <- data.frame(DROP13)
rownames(DROP13) <- c("Drop1","Drop2", "Drop3", "Drop4","Drop5", "Drop6", "Drop7","Drop8", "Drop9", "Drop10","Drop11","Drop12", "Drop13")

# set Variables to NULL for the forLoops
result =NULL
directmatchResult=NULL

############# 
#####
# Result for ONE profile with 1 to 13 alleles drop

for (j in c(1:13)){ 
    dropout = function(DDprofile) {
                                    LOCIDROP <- sample(13,j,F)
                                    ALLEDROP <- sample(2,j,T) 
                                    for (i in c(1:j))
                                    if(DDprofile[ALLEDROP[i], LOCIDROP[i]] <- 0) {
                                                                                  print(DDprofile)
                                                                                  }
                                    return(DDprofile)
                                    }
                  result = dropout(CONTROL) 
                  print(result)
                  directmatchResult = unlist(directMatch(result, genos, ngenos)) -> DROP13[j,]
                  #print(unlist(directMatch(result, genos, ngenos)))  #remove # to have the 
                  }

# Creating a empty DataFrame name DROP_DIFF
DROP_DIFF <- matrix(0, ncol = 27, nrow = 13)
DROP_DIFF <- data.frame(DROP_DIFF)
rownames(DROP_DIFF) <- c("Drop1","Drop2", "Drop3", "Drop4","Drop5", "Drop6", "Drop7","Drop8", "Drop9", "Drop10","Drop11","Drop12", "Drop13")


############# 
#####
# Result for DIFFERENT profiles with 1 to 13 alleles drop

# Creating a empty DataFrame name DROP_DIFF
DROP_DIFF <- matrix(0, ncol = 27, nrow = 13)
DROP_DIFF <- data.frame(DROP_DIFF)
rownames(DROP_DIFF) <- c("Drop1","Drop2", "Drop3", "Drop4","Drop5", "Drop6", "Drop7","Drop8", "Drop9", "Drop10","Drop11","Drop12", "Drop13")


result = NULL
directmatchResult=NULL
  
  for (j in c(1:13)){ 
    
    sample(ngenos,1) -> PROFILE  # choosing # of RANDOM profile(s)
    PROFILE
    genos[PROFILE,,] -> CONTROL # this require a manual insert of the choosen profile # 
    CONTROL
    
    dropout = function(DDprofile) {
      LOCIDROP <- sample(13,j,F)
      ALLEDROP <- sample(2,j,T) 
      for (i in c(1:j))
        if(DDprofile[ALLEDROP[i], LOCIDROP[i]] <- 0) {
          print(DDprofile)
        }
      return(DDprofile)
    }
    result = dropout(CONTROL) 
    print(result)
    #print(unlist(directMatch(result, genos, ngenos)))
    directmatchResult = unlist(directMatch(result, genos, ngenos)) -> DROP_DIFF[j,]
    print(unlist(directMatch(result, genos, ngenos)))
  }




############# 
#####
# Running 50 DIFFERENT profiles with 13 alleleDROP 

# Creating a empty DataFrame name DROP_DIFF
SET50 <- matrix(0, ncol = 27, nrow = 50)
SET50 <- data.frame(SET50)


result =NULL
directmatchResult=NULL
DROP=NULL

for (k in c(1:50)){   # Chossing 1-50 profiles
  ngenos = dim(genos)[1]
  sample(ngenos,1) -> PROFILE  # choosing # of RANDOM profile(s)
  PROFILE
  genos[PROFILE,,] -> CONTROL # this require a manual insert of the choosen profile # 
  CONTROL
  
  for (j in c(13)){   # Choosing # of allele
    dropout = function(DDprofile) {
      LOCIDROP <- sample(13,j,F)
      ALLEDROP <- sample(2,j,T) 
      for (i in c(1:j))
        if(DDprofile[ALLEDROP[i], LOCIDROP[i]] <- 0) {
          print(DDprofile)
        }
      return(DDprofile)
    }
    result = dropout(CONTROL) 
    print(result)
    #print(unlist(directMatch(result, genos, ngenos)))
    directmatchResult = unlist(directMatch(result, genos, ngenos)) -> SET50[k,]
    print(unlist(directMatch(result, genos, ngenos)))
  }
}


df <- read.table(text="        X1   X2   X3    X4    X5    X6    X7    X8    X9   X10   X11   X12   X13  X14  X15 X16 X17 X18 X19 X20 X21 X22 X23 X24 X25 X26 X27
Drop1    0    4   52   280  1397  4668 11751 23069 33768 39255 35634 25657 14703 6525 2400 660 147  28   1   0   0   0   0   0   0   1   0
                 Drop2    0    3   43   248  1217  4221 11005 21794 33402 39223 36289 26763 15372 7020 2474 725 175  22   3   0   0   0   0   0   1   0   0
                 Drop3    0    5   55   377  1687  5652 13410 25289 35660 39699 34524 23374 12588 5381 1740 462  84  11   1   0   0   0   0   1   0   0   0
                 Drop4    0   10  103   629  2805  8631 18973 31888 40276 38913 29300 17250  7635 2637  759 158  30   2   0   0   0   0   1   0   0   0   0
                 Drop5    1   35  248  1227  4922 12911 25334 37014 41410 35233 23450 11828  4642 1345  342  50   7   0   0   0   0   1   0   0   0   0   0
                 Drop6    0   26  232  1365  5009 13906 26666 38362 41894 34706 22052 10522  3889 1083  239  45   3   0   0   0   1   0   0   0   0   0   0
                 Drop7    7   99  620  3091  9705 21711 35560 42481 38676 26694 13851  5476  1598  368   52   9   0   1   0   1   0   0   0   0   0   0   0
                 Drop8   13  276 1791  6955 18113 32620 42356 41367 30149 16590  7026  2112   519   98   13   1   0   0   1   0   0   0   0   0   0   0   0
                 Drop9   12  147 1107  5028 14552 28462 41039 42964 33527 20100  9100  2999   769  172   21   0   0   1   0   0   0   0   0   0   0   0   0
                 Drop10  56  678 3937 13160 28293 42103 44595 34983 20150  8736  2607   584    99   17    1   0   1   0   0   0   0   0   0   0   0   0   0
                 Drop11  70  773 4243 14252 30340 44035 45127 33677 17966  7015  2028   416    52    5    0   1   0   0   0   0   0   0   0   0   0   0   0
                 Drop12 247 2236 9628 24968 41390 47146 38711 22264  9500  3104   689   111     4    1    1   0   0   0   0   0   0   0   0   0   0   0   0
                 Drop13  93 1104 5633 16586 32535 44906 44113 31081 16184  5964  1534   246    19    2    0   0   0   0   0   0   0   0   0   0   0   0   0", header=T)

library(reshape2)
df <- melt(df)  #the function melt reshapes it from wide to long
df
df$rowid <- 1:13  #add a rowid identifying variable
head(df, 20)
ggplot(df, aes(variable, value, group=factor(rowid))) + geom_line(aes(color=factor(rowid)))


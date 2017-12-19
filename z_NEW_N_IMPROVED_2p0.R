#getting database "genos" by setting the directory and loading it into the environment
getwd()
setwd("/Users/shantothenel/Desktop/R")
load("/Users/shantothenel/Desktop/R/Rori's Database.Rdata")

#getting allele frequency table
freq = read.table("AFs_AfA.dat.txt") #reads in data
t(freq) #creates table to easily show allele frequencies



#PURPOSE: getting function that counts how many alleles match between two genotype arrays
  #INPUT: two arrays, where arr1 and arr2 are the profiles to be compared
  #OUTPUT: data that shows how many alleles matched in the same loci 
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




#PURPOSE getting function that counts the allelic matches between a profile and the first n profiles in a database
  #INPUT: an array (the control), the database that contains all the genotypes, the number of genotypes in the database
  #OUTPUT: a vector that contains how many profiles matched the array at 0 alleles to 26 alleles
directMatch = function(array, database, nProfiles) {
  directmatches = rep(0,27)
  for (i in 1:nProfiles) {
    directmatches0 = alleleMatch(array,database[i,,])
    directmatches[directmatches0 + 1] = directmatches[directmatches0 + 1] + 1
  }
  
  return(directmatches)
}

#getting random profile for experiment
dim(genos)  # 200K profiles, 2 allelles per loci w/ 13 loci in total
dim(genos)[1] # 200k
ngenos = dim(genos)[1]
sample(ngenos,1) -> PROFILE  # chose a single profile at random
PROFILE
genos[PROFILE,,] -> CONTROL # this requires a manual insert of the choosen profile # 
CONTROL

#entering the control genotype into function directMatch
directMatch(CONTROL, genos, ngenos) #works!




############# 
###########
#######
#####
##
#   Function for investigating ONE profile with a range of 1 to 13 alleles drop

# Creating a empty data frame named DROP13
DROP13 <- matrix(0, ncol = 27, nrow = 13)
DROP13 <- data.frame(DROP13)
rownames(DROP13) <- c("Drop1","Drop2", "Drop3", "Drop4","Drop5", "Drop6", "Drop7","Drop8", "Drop9", "Drop10","Drop11","Drop12", "Drop13")

# set variables to NULL for the for loops
result=NULL
directmatchResult=NULL

#creating a for loop that will dropout 1 alleles to 13 alleles for the same profile and printing directMatch into the empty dataframe
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
}


# Creating a Plot for the DROP13 data

write.table(DROP13, "~/Desktop/Shannel-project/DROP13.txt", sep="\t") # Exporting data into a table in .txt format
DROP13df = read.table("~/Desktop/Shannel-project/DROP13.txt") # Inserting the table and set it as variable

library(ggplot2)
library(reshape2)
DROP13df <- melt(DROP13df)  #the function melts and reshapes it from wide to long
DROP13df$rowid <- 1:13  #add a rowid identifying variable (# of rows in orginal dataset)
head(DROP13df, 30)
p <- ggplot(DROP13df, aes(variable, value, group=factor(rowid))) + geom_line(aes(color=factor(rowid)))

p + labs(colour = "# of Allele Dropout") + ylab("Profile Matches") + xlab("Alleles") + ggtitle("Dropout at One to Thirteen Random Allele(s) ( N = 1 )", subtitle = NULL)



############# 
###########
#######
#####
##
#   Function for investigating 50 DIFFERENT profiles with 13 alleleDROP 


# Creating a empty data frame named SET50
SET50 <- matrix(0, ncol = 27, nrow = 50)
SET50 <- data.frame(SET50)


result =NULL
directmatchResult=NULL
DROP=NULL

for (k in c(1:50)){   # Chossing 1-50 profiles
  ngenos = dim(genos)[1]
  sample(ngenos,1) -> PROFILE  # choosing # of RANDOM profile(s)
  PROFILE
  genos[PROFILE,,] -> CONTROL # this requires a manual insert of the choosen profile # 
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
    directmatchResult = unlist(directMatch(result, genos, ngenos)) -> SET50[k,]
    print(unlist(directMatch(result, genos, ngenos)))
  }
}

# Creating a plot for the SET50 data

write.table(SET50, "~/Desktop/Shannel-project/SET50.txt", sep="\t") # Exporting data into a table in txt format
SET50df = read.table("~/Desktop/Shannel-project/SET50.txt") # Inserting the table and set it as variable

library(ggplot2)
library(reshape2)
SET50df <- melt(SET50df)  #the function melts and reshapes it from wide to long
SET50df
SET50df$rowid <- 1:50  #add a rowid identifying variable ( # of rows in orginal dataset)
head(SET50df, 60)
p3 <- ggplot(SET50df, aes(variable, value, group=factor(rowid))) + geom_point(aes(color=factor(rowid)))
p3 + labs(colour = "Profiles") + ylab("Profile Matches") + xlab("Alleles") + ggtitle("Dropout at Thirteen Random Allele(s) ( N = 50 )", subtitle = NULL)

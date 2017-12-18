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




############# 
###########
#######
#####
##
#   Function for investigating ONE profile with a range of 1 to 13 alleles drop

# Creating a empty DataFrame name DROP13
DROP13 <- matrix(0, ncol = 27, nrow = 13)
DROP13 <- data.frame(DROP13)
rownames(DROP13) <- c("Drop1","Drop2", "Drop3", "Drop4","Drop5", "Drop6", "Drop7","Drop8", "Drop9", "Drop10","Drop11","Drop12", "Drop13")

# set Variables to NULL for the forLoops
result =NULL
directmatchResult=NULL

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
                  #print(unlist(directMatch(result, genos, ngenos)))  #remove # to have it print
                  }



# Creating a Plot for the DROP13 data

write.table(DROP13, "~/Desktop/Shannel-project/DROP13.txt", sep="\t") # Exporting data into a table in txt format
DROP13df = read.table("~/Desktop/Shannel-project/DROP13.txt") # Inserting the table and set it as variable

library(ggplot2)
library(reshape2)
DROP13df <- melt(DROP13df)  #the function melt reshapes it from wide to long
DROP13df$rowid <- 1:13  #add a rowid identifying variable ( # of rows in orginal dataset)
head(DROP13df, 30)
p <- ggplot(DROP13df, aes(variable, value, group=factor(rowid))) + geom_line(aes(color=factor(rowid)))

p + labs(colour = "# of Allele Dropout") + ylab("Profile Matches") + xlab("Alleles") + ggtitle("Dropout at One to Thirteen Random Allele(s) ( N = 1 )", subtitle = NULL)




############# 
###########
#######
#####
##
#   Function for investigating 13 DIFFERENT profiles with a range of 1 to 13 alleles drop


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

# Creating a Plot for the DROP_DIFF data

write.table(DROP13, "~/Desktop/Shannel-project/DROP_DIFF.txt", sep="\t") # Exporting data into a table in txt format
DROP_DIFFdf = read.table("~/Desktop/Shannel-project/DROP_DIFF.txt") # Inserting the table and set it as variable

library(ggplot2)
library(reshape2)
DROP_DIFFdf <- melt(DROP_DIFFdf)  #the function melt reshapes it from wide to long

DROP_DIFFdf$rowid <- 1:13  #add a rowid identifying variable ( # of rows in orginal dataset)
head(DROP_DIFFdf, 20)
p2 <- ggplot(DROP_DIFFdf, aes(variable, value, group=factor(rowid))) + geom_line(aes(color=factor(rowid)))

p2 + labs(colour = "# of Allele Dropout") + ylab("Profile Matches") + xlab("Alleles") + ggtitle("Dropout at One to Thirteen Random Allele(s) ( N = 13 )", subtitle = NULL)





############# 
###########
#######
#####
##
#   Function for investigating 50 DIFFERENT profiles with 13 alleleDROP 


# Creating a empty DataFrame name SET50
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

# Creating a Plot for the SET50 data

write.table(SET50, "~/Desktop/Shannel-project/SET50.txt", sep="\t") # Exporting data into a table in txt format
SET50df = read.table("~/Desktop/Shannel-project/SET50.txt") # Inserting the table and set it as variable

library(ggplot2)
library(reshape2)
SET50df <- melt(SET50df)  #the function melt reshapes it from wide to long
SET50df
SET50df$rowid <- 1:50  #add a rowid identifying variable ( # of rows in orginal dataset)
head(SET50df, 60)
p3 <- ggplot(SET50df, aes(variable, value, group=factor(rowid))) + geom_point(aes(color=factor(rowid)))
p3 + labs(colour = "Profiles") + ylab("Profile Matches") + xlab("Alleles") + ggtitle("Dropout at Thirteen Random Allele(s) ( N = 50 )", subtitle = NULL)




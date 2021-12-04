##=################################################=##
##                                                  ##
##       Examples in Individual Repeatability       ##
##       FNR647 - Quantitative Methods              ##
##              Spencer T. Gardner                  ##
##                  11/5/2021                       ##
##                                                  ##
##=################################################=##


# Intraclass Correlation Coefficient (Sokal and Rohlf, 1981)
# Ind. Repeatability (R) = proportion of observed variance that results from 
#    among-individual differences.

# Gaussian Data
# Difference between correlation-based, ANOVA-based, and LMM-based, and negative 
# repeatabilities

# Non-Gaussian Data


# Total Sum of Squares (TSS) = how much variation there is in the dependent variable
# Among-group Sum of Squares (SS_a) = variation in each mean to the grand mean
# Within-group Sum of Squares (SS_w) = variation due to individual differences

rm(list = ls(all.names = TRUE))              # Clear all objects from R

getwd()
#"C:/Users/Spencer Gardner/OneDrive/Documents/FNR647 - Quantitative methods for ecologists/QuantMethods_R_script/FNR647"


## === Install.packages ===
library(tidyverse)
library(ggplot2)
library(lme4)                 # linear mixed effects models
library(psych)                # describeBy table toolkit


## === Read data (.csv) check structure
trout <- read.csv("SeaTrout_migration_repeatability.csv", header = TRUE)
str(trout)     #401 obs. of 19 variables


#+############################+#
## ===== Summarize Data ======== 
#-############################-#

# isolate tag year
trout$tag.year <- format(as.Date(trout$Date_tagged, format = "%m/%d/%Y"), "%Y")

# isolate tag month and change to season
trout$tag.month <- format(as.Date(trout$Date_tagged, format = "%m/%d/%Y"), "%m")
unique(trout$tag.month)
     # reclassify as:
     #   Spring("03" "04" "05" "06")
     #   Fall("11" "12" "09" "10")
trout$tag.month <- recode(trout$tag.month, "03"="Spring","04"="Spring","05"="Spring",
                      "06"="Spring","09"="Fall","10"="Fall","11"="Fall","12"="Fall")

# for simplicity separate by season
Spring <- subset(trout, trout$tag.month == "Spring")
Fall <- trout[(trout$tag.month == "Fall"), ]

## === FOR SPRING-TAGGED FISH
# summarize data by season and life stage
sum.table1 <- describeBy(Spring$Length, list(Spring$tag.year,Spring$Lifestage_at_tagging), 
                         mat=TRUE, digits=2)

## === FOR FALL-TAGGED FISH
sum.table.f <- describeBy(Fall$Length, list(Fall$tag.year, Fall$Lifestage_at_tagging), 
                         mat = TRUE, digits = 2)



#+#########################+#
## ===== Format Date ========
#-#########################-#

# convert: chr(date) --> num(d.o.y)
# select columns to be formatted
columns <- c("Date_tagged","First.exit","First.entry","First.recap","Second.exit",
             "Second.entry","Second.recap", "Third.exit","Third.entry","Third.recap")

# convert date to day of year
for (i in columns)trout[,i] <- lubridate::yday(as.Date(trout[,i], format = '%m/%d/%Y'))

#check structure
str(trout)


#+#########################################+#
## ===== Visualize Repeated Behavior ===== ##
#-#########################################-#

# run below function to display multiple plots as one
RUN.FIRST <- {# source: Cookbook for R (http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/)
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
     library(grid)
     
     # Make a list from the ... arguments and plotlist
     plots <- c(list(...), plotlist)
     
     numPlots = length(plots)
     
     # If layout is NULL, then use 'cols' to determine layout
     if (is.null(layout)) {
          # Make the panel
          # ncol: Number of columns of plots
          # nrow: Number of rows needed, calculated from # of cols
          layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                           ncol = cols, nrow = ceiling(numPlots/cols))
     }
     
     if (numPlots==1) {
          print(plots[[1]])
          
     } else {
          # Set up the page
          grid.newpage()
          pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
          
          # Make each plot, in the correct location
          for (i in 1:numPlots) {
               # Get the i,j matrix positions of the regions that contain this subplot
               matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
               
               print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                               layout.pos.col = matchidx$col))
          }
     }
}}

# graphing commands
a <- {ggplot(data = trout, aes(x = First.entry, y = Second.entry))+
     geom_point()+
     geom_abline(slope = 1)+
     scale_y_continuous(limits = c(145,350))+
     scale_x_continuous(limits = c(145,350))+
     labs(x = "River entry season 1", y = "River entry season 2")+ 
     theme_bw()}
b <- {ggplot(data = trout, aes(x = Second.entry, y = Third.entry))+
     geom_point()+
     geom_abline(slope = 1)+
     scale_y_continuous(limits = c(145,350))+
     scale_x_continuous(limits = c(145,350))+
     labs(x = "River entry season 2", y = "River entry season 3")+ 
     theme_bw()}
c <- {ggplot(data = trout, aes(x = First.entry, y = Third.entry))+
     geom_point()+
     geom_abline(slope = 1)+
     scale_y_continuous(limits = c(145,350))+
     scale_x_continuous(limits = c(145,350))+
     labs(x = "River entry season 1", y = "River entry season 3")+ 
     theme_bw()}
multiplot(a,b,c, cols = 3)


#+#########################################+#
## ===== Calculate Repeated Behavior ===== ##
#-#########################################-#

# calculate group means
group.mean.one = mean(trout$First.entry, na.rm = TRUE)
group.mean.two = mean(trout$Second.entry, na.rm = TRUE)
group.mean.three = mean(trout$Third.entry, na.rm = TRUE)
     
# calculate grand means
grand.mean.1n2 = mean(c(trout$First.entry,trout$Second.entry), na.rm = TRUE)
grand.mean.2n3 = mean(c(trout$Second.entry,trout$Third.entry), na.rm = TRUE)
grand.mean.123 = mean(c(trout$First.entry,trout$Second.entry, trout$Third.entry), na.rm = TRUE)
     
# calculate Within-Group Sum of Squares
year.one.SS_w = ((trout$First.entry-group.mean.one)^2)
year.two.SS_w = ((trout$Second.entry-group.mean.two)^2)
year.three.SS_w = ((trout$Third.entry-group.mean.three)^2)

# calculate Between/Among - Group Sum of Squares
year.one.SS_a.1n2 = (((group.mean.one-grand.mean.1n2)^2)+((group.mean.two-grand.mean.1n2)^2))
year.two.SS_a.2n3 = (((group.mean.two-grand.mean.2n3)^2)+((group.mean.three-grand.mean.2n3)^2))
year.three.SS_a.123 = (((group.mean.one-grand.mean.123)^2)+((group.mean.two-grand.mean.123)^2)+
                            ((group.mean.three-grand.mean.123)^2))

k.1 = length(which(trout$First.entry >= 0))
k.2 = length(which(trout$Second.entry >= 0))
k.3 = length(which(trout$Third.entry >= 0))


trout$k.1n2.var <- (trout$First.entry-trout$Second.entry)
mean(trout$k.1n2.var, na.rm = TRUE)
test <- lm(k.1n2.var~id, data = trout)
anova(test)


## ==== TEST yr2 vs yr1 ==== ##
lm.1n2 <- lm(Second.entry~id, data = trout)
anova(lm.1n2)

#OR

test.1n2 = aov(Second.entry~id, data = trout)
summary(test.1n2)

test.1n2$terms

n = (1/(k-1))*(k.2-(2))
R = (21166-3876)/(21166 + (n) * 3876)
R


Repeatability (R) = (MSA - MSW)/(MSA + (n0 - 1) * MSW)
n0 = 

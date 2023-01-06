
#### R Basics ####

#Case sensative "A" vs "a" or "C" vs "c"
#Seperate commands on new line or ;
#Ctrl + 1 jump to source editor
#Ctrl + 2 jump to console
#Four hash marks to organize sections
#Run command with Run or Ctrl + Enter
#Set your working directory through session -> Set working directory -> Choose directory (ctrl+shift+H)
#You can also set it with setwd() and check using getwd()
getwd()

#swirl package for addtitional learning 

#Install package
install.packages("swirl")
#Load package into library
library(swirl)
#Run swirl
swirl()
#Check documentation
?swirl

#### ANOVA ####
install.packages("gge")
install.packages("agricolae")
install.packages("readxl") 
#or read.table() 
?read.table

library(gge)
library(agricolae)
library(readxl)

#import the experimental data
yielddata <- read_xlsx("yielddata.xlsx")
#I have decided to change the names of my columns
colnames(yielddata) <- c("gen", "env", "rep", "yield")
#check entries at top of list
head(yielddata)
#or look at a specific number of entries
head(yielddata, 12)
#check end of entries
tail(yielddata)

#ANOVA
?aov

anova <- aov(yield ~ gen + env:rep + env + env:gen, data = yielddata)

summary(anova)


View(anova)
#or call individual elements with
anova[1]

#Check coefficient of variance
cv<-cv.model(anova)

?LSD.test

meangenotypes<-LSD.test(anova, "gen", alpha = 0.05, group=TRUE)
meangenotypes$groups

meanenvironment<-LSD.test(anova, "env", group=TRUE)
meanenvironment$groups

#### Biplot ####
?gge
m <- gge(yielddata, yield ~ gen*env, scale=FALSE, gen.group= gen, env.group = env)
?biplot
biplot(m, main="GGE biplot", flip=FALSE, origin=0, hull=TRUE,
       col.gen= c(rep("deepskyblue4", 8), rep("darkred", 8)),
       pch.gen=rep(16, 16), col.env = ("Darkgreen"), cex.gen = 0.8,
       cex.env = .5)

####GRAPHICS####
library(datasets)
library(help = "datasets")
data(iris)
summary(iris)
head(iris)


#Graphics using GGPlot
install.packages("ggplot2")
#Load package into library
library(ggplot2)
library(help = "ggplot2")

#Visualize data

#scatterplot
ggplot(iris, aes(Sepal.Length, Sepal.Width, color=Species, shape=Species)) +
    geom_point() +
    xlab("Sepal Length") +  
    ylab("Sepal Width") +
    ggtitle("Sepal Length-Width") +
    theme(plot.title = element_text(hjust = 0.5))

#scatterplot facet wrap
ggplot(iris, aes(Sepal.Length, Sepal.Width, color=Species, shape=Species)) +
    geom_point() +
    xlab("Sepal Length") +  
    ylab("Sepal Width") +
    ggtitle("Sepal Length-Width") +
    theme(plot.title = element_text(hjust = 0.5))+
    facet_wrap(~ Species)

#box plot
ggplot(iris, aes(Species, Sepal.Width, fill=Species)) +
    geom_boxplot() +
    xlab("Species") +  
    ylab("Sepal Width") +
    ggtitle("Sepal Width") +
    theme_dark() +
    theme(plot.title = element_text(hjust = 0.5))

#histogram
ggplot(iris, aes(Sepal.Length, color = Species, fill = Species)) +
    geom_histogram(position = "identity", alpha = 0.5, binwidth = 0.2) +
    xlab("Sepal Length") +  
    ylab("Freq") +
    ggtitle("Sepal Length Histogram") +
    theme(plot.title = element_text(hjust = 0.5))

#density plot
ggplot(iris, aes(Sepal.Length, fill = Species)) +
    geom_density(alpha = 0.5) +
    xlab("Sepal Length") +  
    ylab("Freq") +
    ggtitle("Sepal Length Density") +
    theme(plot.title = element_text(hjust = 0.5))




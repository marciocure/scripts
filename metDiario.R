# 10 jun 2016
# marcio cure & marina hirota
# dados diarios da estacao de sao jose
dev.off()
library(SPEI)
library(zoo)
library(xts)
library(climdex.pcic)
citation("")
# working with station daily data************************************************
# loading data
?Climdex.pcic
?read.table
setwd("C:/Users/marcio/Documents/Mestrado/dados clima/dados inmet")

dat <- read.table("C:/Users/marcio/Documents/Mestrado/dados clima/dados inmet/estacao_sao_jose.txt",
                  sep = ";", header = T)
dat <- as.matrix(dat[, 1:10])
#head(dat)
# converting NAs into 0 (to sum) - not used
dat[is.na(dat)] <- 0

#dat <- na.omit(datt)
#na.action(dat)
#head(dat)


# testing the number of values
#teste <- vector("numeric", length = round(dim(dat)[1]/2))
#for(i in 1:dim(dat)[1]){
#  teste[i] <- as.numeric(dat[2*i, 3]) - as.numeric(dat[((2*i)-1), 3])
#}
#err <- which(teste ==0, arr.ind = T)
# summing rows 2 by 2 to eliminate double dates (different hours)
# when there are 2 rows for a the same date
# we will have daily data
uni <- unique(dat[,2])

datN <- matrix(NA, nrow = length(uni), ncol = dim(dat)[2]-3)
for(i in 1:length(uni)){
  aux <- which(dat[ ,2] == uni[i])
  
  if (length(aux) > 1){
    #print("hello 1")
    datN[i,] <- as.numeric(dat[aux[1], 4:10]) + as.numeric(dat[aux[2], 4:10])
  } else {
    #print("hello 2")
    datN[i,] <- as.numeric(dat[aux[1], 4:10])
  } 

}

# filling in the matrix with dates not present=========================
# dates
datesFull <- seq(as.Date("1jul1961", format = "%d%b%Y"), 
                 as.Date("30abr2016", format = "%d%b%Y"), by = "days")
datesUni <- as.Date(uni, format = "%d/%m/%Y")
# head(datesUni)
plot(datesFull,cex = 0.3, lwd = 2, main="Daily precipitation")
abline(lsfit(1:length(prec), prec), col = "red")
MannKendall(prec)
mtext("tau = 0.0405, p = 2.22e-16", side = 3, adj = 0.3, padj = 4)
MannKendall(prec)

# initializing full matrix
datNFull <- matrix(NA, nrow = length(datesFull), ncol = dim(datN)[2])
dev.off()
# finding the matching dates
for (i in 1:length(datesUni)){
  aux <- which(datesUni[i] == datesFull)
  datNFull[aux, ] <- datN[i, ]
}

# col names
colnames(datNFull) <- colnames(dat[ ,4:10])
#head(datNFull)
# calculating indices==============================================
# vars
prec <- zoo(datNFull[,1], datesFull)
temp <- zoo(datNFull[,6], datesFull)
tmin <- zoo(datNFull[,3], datesFull)
tmax <- zoo(datNFull[,2], datesFull)

# using climdex
#Sys.timezone()
auxDate <- seq(as.POSIXct("1jul1961", format = "%d%b%Y", tz = "America/Sao_Paulo"), 
               as.POSIXct("30abr2016", format = "%d%b%Y", tz = "America/Sao_Paulo"),
               by = "days")
climdex.dates <- as.PCICt(auxDate, cal = "gregorian",
                          tz = "BRT")
# data: tmax, tmin, prec, dates)
ci <- climdexInput.raw(datNFull[,2], datNFull[,3],
                       datNFull[,1], climdex.dates,
                       climdex.dates, climdex.dates,
                       base.range = c(1971,2000),
                       northern.hemisphere = FALSE)

length(climdex.dates)
cdd <- climdex.cdd(ci)

cwd <- climdex.cwd(ci)

precTot <- climdex.prcptot(ci)

r10 <- climdex.r10mm(ci)
MannKendall(r10)
r20 <- climdex.r20mm(ci)

r95Tot <-climdex.r95ptot(ci)

rx1M <- climdex.rx1day(ci, freq = "annual")

sdii <- climdex.sdii(ci)

su <- climdex.su(ci)

tr <- climdex.tr(ci)
plot(r10)
# playing around===================================================
?`climdex.pcic-package`
library(hydroTSM)

# Precipitation --------

var = prec

# plotting (general)
hydroplot(var, ptype="ts", pfreq="o", var.unit="mm", na.rm = T)

hydroplot(var, pfreq = "ma", var.type = "Precipitation", 
          ptype = "ts+boxplot+hist", var.unit = "mm", na.rm = T) #, from = "1961-jul-1", to = "1971-jul-1")
?hydroplot
# extracting the seasonal values for each year
DJF <- dm2seasonal(var, season="DJF", FUN=sum)
MAM <- dm2seasonal(var, season="MAM", FUN=sum)
JJA <- dm2seasonal(var, season="JJA", FUN=sum)
SON <- dm2seasonal(var, season="SON", FUN=sum)

boxplot(ts(DJF), ts(MAM), ts(JJA), ts(SON), main = "Precipitation (base range: 1971 - 2000)",
      ylab = "Precipitation (mm)" , col = "lightblue",names = c("summer", "autumn", "winter", "spring"))

hydroplot(var, pfreq="seasonal", FUN=mean, stype="default",
          season.names = c("summer", "autumn", "winter", "spring"))


# Temperature --------
# teste
var = temp

# plotting (general)
hydroplot(temp, ptype="ts", pfreq="o", var.unit="ºC", na.rm = T)

hydroplot(var, pfreq = "ma", var.type = "Temperature", 
          ptype = "ts+boxplot+hist", var.unit = "ºC", na.rm = T)
?hydroplot
# extracting the seasonal values for each year
DJF <- dm2seasonal(var, season="DJF", FUN=mean)
MAM <- dm2seasonal(var, season="MAM", FUN=mean)
JJA <- dm2seasonal(var, season="JJA", FUN=mean)
SON <- dm2seasonal(var, season="SON", FUN=mean)

boxplot(ts(DJF), ts(MAM), ts(JJA), ts(SON), main = "Temperature",
        names = c("summer", "autumn", "winter", "spring"))

hydroplot(var, pfreq="seasonal", FUN=mean, stype="default",
          season.names = c("summer", "autumn", "winter", "spring"))
daily2annual(prec, FUN=sum, na.rm = T)
daily2monthly()
??daily2monthly
library(zoo)
# daily zoo to monthly zoo
precM <- daily2monthly(prec, FUN=sum, na.rm=TRUE)
# creating a matrix with monthly values per year in each column
M <- matrix(precM, ncol = 12, byrow = F)
colnames(M) <- month.abb
rownames(M) <- x  #unique(format(time(x), trim = F ))
?as.Date
x <- c(1961:2015)
# plotting the monthly precipitation values
library(lattice)
print(matrixplot(M, ColorRamp="Precipitation",
   main="Monthly precipitation at São José Station, [mm/month]"))

## Para temperatura

# daily zoo to monthly zoo
tempM <- daily2monthly(temp, FUN=mean, na.rm=TRUE, out.fmt = "%Y")
# creating a matrix with monthly values per year in each column
Tm <- matrix(tempM, ncol = 12, byrow = F)
colnames(Tm) <- month.abb
rownames(Tm) <- x #unique(format(time(x), trim = F ))
?as.Date 
a <- c(1961:2016)
# plotting the monthly precipitation values
library(lattice)
print(matrixplot(Tm, ColorRamp="Temperature",
                 main="Monthly temperature at São José Station, [ºC/month]"))
??daily2annual

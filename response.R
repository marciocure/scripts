# 20 apr 2018
# response to reviewer 3, question 5, for the nerc/fapesp proposal
# written by rafa and lucy
# the idea is to use SPI values to support interannual variability
# and extremes over the last decades

library(zoo)
library(xts)
library(SCI)
library(hydroTSM)
library(Hmisc)
library(ncdf4)
library(cruts)
library(raster)
library(rgdal)

setwd("~/Documents/lucy/")

# opening cru netcdf dataset for reading
#ncCru <- nc_open("prec_1961_2015_SA.nc", readunlim = FALSE)

# getting prec values
#precCru <- ncvar_get(ncCru, varid = "pre")

cruRaster <- cruts2raster("prec_1961_2015_SA.nc",
                          timeRange = c("1961-01-01", "2015-12-16"))

# coordinates to extract time series
# veadeiros c(-47.619555, -14.065355)
# açu c(-36.946681, -5.579080)
lons <- c(-47.619555, -36.946681)
lats <- c(-14.065355, -5.579080)
coords <- cbind(lons, lats)
coordsProj <- SpatialPoints(coords, 
                            #proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
                            proj4string = CRS("+init=epsg:4326 +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

# time series => [1, ] = veadeiros; [2, ] = açu
tsP <- extract(cruRaster, coordsProj)

# apply spi
# scale to be calculated
tscale <- 6
months <- seq(as.Date("1jan1971", format = "%d%b%Y"), 
              as.Date("1dec2015", format = "%d%b%Y"), by = "months")

# VEADEIROS====================
# select period of 1971 - 2000 [121:480] as baseline
spi.ref.par <- fitSCI(as.vector(tsP[1,121:480]), first.mon = 1, distr = "gamma",
                      time.scale = tscale, p0 = TRUE) 
## apply the parameters of the reference period to the control plot (plot A)
spi.vead <- transformSCI(as.vector(tsP[1,121:660]), first.mon = 1, obj = spi.ref.par) 

# plotting
spi.vead.zoo <- zoo(spi.vead,months)
# finding extreme events
extremes <- which(coredata(spi.vead.zoo) <= -1.5, arr.ind = TRUE) # finding extremes
# accounting for separate (not followed extremes)
q <- which(c(0, diff(extremes)) != 1) 
extremes <- extremes[q]

plot(spi.vead.zoo, type = "b", xlab = "Time (months)",
     ylab = "SPI")
abline(h = 0)
abline(h = -1.5, col = "red", lwd = 3, lty = 2)
abline(h = 1.5, col = "darkblue", lwd = 3, lty = 2)
for (i in 1:length(extremes)){  # extreme dates
  abline(v = as.numeric(index(spi.vead.zoo[extremes[i]])), 
                        col = "gray", lty = "dotdash", lwd = 1)
  axis(3, at = as.numeric(index(spi.vead.zoo[extremes[i]])), 
       labels = index(spi.vead.zoo[extremes[i]]), las = 2,
       col.axis = "darkgreen", cex.axis = 0.7)
}


# ACU==========================
# select period of 1971 - 2000 [121:480] as baseline
spi.ref.par.acu <- fitSCI(as.vector(tsP[2,121:480]), first.mon = 1, distr = "gamma",
                      time.scale = tscale, p0 = TRUE) 
## apply the parameters of the reference period to the control plot (plot A)
spi.acu <- transformSCI(as.vector(tsP[2,]), first.mon = 1, obj = spi.ref.par.acu) 

# plotting
spi.acu.zoo <- zoo(spi.acu,months)
# finding extreme events
extremes <- which(coredata(spi.acu.zoo) <= -1.5, arr.ind = TRUE) # finding extremes
# accounting for separate (not followed extremes)
q <- which(c(0, diff(extremes)) != 1) 
extremes <- extremes[q]

plot(spi.acu.zoo, type = "b", xlab = "Time (months)",
     ylab = "SPI")
abline(h = 0)
abline(h = -1.5, col = "red", lwd = 3, lty = 2)
abline(h = 1.5, col = "darkblue", lwd = 3, lty = 2)
for (i in 1:length(extremes)){  # extreme dates
  abline(v = as.numeric(index(spi.acu.zoo[extremes[i]])), 
         col = "gray", lty = "dotdash", lwd = 1)
  axis(3, at = as.numeric(index(spi.acu.zoo[extremes[i]])), 
       labels = index(spi.acu.zoo[extremes[i]]), las = 2,
       col.axis = "darkgreen", cex.axis = 0.7)
}

# boxplots to evaluate interannual variability===============
# general
boxplot(tsP[1,121:660])
boxplot(tsP[2,121:660])

# boxplot per month
# VEADEIROS
cmonth <- format(time(zoo(tsP[1, 121:660],months)), "%b")
# Creating ordered monthly factors
mon <- factor(cmonth, levels=unique(cmonth), ordered=TRUE)
#
boxplot(tsP[1, 121:660] ~ mon, ylab = "Monthly precipitation (mm)",
        xlab = "Time (months)")
abline(h = 100, col = "red", lwd = 3, lty = 2)

# ACU
boxplot(tsP[2, 121:660] ~ mon, ylab = "Monthly precipitation (mm)",
        xlab = "Time (months)")
abline(h = 100, col = "red", lwd = 3, lty = 2)

# interannual variation of msi an total annual precipitation
source('~/Documents/scripts R/climatic indices/msi2D.R')
source('~/Documents/scripts R/climatic indices/ap2D.R')

# MSI=========
yr <- seq(as.Date("1jan1971", format = "%d%b%Y"), 
          as.Date("1dec2015", format = "%d%b%Y"), by = "year")
msi2 <- msi2D(tsP[, 121:660])
plot(zoo(msi2[1,], yr), xlab = "Time (years)", 
     ylab = "Markham seasonality index (MSI)",
     ylim = c(0,1), col = "blue", lwd = 2, type = "b")
abline(h = max(msi2[1,]), lty = "dotdash", col = "lightblue")
abline(h = min(msi2[1,]), lty = "dotdash", col = "lightblue")

lines(zoo(msi2[2,], yr), col = "red", lwd = 2, type = "b")
abline(h = max(msi2[2,]), lty = "dotdash", col = "tomato")
abline(h = min(msi2[2,]), lty = "dotdash", col = "tomato")

legend("bottomleft", legend = c("Veadeiros", "Açu"),
       col = c("blue", "red"), lwd = c(2, 2), 
       lty = c(1, 1))

# ANNUAL PRECIP=========
ap2 <- ap2D(tsP[, 121:660])
plot(zoo(ap2[1,], yr), xlab = "Time (years)", 
     ylab = "Annual precipitation total (mm)",
     ylim = c(0, max(ap2)), 
     col = "blue", lwd = 2, type = "b")
abline(h = max(ap2[1,]), lty = "dotdash", col = "lightblue")
abline(h = min(ap2[1,]), lty = "dotdash", col = "lightblue")

lines(zoo(ap2[2,], yr), col = "red", lwd = 2, type = "b")
abline(h = max(ap2[2,]), lty = "dotdash", col = "tomato")
abline(h = min(ap2[2,]), lty = "dotdash", col = "tomato")

legend("bottomleft", legend = c("Veadeiros", "Açu"),
       col = c("blue", "red"), lwd = c(2, 2), 
       lty = c(1, 1))

# coefficient of variation (annual and monthly)
# annual
cv <- apply(ap2, 1, function(x) {sd(x)/mean(x)})
#> cv
#[1] 0.1477075 0.3963311

# monthly
dfP1 <- as.data.frame(matrix(tsP[1,121:660], 
                             nrow = length(tsP[1,121:660])/12, 
                             ncol = 12, byrow = TRUE))
colnames (dfP1) <- c("Jan", "Feb", "Mar", "Apr",
                     "May", "Jun", "Jul", "Aug",
                     "Sep", "Oct", "Nov", "Dec")
dfP2 <- as.data.frame(matrix(tsP[2,121:660], 
                             nrow = length(tsP[2,121:660])/12, 
                             ncol = 12, byrow = TRUE))
colnames (dfP2) <- c("Jan", "Feb", "Mar", "Apr",
                     "May", "Jun", "Jul", "Aug",
                     "Sep", "Oct", "Nov", "Dec")

cv1 <- apply(dfP1, 2, function(x) {sd(x)/mean(x)}) # vead
#Jan       Feb       Mar       Apr       May       Jun       Jul       Aug       Sep       Oct       Nov       Dec 
#0.4550251 0.4115408 0.4044015 0.4415406 0.7148514 1.4937996 2.1207726 1.2899709 0.7181939 0.4707078 0.2849866 0.3525582 
cv2 <- apply(dfP2, 2, function(x) {sd(x)/mean(x)}) # acu
#Jan       Feb       Mar       Apr       May       Jun       Jul       Aug       Sep       Oct       Nov       Dec 
#0.8372159 0.5599876 0.5381032 0.6576584 0.6249922 0.6990305 0.8998591 0.9749686 1.3676851 1.7566106 0.9778074 0.9672174 
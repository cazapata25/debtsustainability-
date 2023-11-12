
## ----------------------------------------------------------------
## Paper
## Fiscal regimes and debt sustainability in Colombia
## ----------------------------------------------------------------

# Ligraries
library(TSA)
library(mFilter)
library(seasonal)
library(car)
library(zoo)
library(tseries)
library(readr)
library(orcutt)
library(ggplot2)
library(MSwM)
library(sandwich)
library(strucchange)

#rm(list=ls())
library(readr)
#data80 <- read_csv("data80.csv")
summary(data80)

#--------------------------------------------------------------------------

# Series de tiempo de datos anuales
datos <- data80
len <- nrow(datos)
dt <- ts(datos$d, start = c(1980), freq = 1)    # Debt
st <- ts(datos$s, start = c(1980), freq = 1)    # Balance total
bpt <- ts(datos$bp, start = c(1980), freq = 1)  # Primary balance 
gt <- ts(datos$g, start = c(1980), freq = 1)    # Expendiature - without interest
it <-  ts(datos$i, start = c(1980), freq = 1)   # Interest (payments)
rt <-  ts(datos$r, start = c(1980), freq = 1)   # Real interest
gdpn <- ts(log(datos$gdpn), start = c(1980), freq = 1) # Log GDP
expend <- ts(log(datos$expend), start = c(1980), freq = 1) # Log expenditures
yt <-  ts(datos$y, start = c(1980), freq = 1)         # Growth rate GDP
radjt <-  ts(datos$radj, start = c(1980), freq = 1)   # Interest adjusted
pi <- ts(datos$pi, start = c(1980), freq = 1)         # Inflation
ti <-  ts(datos$ti, start = c(1980), freq = 1)        # Terms of trade index 
oil <-  ts(datos$oil, start = c(1980), freq = 1)      # Growth oil prices 
years <- 1980:2021

## ----------------------------------------------------------------

# Figure 1. Public debt and fiscal deficit in Colombia

df1 <- data.frame(dt,years)
df2 <- data.frame(bpt,years)
ylim.prim <- c(0, 0.6) 
ylim.sec <- c(-0.06,0.02)  
b1 <- diff(ylim.prim)/diff(ylim.sec)
a1 <- ylim.prim[1] - b1*ylim.sec[1] 

windows()
ggplot() + 
    geom_area(aes(x=years, y=dt),fill = "grey",alpha=0.8) + 
    geom_line(aes(x=years, y=a1+bpt*b1)) +
    scale_y_continuous(name = "Public debt", 
                       sec.axis = sec_axis(~(.-a1)/b1, name="Fiscal deficit")) +
    ggtitle("Fig 1. Output gap and cyclical component of spending ") +
    labs(x = "Years") 


## ----------------------------------------------------------------

# HP Filter

gdp.hp <- hpfilter(gdpn, freq = 100)
gdp.gap <- gdp.hp$cycle

windows()
par(mfrow = c(1, 2), cex = 0.8)
plot.ts(gdpn, ylab = "Log GDP")  # plot time series
lines(gdp.hp$trend, col = "red")  # include HP trend
legend("topleft", legend = c("Log GDP", "Trend"), lty = 1, 
       col = c("black", "red"), bty = "n")
plot.ts(gdp.hp$cycle, ylab = "", main="Output Gap")  # plot cycle
legend("topright", legend = c("Cycle"), lty = 1, col = c("black"), 
       bty = "n")

g.hp <- hpfilter(expend, freq = 100)
g.cycle <- g.hp$cycle

windows()
par(mfrow = c(1, 2), cex = 0.8)
plot.ts(expend, ylab = "Log Gasto")  # plot time series
lines(g.hp$trend, col = "red")  # include HP trend
legend("topleft", legend = c("Log Expend.", "Trend"), lty = 1, 
       col = c("black", "red"), bty = "n")
plot.ts(g.hp$cycle, ylab = "",main="cyclical component of government expenditure ")  # plot cycle
legend("topright", legend = c("Cycle"), lty = 1, col = c("black"), 
       bty = "n")


df3 <- data.frame(gdp.gap,years)
df4 <- data.frame(g.cycle,years)
ylim.prim2 <- c(-0.2, 0.2)   
ylim.sec2 <- c(-0.2,0.2)  
d <- diff(ylim.prim2)/diff(ylim.sec2)
c <- ylim.prim2[1] - d*ylim.sec2[1] 


# Figure 2. Public debt and fiscal deficit in Colombia

windows()
ggplot(df3) + 
    geom_line(aes(x=years, y=gdp.gap)) + 
    geom_line(aes(x=years, y=c+g.cycle*d),linetype=2) +
    scale_y_continuous(name = "Output Gap", 
                       sec.axis = sec_axis(~(.-c)/d, name="cyclical public spending")) +
    ggtitle("Fig 2. HP filtered output gap and cyclical public spending") +
    theme(legend.position="left") + labs(x = "Years") 

# -----------------------------------------------------------------

# Adjustment for cyclical components 

gdp.gapt <- ts(gdp.gap, start = c(1980), freq = 1)
g.cyclet <- ts(g.cycle, start = c(1980), freq = 1)
modbp_adj <- lm(bpt ~ gdp.gapt + g.cyclet)
summary(modbp_adj)

(alpha <- modbp_adj[["coefficients"]][1])
(alphagap <- modbp_adj[["coefficients"]][2])
(alphaexp <- modbp_adj[["coefficients"]][3])

(bp.adj_cycle <- bpt -alpha + alphagap*gdp.gapt + alphaexp*g.cyclet)

# Figure 3. Relationship between fiscal balance and debt 

windows()
plot(dt[1:41],bpt[2:42], ylab ="Primary balance /GDP", xlab="Lag Debt/GDP",
     pch=20,col="darkblue") 
abline(lm(bpt[2:42]~dt[1:41]),col="darkblue")


windows()
plot(dt[1:41],bp.adj_cycle[2:42],ylab ="Primary balance adjusted/GDP", xlab="Lag Debt/GDP",
     pch=20,col="darkblue") #
abline(lm(bp.adj_cycle[2:42]~dt[1:41]),col="darkblue")

## ----------------------------------------------------------------
## ----------------------------------------------------------------

## Model results
# Bohn's FRF
# Var. control: gdp.gap, g.cycle

bptt <- bpt[2:len]
bpts <- bp.adj_cycle[2:len]
lagdt <- dt[1:len-1]
gdp.gapt <- gdp.gap[2:len]
g.cyclet <- g.cycle[2:len]

# 1. Model FRF lineal
model1 <- lm(bpts ~ lagdt + gdp.gapt + g.cyclet) # Ec. (X)
summary(model1)

NLLSco1 = cochrane.orcutt(model1)
summary(NLLSco1)

#-----------------------------------------------------------
# Model FRF lineal 
# Var. control: pib.gap, g.gap, pi, ti, oil

pit <- pi[2:len]
tit <- ti[2:len]
oilt <- oil[2:len]

model2 <-lm(bpts ~ lagdt+gdp.gapt+g.cyclet+pit+tit+oilt)
summary(model2)

NLLSco2 = cochrane.orcutt(model2)
summary(NLLSco2)

#-----------------------------------------------------------

# Model FRF nonlineal 
# Var. de control: gdp.gap, g.cycle
lagdt2 <- lagdt^2
lagdt3 <- lagdt^3

model3 <- lm(bpts ~ lagdt+lagdt2+lagdt3+gdp.gapt+g.cyclet)
summary(model3)

NLLSco3 = cochrane.orcutt(model3)
summary(NLLSco3)

# Var. control: pib.gap, g.gap, pi, ti, oil

model4 <- lm(bpts ~ lagdt+lagdt2+lagdt3+gdp.gapt+g.cyclet+pit+tit+oilt)
summary(model4)

NLLSco4 = cochrane.orcutt(model4)
summary(NLLSco4)

## ----------------------------------------------------------------
## ----------------------------------------------------------------

## Markov-Switching for model 1

model1.MS <- msmFit(model1, k=2, sw=c(T,T,T,T,T), p=0)
summary(model1.MS)

windows()
plotProb(model1.MS,which=2)


#slotNames(mod@Fit)
filtprob=(model1.MS@Fit@filtProb)  # Filtered probability
filtprob.res = ts(filtprob[,2],start=c(1980),end=c(2021),freq=1)
smoprob=(model1.MS@Fit@smoProb) # Smoothed probability
smoprob.res = ts(smoprob[,2],start=c(1980),end=c(2021),freq=1)


# Figure 4. Probability and filtered probability for unsustainable regime

windows()
plot(filtprob.res,type="h",lwd =10,col="grey", ylab= "Probability",
     main="Fig. 4. Probability and filtered probability for unsustainable regime",cex=0.7)
lines(smoprob.res,col="black")
legend("top",legend = c("Smoothed probability","Filtered probability"),
       lty=c(1,1),cex=0.8,col=c("gray","black"))

windows()
ggplot() + 
    geom_area(aes(x=years, y= filtprob.res),fill = "grey") + 
    geom_line(aes(x=years, y= smoprob.res)) +
    ggtitle("Probability") +
    labs(x = "years", y="Probability") + 
    ggtitle("Fig 4. Probability and filtered probability for unsustainable regime")

## ------------------------------------------------------------
## ------------------------------------------------------------

# Debt-Stabilizing condition and Structural breaks
bp.ri <- breakpoints(radjt ~ 1, h = 2)
summary(bp.ri)
bp.ri2 <- lm(radjt ~ breakfactor(bp.ri, breaks = 2))
summary(bp.ri2)

bp.rits <- ts(fitted(bp.ri2, start = 1980))

# Figure 5. Growth-adjusted real interest rate
windows()
ggplot() + 
    geom_area(aes(x=years, y=radjt),fill = "grey") + 
    geom_line(aes(x=years, y=bp.rits)) + labs(x = "Years") +
    ggtitle("Growth-adjusted real interest rate")


# ------------------------------------------------------------------------
# ------------------------------------------------------------------------

## Robustness Checks

## Markov-Switching for model 2

model2.MS <- msmFit(model2, k=2, sw=c(T,T,T,T,T,T,T,T), p=0)
summary(model2.MS)

# --------------------------------------------------------------------
# datosf <- cbind(bpts,lagdt,gdp.gapt,g.cyclet,pit,tit,oilt,radjt[2:42],bp.rits[2:42],filtprob[,1],smoprob[,1])
# 
# library(writexl)
# datosdf <- data.frame(datosf)
# write_xlsx(datosdf,"C:/Users/data.xlsx")



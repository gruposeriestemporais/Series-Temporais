### Valor em Risco

rm(list=ls())

### Obtendo os dados

library(quantmod)
library(ggplot2)

TS<- getSymbols("^SSMI", src = "yahoo", from = "1995-11-20", to = "2019-04-17", auto.assign = FALSE)

# Como a série possui NAs podemos utilizar o pacote imputeTS.

library(imputeTS)
TS<-na.kalman(TS)
TS

head(TS)
tail(TS)
summary(TS)
str(TS)

### Gráfico da série

ggplot(TS, aes(x = index(TS), y = TS[,6])) + geom_line(color = "darkblue") +
  ggtitle("Série de preços do SMI") +
  xlab("Data") + ylab("Preço ($)") + theme(plot.title = element_text(hjust = 0.5)) + 
  scale_x_date(date_labels = "%b %y", date_breaks = "6 months")

### Retornos

ret <- diff(log(TS$SSMI.Adjusted))
ret <- ret[-1,]

summary(ret) # resumo
library(fBasics)
basicStats(ret)

sd(ret) # desvio padrão
mean(ret) # média

### Gráfico série de retornos

ggplot(ret, aes(x = index(ret), y = ret)) + geom_line(color = "deepskyblue4") +
  ggtitle("Série de retornos") + xlab("Data") + ylab("Retorno") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_date(date_labels = "%b %y", date_breaks = "6 months")

### Histograma

hist(ret,nclass=100,probability=T,col="lightblue3",main="",xlab= "Retorno SMI")
LX <- seq(from=min(ret),to=max(ret),length=1000)
lines(LX,dnorm(LX,mean=mean(ret),sd=sd(ret)))


### VaR usando o modelo GARCH

length(ret)
library(rugarch)
library(zoo)

retdf <- as.data.frame(ret)
hv <- rollapply(retdf, 150, sd)
l <- length(hv)
hv <- hv[(l-2000+1):l] # serão utilizadas para validação do modelo

gspec11 <- ugarchspec(variance.model = list(model = "sGARCH", 
                     garchOrder = c(1, 1)),
                      mean.model=list(armaOrder=c(0,0), 
                      include.mean = FALSE), distribution="norm")

roll11 <- ugarchroll(gspec11,retdf, n.start=3990,
                     refit.every = 25, refit.window = "moving",
                     VaR.alpha = c(0.025, 0.05))

rt <- ret[3990:5950]  # retornos para validação

# VaR 5%

VaR <- sd(rt) * qnorm(0.05)
VaRMA <- hv * qnorm(0.05) # Média móvel
VaRGARCH <- roll11@forecast$VaR[,2]  

xaxis <- rownames(ret)
xaxis <- xaxis[3990:5989]
xaxis <- c(xaxis[1], xaxis[500], xaxis[1000], xaxis[1500], xaxis[2000])


# VaR Plot

plot(rt, type = "l", pch = 17, cex = 0.8,  col = gray(0.4, 0.7),
     ylab = "Retornos", main = "95%  Previsão VaR", xaxt = "n")
axis(1, at=c(1, 500, 1000, 1500,2000), labels=xaxis)
lines(VaRGARCH, col = 1)
lines(VaRMA, col = 4)
abline(h=VaR, col = 2)
legend('left', c("GARCH(1,1)", "MA 150 days", "Static VaR") , 
       lty=1, col=c(1,4,2), bty='n', cex=0.7)

# VaR volatilidade

plot(abs(rt), type = "l", col = grey(0.4, 0.5), 
     ylab = "Retornos absolutos", main = "Previsão Volatilidade VaR", xaxt = "n")
axis(1, at=c(1, 500, 1000, 1500, 2000), labels=xaxis)
Sigma11 <- roll11@forecast$density$Sigma
lines(Sigma11, col = 1)
lines(hv, col = 4)
abline(h=sd(rt), col = 2)
legend('topleft', c("GARCH(1,1)", "MA 150 days", "Unconditional") , 
       lty=1, col=c(1,4, 2), bty='n', cex=.75)

##### VaR para ARMA(p,q)-GARCH(m,n)

rm(list=ls())

# Simulando um ARMA(1,1)-GARCH(1,1) com inovação t

nu <- 3
fixed.p <- list(mu = 0, 
          ar1 = 0.5, 
          ma1 = 0.3, 
          omega = 4, 
          alpha1 = 0.4, 
          beta1 = 0.2, 
          shape = nu) 
          armaOrder <- c(1,1)
          garchOrder <- c(1,1) 
varModel <- list(model = "sGARCH", garchOrder = garchOrder)
spec <- ugarchspec(varModel, mean.model = list(armaOrder = armaOrder),
              fixed.pars = fixed.p, distribution.model = "std")


n <- 1000 
x <- ugarchpath(spec, n.sim = n, m.sim = 1, rseed = 271) 
X <- fitted(x) 
sig <- sigma(x)
eps <- x@path$residSim
stopifnot(all.equal(X,   x@path$seriesSim, check.attributes = FALSE),
          all.equal(sig, x@path$sigmaSim,  check.attributes = FALSE))

plot(X,   type = "l", xlab = "t", ylab = expression(X[t]))
plot(sig, type = "h", xlab = "t", ylab = expression(sigma[t]))

# Modelo ARMA - GARCH

spec <- ugarchspec(varModel, mean.model = list(armaOrder = armaOrder),
                   distribution.model = "std") 
fit <- ugarchfit(spec, data = X) 

mu. <- fitted(fit) 
sig. <- sigma(fit)

stopifnot(all.equal(as.numeric(mu.),  fit@fit$fitted.values),
          all.equal(as.numeric(sig.), fit@fit$sigma))

plot(X, type = "l", xlab = "t",
     ylab = expression("Data"~X[t]~"and fitted values"~hat(mu)[t]))
lines(as.numeric(mu.), col = adjustcolor("blue", alpha.f = 0.5))
legend("bottomright", bty = "n", lty = c(1,1),
       col = c("black", adjustcolor("blue", alpha.f = 0.5)),
       legend = c(expression(X[t]), expression(hat(mu)[t])))

### Cálculo do VaR

alpha <- 0.95

VaR. <- as.numeric(quantile(fit, probs = alpha))

nu. <- fit@fit$coef["shape"] 
VaR.. <- as.numeric(mu. + sig. * sqrt((nu.-2)/nu.) * qt(alpha, df = nu.))
stopifnot(all.equal(VaR.., VaR.))


fspec <- getspec(fit) 
setfixed(fspec) <- as.list(coef(fit)) 
m <- ceiling(n / 10) 
pred <- ugarchforecast(fspec, data = X, n.ahead = 1, n.roll = m-1, out.sample = m)

mu.predict <- fitted(pred) 
sig.predict <- sigma(pred) 
VaR.predict <- as.numeric(quantile(pred, probs = alpha)) 


stopifnot(all.equal(mu.predict, pred@forecast$seriesFor, check.attributes = FALSE),
          all.equal(sig.predict, pred@forecast$sigmaFor, check.attributes = FALSE))

VaR.predict. <- as.numeric(mu.predict + sig.predict * sqrt((nu.-2)/nu.) *
                             qt(alpha, df = nu.)) 
stopifnot(all.equal(VaR.predict., VaR.predict))

B <- 1000
set.seed(271)
X.sim.obj <- ugarchpath(fspec, n.sim = m, m.sim = B)

X.sim <- fitted(X.sim.obj) 
sig.sim <- sigma(X.sim.obj)
eps.sim <- X.sim.obj@path$residSim 
VaR.sim <- (X.sim - eps.sim) + sig.sim * sqrt((nu.-2)/nu.) * qt(alpha, df = nu.)
VaR.CI <- apply(VaR.sim, 1, function(x) quantile(x, probs = c(0.025, 0.975)))

### Gráfico

yran <- range(X, 
              mu., VaR., 
              mu.predict, VaR.predict, VaR.CI) 
myran <- max(abs(yran))
yran <- c(-myran, myran) 
xran <- c(1, length(X) + m) 

btest <- VaRTest(1-alpha, actual = -X,
                 VaR = quantile(ugarchfit(spec, data = -X), probs = 1-alpha))
btest$expected.exceed 

plot(X, type = "l", xlim = xran, ylim = yran, xlab = "Time t", ylab = "",
     main = "Simulated ARMA-GARCH, fit, VaR, VaR predictions and CIs")
lines(as.numeric(mu.), col = adjustcolor("darkblue", alpha.f = 0.5))
lines(VaR., col = "darkred")
mtext(paste0("Expected exceed.: ",btest$expected.exceed,"   ",
             "Actual exceed.: ",btest$actual.exceed,"   ",
             "Test: ", btest$cc.Decision),
      side = 4, adj = 0, line = 0.5, cex = 0.9) 

## Previsão

t. <- length(X) + seq_len(m) 
lines(t., mu.predict, col = "blue") 
lines(t., VaR.predict, col = "red") 
lines(t., VaR.CI[1,], col = "orange")
lines(t., VaR.CI[2,], col = "orange") 
legend("bottomright", bty = "n", lty = rep(1, 6), lwd = 1.6,
       col = c("black", adjustcolor("darkblue", alpha.f = 0.5), "blue",
               "darkred", "red", "orange"),
       legend = c(expression(X[t]), expression(hat(mu)[t]),
                  expression("Predicted"~mu[t]~"(or"~X[t]*")"),
                  substitute(widehat(VaR)[a], list(a = alpha)),
                  substitute("Predicted"~VaR[a], list(a = alpha)),
                  substitute("95%-CI for"~VaR[a], list(a = alpha))))
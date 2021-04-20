# ZHENG XIN, r0766879, KU LEUVEN
# R version: R version 3.5.0 (2018-04-23) -- "Joy in Playing"
# R packages and Data preparation ====
library(lmtest)
library(MASS)
library(gvlma)
library(rstatix)
library(psych)
library(DescTools)
library(performance)
library(car)
library(robustbase)
library(caret)
library(TeachingDemos)
library(segmented)
library(nortest)
# Data preparation
rm(list = ls())
data.full <- read.table('invertebrate.txt', header = T)
set.seed(0766879)
d.test <- sample(1:dim(data.full)[1], 200 )
data.test <- data.full[d.test, ]
data.training <- data.full[-d.test, ]
# Q1====
# Perform an exploratory analysis of the variables 
# (compute descriptive statistics and make histograms, boxplots, scatter plots, . . . )
attach(data.training)
# Descriptive statistics
str(data.training)
summary(data.training)
# Correlation Matrix with P-values
cor_mat(data.training)
cor_pmat(data.training)
cor <- cor(data.training[, !names(data.training) == 'SWI']) # correlation between predictor variables
cor # High correlation bwtween duration and temperature (Question 6)
dim(data.training)

# Exploratory analysis
histNorm <- function(x, densCol = "darkblue", xlab = ''){
  m <- mean(x)
  std <- sqrt(var(x))
  h <- max(hist(x,plot=FALSE)$density)
  d <- dnorm(x, mean=m, sd=std)
  maxY <- max(h,d)
  hist(x, prob=TRUE,
       xlab = xlab, ylab="Frequency", ylim=c(0, maxY),
       main="Histogram")
  curve(dnorm(x, mean=m, sd=std),
        col=densCol, lwd=2, add=TRUE)
}

par(mfrow = c(3,2))
histNorm(data.training$SWI, xlab = "SWI")
histNorm(data.training$SWF, xlab = "SWF")
histNorm(data.training$temperature, xlab = "temperature")
histNorm(data.training$size, xlab = "size")
histNorm(data.training$management, xlab = "management")
# management as a Categorical predictor, not normally distributed
histNorm(data.training$duration, xlab = "duration")

# boxplots
par(mfrow = c(3,2))
boxplot(SWI, main = "Boxplot of SWI") # two outliers, both smaller than 4.5.
boxplot(SWF, main = "Boxplot of SWF") # three outliers
boxplot(temperature, main = "Boxplot of temperature") # three outliers
boxplot(size, main = "Boxplot of size")
boxplot(management, main = "Boxplot of management")
boxplot(duration, main = "Boxplot of duration") # one outlier
lab_y <- seq(1,200)
source("https://raw.githubusercontent.com/talgalili/R-code-snippets/master/boxplot.with.outlier.label.r")
# Load the function to label all the outliers in a boxplot
par(mfrow = c(2,3))
# OBS 51, 286
boxplot.with.outlier.label(SWI, row.names(data.training), main = "Boxplot of SWI")
# OBS 51, 139, 351
boxplot.with.outlier.label(SWF, row.names(data.training), main = "Boxplot of SWF")
# OBS 1, 3, 397
boxplot.with.outlier.label(temperature, row.names(data.training), main = "Boxplot of temperature")
boxplot.with.outlier.label(size, row.names(data.training), main = "Boxplot of size")
boxplot.with.outlier.label(management, row.names(data.training), main = "Boxplot of management")
# OBS 7
boxplot.with.outlier.label(duration, row.names(data.training), main = "Boxplot of duration")
# ~ 1, 3, 7, 51, 139, 286, 351, 397

# Scatter plot
par(mfrow = c(1,1))
pairs(data.training)
pairs(data.training, panel = function(x,y) {points(x,y); lines(lowess(x,y), col = "red")})
pairs.panels(data.training) # Robust fitting is done using lowess regression.
pairs.panels(data.training, lm=TRUE) # lm=TRUE, linear regression fits are shown for both y by x and x by y. 

# Q2====
# fit a linear first-order regression model with SWI as outcome
# and SWF, temperature, size and management (not duration!) as predictors.
data.training <- data.training[, !names(data.training) == 'duration']
n_test <- dim(data.training)[1]
n_test 
p <- dim(data.training)[2]
p

# linear first-order regression model
fit <- lm(SWI ~ SWF+temperature+size+management, data = data.training)
summary(fit) # test whether a particular regression coefficient is significantly different from zero.
# ANOVA, test whether the regression model as a whole is performing significantly better than a null model
anova(fit) 
# SWF, temperature, management are significant, size non-significant
# Individual confidence intervals
alpha <- 0.05
confint(fit, level = 1 - alpha)
# Simultaneous confidence intervals with Bonferroni correction
alpha <- 0.05
confint(fit, level = 1 - alpha / 2)

# (a) Check whether a first-order model adequately captures the variability in the outcome====
# Multiple R-squared:  0.5805,	Adjusted R-squared:  0.5719
# R^2
summary(fit)
summary(fit)$r.squared
# 58% of the total variance in the outcome is explained by the first-order model

# (b) Check the Gauss-Markov conditions====
# Check model assumptions
fit.res <- residuals(fit)
fit.stdres <- stdres(fit)
fit.fittedvalues <- fitted.values(fit)
par(mfrow = c(2,2))
qqnorm(fit.stdres, main="")
qqline(fit.stdres)
plot(fit.res, xlab = "Index", ylab = "Residual")
plot(fit.fittedvalues, fit.res, xlab = "Fitted value", ylab = "Residual")
lines(lowess(fit.res ~ fit.fittedvalues), col = "red")
plot(fit.stdres, xlab = "Index", ylab = "Standardized residual", ylim = c(-3,3))
abline(h = -2.5, lty = 2)
abline(h = 2.5, lty = 2)
# UL: small deviations from normal distributed residuals
# UR: pattern indicates Homoscedasticity (no heteroscedastic errors)
# BL: curved band suggests linearity assumption is not satisfied
# BR: outliers
par(mfrow = c(2,2))
plot(SWF, fit$residuals, ylab = "Residual")
lines(lowess(fit$residuals ~ SWF), col = "red")
# plot indicates the linear model is defective (add quadratic terms)
plot(temperature, fit$residuals, ylab = "Residual")
lines(lowess(fit$residuals ~ temperature), col = "red")
# plot indicates the errors are heteroscedastic
plot(size, fit$residuals, ylab = "Residual")
lines(lowess(fit$residuals ~ size), col = "red")
plot(management, fit$residuals, ylab = "Residual")
lines(lowess(fit$residuals ~ management), col = "red")
par(mfrow = c(1,1))

# Gauss-Markov conditions tests
summary(gvlma.lm(fit))

# Checking the normality of the residuals
plot(fit, which = 2)
# Shapiro-Wilk test and Kolmogorov-Smirnov test Testing Normality
shapiro.test(residuals(fit))
LillieTest(residuals(fit))
check_normality(fit) # OK: Residuals appear as normally distributed 

# Checking the linearity of the relationship
plot(fit, which = 1)
# plot the relationship between the fitted values and the observed values for the outcome variable
# A straight line suggests that there’s nothing grossly wrong
plot(fit.fittedvalues, SWI, xlab = "Fitted Values", ylab = "Observed Values")
lines(lowess(SWI ~ fit.fittedvalues), col = 'red')
# for each individual predictor
par(mfrow = c(2,2))
# partial-residual plots, cannot contain interactions
termplot(fit, partial.resid = TRUE) 
crPlots(fit)
ceresPlots(fit) # less prone to leakage of nonlinearity among the predictors.
residualPlots(model = fit) # Adding SWF^2
# this function also reports the results of a bunch of curvature tests. 
# For a predictor variable X, this test is equivalent to adding a new predictor 
# to the model corresponding to X^2. If it comes up significant, it implies that 
# there is some nonlinear relationship between the variable and the residuals.
par(mfrow = c(1,1))

# Checking the homogeneity of variance
plot(fit, which = 3)
ncvTest(fit)
bptest(fit, ~ SWF + temperature + size + management) # there’s no violation of heteroskedasticity
coeftest(fit, vcov= hccm)
# if homogeneity of variance is violated, sandwich estimators is applied.
# Because the homogeneity of variance assumption wasn’t violated, 
# these t tests are pretty much identical to the former ones in the summary(fit)

# Checking independence, which we assume to be met
DurbinWatsonTest(fit, alternative="two.sided", data=data.training)
durbinWatsonTest(fit, alternative="two.sided", data=data.training)

# (c) Check whether there is (severe) multicollinearity====
# Correlation
corx <- cor
# small correlations between variables
# VIF: the largest VIF is larger than 10, or 
# if the mean of the VIF values is considerably larger than 1.
VIF <- diag(solve(corx))
max(VIF)
mean(VIF)
# Eigenvalues: A condition number nj > 30 is an indication for multicollinearity.
corx.eig <- eigen(corx)$values
corx.eig
sqrt(max(corx.eig)/corx.eig)
# indicating no multicollinearity

# (d) Check whether there are influential outliers====
plot(fit, which = 4)
plot(fit, which = 5)
plot(fit, which = 6)
# This function creates a “bubble” plot of Studentized residuals versus hat values,
# with the areas of the circles representing the observations proportional to Cook's distance. 
# Vertical reference lines are drawn at twice and three times the average hat value, 
# horizontal reference lines at -2, 0, and 2 on the Studentized-residual scale.
influencePlot(fit, main="influence Plot", sub="cook's distance")
# added-variable partial-regression plots
# identify data points with high leverage and 
# influential data points that might not have high leverage. 
par(mfrow = c(2,2)) 
avPlot(fit, variable = 'SWF')
avPlot(fit, variable = 'management')
avPlot(fit, variable = 'temperature')
avPlot(fit, variable = 'size')
par(mfrow = c(1,1)) 

# Classical approaches to find the vertical outliers and the leverage points====
# Standardized residuals
par(mfrow = c(1,1)) 
fit.stdres <- stdres(fit)
plot(fit.stdres, ylim = c(-4,4), ylab = "Standardized residuals")
abline(h = c(-2.5,2.5), col = "red")
label_x <- seq(1,200)
text(subset(fit.stdres,fit.stdres >2.5), labels=row.names(subset(data.training, fit.stdres > 2.5)),
     x = as.character(label_x[fit.stdres >2.5]), cex = 0.7, pos = 1)
text(subset(fit.stdres,fit.stdres < -2.5), labels=row.names(subset(data.training, fit.stdres < -2.5)),
     x = as.character(label_x[fit.stdres < -2.5]), cex = 0.7, pos = 1)
which(fit.stdres > 2.5 | fit.stdres < -2.5) 
# Studentized residuals
fit.studres <- studres(fit)
plot(fit.studres, ylim = c(-4,4), ylab = "Studentized residuals")
abline(h = c(-2.5,2.5), col = "red")
text(subset(fit.studres,fit.studres >2.5), labels=row.names(subset(data.training, fit.studres > 2.5)),
     x = as.character(label_x[fit.studres >2.5]), cex = 0.7, pos = 1)
text(subset(fit.studres,fit.studres < -2.5), labels=row.names(subset(data.training, fit.studres < -2.5)),
     x = as.character(label_x[fit.studres < -2.5]), cex = 0.7, pos = 1)
which(fit.studres > 2.5 | fit.studres < -2.5) 
# Classical approaches to find the leverage points
# Diagonal elements of hat matrix
fit.influence <- influence(fit)
plot(fit.influence$hat, ylab = "Diagonal elements of hat matrix")
h = 2*p/n_test
abline(h = h, col = "red")
text(subset(fit.influence$hat,fit.influence$hat > h), 
     labels=row.names(subset(data.training, fit.influence$hat > h)),
     x = as.character(label_x[fit.influence$hat > h]), cex = 0.7, pos = 1)
which(fit.influence$hat > h) 
# measures of influence
# DFFITS
fit.dffits <- dffits(fit)
plot(fit.dffits, ylab = "DFFITS")
h = 2*sqrt(p/n_test)
abline(h = h, col = "red")
text(subset(fit.dffits,fit.dffits > h), labels=row.names(subset(data.training, fit.dffits > h)),
     x = as.character(label_x[fit.dffits > h]), cex = 0.7, pos = 1)
which(fit.dffits > h) 
# Cook's distance
fit.Cd <- cooks.distance(fit)
plot(fit.Cd, ylab = "Cook's distance")
abline(h = 1, col = "red")
which(fit.Cd > 1)
# DFBETAS
fit.dfbetas <- dfbetas(fit)
plot(fit.dfbetas, ylab = "DFBETAS")
h = 2/sqrt(n_test)
abline(h = h, col = "red")
x = fit.dfbetas[,2] > h
text(subset(fit.dfbetas[,2], x), labels=row.names(subset(data.training, x)),
     x = data.frame(fit.dfbetas)[,1][x], cex = 0.7, pos = 4)
which(fit.dfbetas[,2] > h) 
# Outliers are not noticed by Cook's distance, but DFFITS and DFBETAS are more powerful.

# Bonferroni Outlier Test
outlierTest(fit) # No outliers with Bonferroni p < 0.05

# robust diagnostic plot====
# Reweighted LTS (maximal (50%) breakdown value)
par(mfrow = c(1,1)) 
RLTS <- ltsReg(SWI ~ SWF+temperature+size+management, data = data.training)
summary(RLTS)
summary(RLTS)$r.squared # 63%
# Note: It is strongly recommend using lmrob() instead of ltsReg!
lmrob <- lmrob(SWI ~ SWF+temperature+size+management, data = data.training)
summary(lmrob)
summary(lmrob)$r.squared # 60%
# rqq: Normal Q-Q plot of the standardized residuals;
# rindex: plot of the standardized residuals versus their index;
# rfit: plot of the standardized residuals versus fitted values;
# rdiag: regression diagnostic plot.
plot(RLTS, which = 'rqq') 
plot(RLTS, which = 'rindex') 
plot(RLTS, which = 'rfit') # No. 113, 151, 190, 198 --> OBS 218, 286, 371, 395 
plot(RLTS, which = 'rdiag')
# Diagnostic plot
RLTS.stdres <- RLTS$residuals/RLTS$scale
plot(RLTS$RD, RLTS.stdres, ylim = c(-5,5), 
     xlab = "Robust distance", ylab = "Standardized 50% LTS residuals", 
     main = 'Regression Diagnostic Plot')
v = sqrt(qchisq(0.975, p - 1))
abline(v = sqrt(qchisq(0.975, p - 1)), col = "red")
abline(h = c(-2.5,2.5), col = "red")
text(subset(RLTS.stdres,RLTS.stdres >2.5), labels=row.names(subset(data.training, RLTS.stdres > 2.5)),
     x = as.character(RLTS$RD[RLTS.stdres >2.5]), cex = 0.7, pos = 2)
text(subset(RLTS.stdres,RLTS.stdres < -2.5), labels=row.names(subset(data.training, RLTS.stdres < -2.5)),
     x = as.character(RLTS$RD[RLTS.stdres < -2.5]), cex = 0.7, pos = 2)
which(RLTS.stdres > 2.5 | RLTS.stdres < -2.5) # vertical outliers: OBS 218 286 371 395 
text(subset(RLTS.stdres, RLTS$RD > v), labels=row.names(subset(data.training, RLTS$RD > v)),
     x = as.character(RLTS$RD[RLTS$RD > v]), cex = 0.7, pos = 1)
which(RLTS$RD > v) # good leverage points: OBS 1, 3, 27
# RLTS (30% breakdown value)
RLTS2 <- ltsReg(SWI ~ I(SWF^2)+temperature+management, data = data.training, alpha = 0.7)
summary(RLTS2)
# Detection of outliers
plot(RLTS2, which = 'rqq') 
plot(RLTS2, which = 'rindex') 
plot(RLTS2, which = 'rfit')
plot(RLTS2, which = 'rdiag')
RLTS2.stdres <- RLTS2$residuals/RLTS2$scale
plot(RLTS2$RD, RLTS2.stdres, ylim = c(-5,5), 
     xlab = "Robust distance", ylab = "Standardized 30% LTS residuals", 
     main = 'Regression Diagnostic Plot')
v = sqrt(qchisq(0.975, p - 1))
abline(v = sqrt(qchisq(0.975, p - 1)), col = "red")
abline(h = c(-2.5,2.5), col = "red")
text(subset(RLTS2.stdres,RLTS2.stdres >2.5), 
     labels=row.names(subset(data.training, RLTS2.stdres > 2.5)),
     x = as.character(RLTS2$RD[RLTS2.stdres >2.5]), cex = 0.7, pos = 2)
text(subset(RLTS2.stdres,RLTS2.stdres < -2.5), 
     labels=row.names(subset(data.training, RLTS2.stdres < -2.5)),
     x = as.character(RLTS2$RD[RLTS2.stdres < -2.5]), cex = 0.7, pos = 2)
which(RLTS2.stdres > 2.5 | RLTS2.stdres < -2.5) 
text(subset(RLTS2.stdres, RLTS2$RD > v), labels=row.names(subset(data.training, RLTS2$RD > v)),
     x = as.character(RLTS2$RD[RLTS2$RD > v]), cex = 0.7, pos = 1)
which(RLTS2$RD > v) 

# Q3====
# Build a good linear regression model may containing higher-order terms, interactions, 
# transformed variables and/or other methods to improve the model assumptions.
# Model 1: Variable selection with Interaction terms====
# First look at the interaction terms. Generally the third and higher order interactions 
# are weak and hard to interpret, so look at the main effects and second order interactions. 
# The R formula syntax using ^2 to mean "all two-way interactions of the variables". 
fit_with <- lm(SWI ~ (SWF + temperature + management + size)^2, data = data.training)
summary(fit_with) # temperature:management interaction term is significant at the 5% level
# Backward elimination based on AIC
fit.full <- lm(SWI ~ (SWF + temperature + management + size)^2, data = data.training)
fit.full
stepAIC(fit.full, scope = list(upper = ~ (SWF + temperature + management + size)^2, lower = ~ 1), 
        direction = "backward")
# AIC=-339.91 to AIC=-348.44
# SWI ~ SWF + temperature + management + temperature:management
# Forward selection based on AIC
fit.null <- lm(SWI ~ 1, data = data.training)
fit.null
stepAIC(fit.null, scope = list(upper = ~ (SWF + temperature + management + size)^2, lower = ~ 1), 
        direction = "forward")
# AIC=-179.47 to  AIC=-348.44
# SWI ~ SWF + temperature + management + temperature:management
# Stepwise selection based on AIC (started at full model)
stepAIC(fit.full, scope=list(upper = ~ (SWF + temperature + management + size)^2, lower = ~ 1), 
        direction = "both")
# AIC=-339.91 to AIC=-348.44
# SWI ~ SWF + temperature + management + temperature:management
# Stepwise selection based on AIC (started at null model)
stepAIC(fit.null, scope=list(upper = ~ (SWF + temperature + management + size)^2, lower = ~ 1), 
        direction = "both")
# AIC=-179.47 to  AIC=-348.44
# SWI ~ SWF + temperature + management + temperature:management

fit_with <- lm(formula = SWI ~ SWF + temperature * management, data = data.training)
summary(fit_with) # temperature:management interaction term is only significant at the 10% level
# Reason 1: Many statisticians use a much larger significance level for the AB interaction F test than 
# what they use for the main effects. The reason is to get a higher chance to detect existing interactions. 
summary(fit_with)$r.squared # 0.5872287, 59% of the total variance in the outcome is explained
# Reason 2: stepAIC is equivalent to applying a hypothesis test with the significance level 0.157.
# According to the stepwise selection procedure, the model with the interaction term has smaller AIC
# as well as larger goodness-of-fit (R^2)
# Reason 3: As for interpretation, the longer being subject to nature management, the higher stability of
# nature area, the less infulence of the temperature change.
relweights <- function(fit, ...) {
  R <- cor(fit$model)
  nvar <- ncol(R)
  rxx <- R[2:nvar, 2:nvar]
  rxy <- R[2:nvar, 1]
  svd <- eigen(rxx)
  evec <- svd$vectors
  ev <- svd$values
  delta <- diag(sqrt(ev))
  # correlations between original predictors and new orthogonal variables
  lambda <- evec %*% delta %*% t(evec)
  lambdasq <- lambda^2
  # regression coefficients of Y on orthogonal variables
  beta <- solve(lambda) %*% rxy
  rsquare <- colSums(beta^2)
  rawwgt <- lambdasq %*% beta^2
  import <- (rawwgt/rsquare) * 100
  lbls <- names(fit$model[2:nvar])
  rownames(import) <- lbls
  colnames(import) <- "Weights"
  # plot results
  barplot(t(import), names.arg = lbls, ylab = "% of R-Square",
          xlab = "Predictor Variables", main = "Relative Importance of Predictor Variables",
          sub = paste("R-Square = ", round(rsquare, digits = 3)),
          ...)
  return(import)
}
relweights(fit, col = "lightgrey")
relweights(fit_with, col = "blue")
# size is dropped

# Model 2: Variable selection without Interaction terms====
# Backward elimination based on F-statistic/t-statistic
dropterm(fit.full, test = "F")
fit_drop <- update(fit.full, ~ . - temperature:size)
dropterm(fit_drop, test = "F")
fit_drop <- update(fit_drop, ~ . - temperature:size - SWF:management)
dropterm(fit_drop, test = "F")
fit_drop <- update(fit_drop, ~ . - temperature:size - SWF:management - SWF:size)
dropterm(fit_drop, test = "F")
fit_drop <- update(fit_drop, ~ . - temperature:size - SWF:management - SWF:size - management:size)
dropterm(fit_drop, test = "F")
fit_drop <- update(fit_drop, ~ . - temperature:size - SWF:management - SWF:size - management:size - size)
dropterm(fit_drop, test = "F")
fit_drop <- update(fit_drop, ~ . - temperature:size - SWF:management - SWF:size - management:size - size - SWF:temperature)
dropterm(fit_drop, test = "F")
fit_drop <- update(fit_drop, ~ . - temperature:size - SWF:management - SWF:size - management:size - size - SWF:temperature - temperature:management)
dropterm(fit_drop, test = "F")
# SWI ~ SWF + temperature + management
# Forward selection based on F-statistic/t-statistic
addterm(fit.null, ~ . + SWF + temperature + size + management + SWF:temperature + SWF:management + SWF:size +
          temperature:management + temperature:size + management:size, test = "F")
fit_add <- update(fit.null, ~ . + SWF)
addterm(fit_add, ~ . + temperature + size + management + SWF:temperature + SWF:management + SWF:size +
          temperature:management + temperature:size + management:size, test = "F")
fit_add <- update(fit_add, ~ . + temperature)
addterm(fit_add, ~. + size + management + SWF:temperature + SWF:management + SWF:size +
          temperature:management + temperature:size + management:size, test = "F")
fit_add <- update(fit_add, ~ . + management)
addterm(fit_add, ~. + size + SWF:temperature + SWF:management + SWF:size +
          temperature:management + temperature:size + management:size, test = "F")
# SWI ~ SWF + temperature + management

fit_without <- lm(formula = SWI ~ SWF + temperature + management, data = data.training)
summary(fit_without)
summary(fit_without)$r.squared # 0.5804086, 58% of the total variance in the outcome is explained
anova(fit_with, fit_without) # Reason 1: P = 0.07421, not significantly different between the two models
# −2log-likelihood+kn, where n represents the number of parameters in the fitted model, and k=2 for the usual AIC, or 
# k=log(N) (N being the number of observations) for the so-called BIC or SBC (Schwarz's Bayesian criterion)
stepAIC(fit_without, scope = list(upper = ~ SWF + temperature * management, lower = ~ 1), direction = "both")
AIC(fit_with, fit_without)
# Reason 2: AIC=-347.16 to AIC=-348.44, little increase in AIC value
# Reason 3: When an interaction isn’t significant, drop it if you are just checking for the presence of an interaction 
# to make sure you are specifying the model correctly. The interaction uses up df,
# changes the meaning of the lower order coefficients and complicates the model.
# But if you actually hypothesized an interaction that wasn’t significant, leave it in the model. 
# The insignificant interaction means something in this case – it helps you evaluate your hypothesis. 
# Taking it out can do more damage in specification error than in will in the loss of df.
#  Reason 4: Overall, no improvement on the assumptions of the model 1, which has a few outliers

# Model 3: Adding the Quadratic term of SWF====
fit3_with <- lm(SWI ~ SWF + I(SWF^2) + temperature * management, data = data.training)
summary(fit3_with) # SWF is non-significant, temperature:management is only significant at the 10% level
summary(fit3_with)$r.squared # 61% of the total variance in the outcome is explained
# Stepwise selection
fit.full <- lm(SWI ~ SWF + I(SWF^2) + temperature * management, 
               data = data.training)
fit.full
fit.null <- lm(SWI ~ 1, data = data.training)
fit.null
# Stepwise selection based on AIC (started at full model)
stepAIC(fit.full, scope=list(upper = ~ SWF + I(SWF^2) + temperature * management, lower = ~ 1), direction = "both")
# AIC=-357.84 to AIC=-359.58
# SWI ~ I(SWF^2) + temperature + management + temperature:management
# Stepwise selection based on AIC (started at null model)
stepAIC(fit.null, scope=list(upper = ~ SWF + I(SWF^2) + temperature * management, lower = ~ 1), direction = "both")
# AIC=-179.47 to AIC=-359.58
# SWI ~ I(SWF^2) + temperature + management + temperature:management
fit3_with <- lm(SWI ~ I(SWF^2) + temperature * management, data = data.training)
summary(fit3_with)
summary(fit3_with)$r.squared # 61% of the total variance in the outcome is explained

fit3_without <- lm(SWI ~ SWF + I(SWF^2) + temperature + management, data = data.training)
summary(fit3_without) # SWF is non-significant
summary(fit3_without)$r.squared # 60% of the total variance in the outcome is explained
# Stepwise selection
fit.full <- lm(SWI ~ SWF + I(SWF^2) + temperature + management, 
               data = data.training)
fit.full
fit.null <- lm(SWI ~ 1, data = data.training)
fit.null
# Stepwise selection based on AIC (started at full model)
stepAIC(fit.full, scope=list(upper = ~ SWF + I(SWF^2) + temperature + management, lower = ~ 1), direction = "both")
# AIC=-357.03 to AIC=-358.68
# SWI ~ I(SWF^2) + temperature + management
# Stepwise selection based on AIC (started at null model)
stepAIC(fit.null, scope=list(upper = ~ SWF + I(SWF^2) + temperature + management, lower = ~ 1), direction = "both")
# AIC=-179.47 to AIC=-358.68
# SWI ~ I(SWF^2) + temperature + management
fit3_without <- lm(SWI ~ I(SWF^2) + temperature + management, data = data.training)
summary(fit3_without)
summary(fit3_without)$r.squared # 60% of the total variance in the outcome is explained

# Model 4: Transformations====
# Box-Cox transformation on Y, one of the solutions to the problem of linearality ====
sum(data.training$SWI <= 0) # response should be strictly positive
par(mfrow = c(1,1))
out_without <- boxcox(SWI ~ I(SWF^2)+temperature+management, plotit = TRUE, data = data.training)
lambda_without <- out_without$x[which(out_without$y == max(out_without$y))]
lambda_without # lambda = 0.7878788
out_with <- boxcox(SWI ~ I(SWF^2)+temperature*management, plotit = TRUE, data = data.training)
lambda_with <- out_with$x[which(out_with$y == max(out_with$y))]
lambda_with # lambda = 0.8282828
# powerTransform uses the maximum likelihood-like approach of Box and Cox (1964) to select a transformatiion
# of a univariate or multivariate response for normality, linearity and/or constant variance.
powerTransform(fit3_with, family="bcnPower")
powerTransform(fit3_without, family="bcnPower")
# lambda is approximately equal to 1, no Box-cox transformation

# X Variable transformation of temperature====

# try segmented linear regression/Piecewise linear regression====
fit_without_segmented <- segmented.lm(fit3_without, seg.Z = ~temperature, psi = c(12, 27), control=seg.control(display=FALSE))
summary(fit_without_segmented)
# Estimated Break-Point(s): 12.748, 27.200
fit_without_segmented.res <- residuals(fit_without_segmented)
fit_without_segmented.stdres <- stdres(fit_without_segmented)
fit_without_segmented.fittedvalues <- fitted.values(fit_without_segmented)
par(mfrow = c(2,2))
qqnorm(fit_without_segmented.stdres, main="")
qqline(fit_without_segmented.stdres)
plot(fit_without_segmented.res, xlab = "Index", ylab = "Residual")
plot(fit_without_segmented.fittedvalues, fit_without_segmented.res, xlab = "Fitted value", ylab = "Residual")
lines(lowess(fit_without_segmented.res ~ fit_without_segmented.fittedvalues), col = "red")
plot(fit_without_segmented.stdres, xlab = "Index", ylab = "Standardized residual", ylim = c(-3,3))
abline(h = -2.5, lty = 2)
abline(h = 2.5, lty = 2)
par(mfrow = c(1,3))
plot(I(SWF^2), fit_without_segmented$residuals, ylab = "Residual")
lines(lowess(fit_without_segmented$residuals ~ I(SWF^2)), col = "red")
# plot indicates the errors are heteroscedastic
plot(temperature, fit_without_segmented$residuals, ylab = "Residual", main = 'Piecewise')
lines(lowess(fit_without_segmented$residuals ~ temperature), col = "red")
# plot indicates the linear model is defective (curve segmentation > 20) and the errors are heteroscedastic
plot(management, fit_without_segmented$residuals, ylab = "Residual")
lines(lowess(fit_without_segmented$residuals ~ management), col = "red")
par(mfrow = c(1,1))

fit_with_segmented <- segmented.lm(fit3_with, seg.Z = ~temperature, psi = c(12, 27), control=seg.control(display=FALSE))
# Estimated Break-Point(s): 12.621, 27.200
summary(fit_with_segmented)
fit_with_segmented.res <- residuals(fit_with_segmented)
fit_with_segmented.stdres <- stdres(fit_with_segmented)
fit_with_segmented.fittedvalues <- fitted.values(fit_with_segmented)
par(mfrow = c(2,2))
qqnorm(fit_with_segmented.stdres, main="")
qqline(fit_with_segmented.stdres)
plot(fit_with_segmented.res, xlab = "Index", ylab = "Residual")
plot(fit_with_segmented.fittedvalues, fit_with_segmented.res, xlab = "Fitted value", ylab = "Residual")
lines(lowess(fit_with_segmented.res ~ fit_with_segmented.fittedvalues), col = "red")
plot(fit_with_segmented.stdres, xlab = "Index", ylab = "Standardized residual", ylim = c(-3,3))
abline(h = -2.5, lty = 2)
abline(h = 2.5, lty = 2)
par(mfrow = c(2,2))
plot(I(SWF^2), fit_with_segmented$residuals, ylab = "Residual")
lines(lowess(fit_with_segmented$residuals ~ I(SWF^2)), col = "red")
# plot indicates the errors are heteroscedastic
plot(temperature, fit_with_segmented$residuals, ylab = "Residual")
lines(lowess(fit_with_segmented$residuals ~ temperature), col = "red")
# plot indicates the linear model is defective (curve segmentation > 20) and the errors are heteroscedastic
plot(management, fit_with_segmented$residuals, ylab = "Residual")
lines(lowess(fit_with_segmented$residuals ~ management), col = "red")
plot(temperature*management, fit_with_segmented$residuals, ylab = "Residual")
lines(lowess(fit_with_segmented$residuals ~ temperature*management), col = "red")
par(mfrow = c(1,1))

# link function====
boxTidwell(SWI ~ I(SWF^2) + temperature + abs(management+0.00000001),
           other.x = ~ size, data = data.training) 
# lambda = 0.5, but not significant with P = 0.2012
plot(SWI ~ temperature, data = data.training)
lines(lowess(SWI ~ temperature))
plot(SWI ~ sqrt(temperature), data = data.training)
lines(lowess(SWI ~ sqrt(temperature)))
plot(SWI ~ logit(temperature), data = data.training)
lines(lowess(SWI ~ logit(temperature)))
plot(SWI ~ I(temperature^2), data = data.training)
lines(lowess(SWI ~ I(temperature^2)))
plot(SWI ~ log(temperature), data = data.training)
lines(lowess(SWI ~ log(temperature)))
# try the logarithms transformation, logit transformation, square-root transformation 
# and even the quadratic term, in order to spread out the tails of the distribution. 
f1 <- lm(SWI ~ I(SWF^2)+I(temperature^2)+management, data = data.training)
f11 <- lm(SWI ~ I(SWF^2)+I(temperature^2)*management, data = data.training)
f2<- lm(SWI ~ I(SWF^2)+logit(temperature)+management, data = data.training)
f22<- lm(SWI ~ I(SWF^2)+logit(temperature)*management, data = data.training)
f3 <- lm(SWI ~ I(SWF^2)+log(temperature)+management, data = data.training)
f33 <- lm(SWI ~ I(SWF^2)+log(temperature)*management, data = data.training)
fit4_without <- lm(SWI ~ I(SWF^2)+sqrt(temperature)+management, data = data.training)
fit4_with <- lm(SWI ~ I(SWF^2)+sqrt(temperature)*management, data = data.training)

# Model 5: Weighted least squares model=====
fit4_without <- lm(formula = SWI ~ I(SWF^2) + temperature + management, data = data.training)
w_without <- 1/lm(abs(stdres(fit4_without)) ~  I(SWF^2)+ temperature +management, data = data.training)$fitted.values^2
fit5_nontrafo <- lm(SWI ~ I(SWF^2)+ temperature +management, weight = w_without, data = data.training)

fit4_without_trafo <- lm(formula = SWI ~ I(SWF^2) + sqrt(temperature) + management, data = data.training)
w_trafo <- 1/lm(abs(stdres(fit4_without_trafo)) ~  I(SWF^2)+ sqrt(temperature) +management, data = data.training)$fitted.values^2
fit5_trafo <- lm(SWI ~ I(SWF^2)+ sqrt(temperature) +management, weight = w_trafo, data = data.training)

fit4_with <- lm(formula = SWI ~ I(SWF^2) + temperature * management, data = data.training)
w_with <- 1/lm(abs(stdres(fit4_with)) ~  I(SWF^2)+temperature*management, data = data.training)$fitted.values^2
fit5_with_nontrafo <- lm(SWI ~ I(SWF^2)+temperature*management, weight = w_with, data = data.training)

fit4_with_trafo <- lm(formula = SWI ~ I(SWF^2) + sqrt(temperature)  * management, data = data.training)
w_with_trafo <- 1/lm(abs(stdres(fit4_with_trafo)) ~  I(SWF^2)+sqrt(temperature) *management, data = data.training)$fitted.values^2
fit5_with_trafo <- lm(SWI ~ I(SWF^2)+sqrt(temperature) *management, weight = w_with_trafo, data = data.training)

# Model 6: Boxcox transformation of Model 5====
out_without <- boxcox(fit5_nontrafo, plotit = TRUE)
lambda_without <- out_without$x[which(out_without$y == max(out_without$y))] 
lambda_without # sqrt(y)
out_without <- boxcox(fit5_trafo, plotit = TRUE)
lambda_without <- out_without$x[which(out_without$y == max(out_without$y))] 
lambda_without # sqrt(y)
fit6_nontrafo <- lm(SWI^0.5 ~ I(SWF^2)+temperature+management, 
                    weight = w_without, data = data.training)
fit6_trafo <- lm(SWI^0.5 ~ I(SWF^2)+sqrt(temperature)+management, 
                 weight = w_trafo, data = data.training)
# Check model assumptions

# Leave-one-out methods: PRESS
# Models with small PRESSp values (or PRESSp/n) are considered good candidate models
PRESS1 <- sum((residuals(fit5_nontrafo) / (1 - lm.influence(fit5_nontrafo)$hat))^2)
PRESS2 <- sum((residuals(fit5_trafo) / (1 - lm.influence(fit5_trafo)$hat))^2)
PRESS3 <- sum((residuals(fit6_nontrafo) / (1 - lm.influence(fit6_nontrafo)$hat))^2)
PRESS4 <- sum((residuals(fit6_trafo) / (1 - lm.influence(fit6_trafo)$hat))^2)
PRESS5 <- sum((residuals(fit5_with_nontrafo) / (1 - lm.influence(fit5_with_nontrafo)$hat))^2)
PRESS6 <- sum((residuals(fit5_with_trafo) / (1 - lm.influence(fit5_with_trafo)$hat))^2)
PRESS <- c(PRESS1, PRESS2,PRESS3, PRESS4,PRESS5, PRESS6)
names(PRESS) <- c("fit5_nontrafo", 'fit5_trafo','fit6_nontrafo',"fit6_trafo",'fit5_with_nontrafo','fit5_with_trafo')
sort(PRESS)

# MSE
MSE1 <- summary(fit5_nontrafo)$sigma^2
MSE2 <- summary(fit5_trafo)$sigma^2
MSE3 <- summary(fit6_nontrafo)$sigma^2
MSE4 <- summary(fit6_trafo)$sigma^2
MSE5 <- summary(fit5_with_nontrafo)$sigma^2
MSE6 <- summary(fit5_with_trafo)$sigma^2
MSE <- c(MSE1, MSE2, MSE3,MSE4, MSE5, MSE6)
names(MSE) <- c("fit5_nontrafo", 'fit5_trafo','fit6_nontrafo',"fit6_trafo",'fit5_with_nontrafo','fit5_with_trafo')
sort(MSE)
detach(data.training)
# Q4====
# Fit both models to the validation data. Investigate and compare their performance.
# model 5 fit5_nontrafo; model 6 fit6_nontrafo.
attach(data.test)

# Fit both models to the validation data. Investigate and compare their performance.
# model 6 fit6_trafo; model 6 fit6_nontrafo.
attach(data.test)

# Model fit6_nontrafo: SWI^0.5 ~ I(SWF^2) + temperature+management, weights = w_without
model5_OLS <- lm(SWI ~ I(SWF^2) + temperature + management, data = data.test)
w5_without <- 1/lm(abs(stdres(model5_OLS)) ~  I(SWF^2)+ temperature +management, data = data.test)$fitted.values^2
model6_nontrafo.val <- lm(SWI^0.5 ~ I(SWF^2) + temperature + management, weights = w5_without, data = data.test)
summary(model6_nontrafo.val); summary(fit6_nontrafo)
summary(model6_nontrafo.val)$r.squared
# Model fit6_trafo: SWI^0.5 ~ I(SWF^2) + sqrt(temperature)+management, weights = w_without
model5_trafo <- lm(SWI ~ I(SWF^2) + sqrt(temperature) + management, data = data.test)
w5_trafo <- 1/lm(abs(stdres(model5_trafo)) ~  I(SWF^2)+ sqrt(temperature) +management, data = data.test)$fitted.values^2
model6_trafo.val <- lm(SWI^0.5 ~ I(SWF^2) + sqrt(temperature) + management, weights = w5_trafo, data = data.test)
summary(model6_trafo.val); summary(fit6_trafo)
summary(model6_trafo.val)$r.squared

# Compare estimated coefficients and standard errors
# Individual confidence intervals
alpha <- 0.05
confint(fit6_trafo, level = 1 - alpha)
confint(model6_trafo.val, level = 1 - alpha)
confint(model6_nontrafo.val, level = 1 - alpha)
confint(fit6_nontrafo, level = 1 - alpha)
# Simultaneous confidence intervals with Bonferroni correction
alpha <- 0.05
confint(fit6_trafo, level = 1 - alpha/2)
confint(model6_trafo.val, level = 1 - alpha/2)
confint(model6_nontrafo.val, level = 1 - alpha/2)
confint(fit6_nontrafo, level = 1 - alpha/2)

# Prediction
# A prediction interval reflects the uncertainty of a single value, 
# while a confidence interval reflects the uncertainty of the predicted mean.
pred_trafo <- predict(fit6_trafo, newdata = data.test, interval = "prediction")
pred_nontrafo <- predict(fit6_nontrafo, newdata = data.test, interval = "prediction")
predict(fit6_trafo, newdata = data.test, interval = "confidence")
predict(fit6_nontrafo, newdata = data.test, interval = "confidence")
# MSEP
MSEP1 <- mean((predict(model6_trafo.val, newdata = data.test) - SWI)^2)
MSEP2 <- mean((predict(model6_nontrafo.val, newdata = data.test) - SWI)^2)
MSEP <- c(MSEP1, MSEP2)
names(MSEP) <- c("model6_trafo.val", "model6_nontrafo.val")
sort(MSEP)

# Leave-one-out methods: PRESS
# Models with small PRESSp values (or PRESSp/n) are considered good candidate models
PRESS1 <- sum((residuals(model6_trafo.val) / (1 - lm.influence(model6_trafo.val)$hat))^2)
PRESS2 <- sum((residuals(model6_nontrafo.val) / (1 - lm.influence(model6_nontrafo.val)$hat))^2)
PRESS <- c(PRESS1, PRESS2)
names(PRESS) <- c('model6_trafo.val',"model6_nontrafo.val")
sort(PRESS)
# MSE
MSE1 <- summary(model6_trafo.val)$sigma^2
MSE2 <- summary(model6_nontrafo.val)$sigma^2
MSE <- c(MSE1, MSE2)
names(MSE) <- c("model6_trafo.val", 'model6_nontrafo.val')
sort(MSE)

detach(data.test)

# Q5====
# fit your ultimate model (fit6_nontrafo) of preference to the full dataset.
attach(data.full)

# Model fit6_nontrafo: SWI^0.5 ~ I(SWF^2) + temperature+management, weights = w_without
fit_full <- lm(SWI ~ I(SWF^2) + temperature + management, data = data.full)
w_full <- 1/lm(abs(stdres(fit_full)) ~  I(SWF^2)+ temperature +management, data = data.full)$fitted.values^2
model6.full <- lm(SWI^0.5 ~ I(SWF^2) + temperature + management, weights = w_full, data = data.full)
summary(model6.full)
summary(model6.full)$r.squared

# Individual confidence intervals
alpha <- 0.05
confint(fit6_nontrafo, level = 1 - alpha)
confint(model6_nontrafo.val, level = 1 - alpha)
confint(model6.full, level = 1 - alpha)
# Simultaneous confidence intervals with Bonferroni correction
alpha <- 0.05
confint(fit6_nontrafo, level = 1 - alpha/2)
confint(model6_nontrafo.val, level = 1 - alpha/2)
confint(model6.full, level = 1 - alpha/2)

# Check model assumptions
model6.full.res <- residuals(model6.full)
model6.full.stdres <- stdres(model6.full)
model6.full.fittedvalues <- fitted.values(model6.full)
par(mfrow = c(2,2))
qqnorm(model6.full.stdres, main="")
qqline(model6.full.stdres)
plot(model6.full.res, xlab = "Index", ylab = "Residual")
plot(model6.full.fittedvalues, model6.full.res, xlab = "Fitted value", ylab = "Residual")
lines(lowess(model6.full.res ~ model6.full.fittedvalues), col = "red")
plot(model6.full.stdres, xlab = "Index", ylab = "Standardized residual", ylim = c(-3,3))
abline(h = -2.5, lty = 2)
abline(h = 2.5, lty = 2)
# UL: small deviations from normal distributed residuals
# UR: pattern indicates no heteroscedastic errors
# BL: linearity assumption is satisfied
# BR: outliers
par(mfrow = c(1,3))
plot(I(SWF^2), model6.full$residuals, ylab = "Residual")
lines(lowess(model6.full$residuals ~ I(SWF^2)), col = "red")
plot(temperature, model6.full$residuals, ylab = "Residual")
lines(lowess(model6.full$residuals ~ temperature), col = "red")
plot(management, model6.full$residuals, ylab = "Residual")
lines(lowess(model6.full$residuals ~ management), col = "red")
par(mfrow = c(1,1))

# Checking the normality of the residuals
plot(model6.full, which = 2)
# Shapiro-Wilk test and Kolmogorov-Smirnov test Testing Normality
shapiro.test(residuals(model6.full))
LillieTest(residuals(model6.full)) 
sf.test(residuals(model6.full))  
check_normality(model6.full)# OK: Residuals appear as normally distributed 
# Checking the linearity of the relationship
plot(model6.full, which = 1)
# plot the relationship between the fitted values and the observed values for the outcome variable
plot(model6.full.fittedvalues, SWI, xlab = "Fitted Values", ylab = "Observed Values")
lines(lowess(SWI ~ model6.full.fittedvalues), col = 'red')
# for each individual predictor
par(mfrow = c(1,3))
# partial-residual plots, cannot contain interactions
termplot(model6.full, partial.resid = TRUE) 
crPlots(model6.full)
par(mfrow = c(1,1))
# Checking the homogeneity of variance
# https://stats.stackexchange.com/questions/193061/what-is-the-difference-between-these-two-breusch-pagan-tests
# In short, the studentized BP test is more robust, usually go with bptest, 
# with studentize = TRUE (default) and varformula = ~ fitted.values(my.lm) as options, 
# for an initial approach for homoskedasticity.
plot(model6.full, which = 3)
ncvTest(model6.full) 
bptest(model6.full, ~ SWF + temperature + management) 
bptest(model6.full, ~ fitted.values(model6.full)) # accepted
coeftest(model6.full, vcov= hccm)
summary(model6.full)
# if homogeneity of variance is violated, sandwich estimators is applied.
# Because the homogeneity of variance assumption wasn’t violated, 
# these t tests are pretty much identical to the former ones in the summary(model6.full)

# outliers
plot(model6.full, which = 4)
plot(model6.full, which = 5)
influencePlot(model6.full, main="influence Plot", sub="cook's distance")
par(mfrow = c(1,3))
avPlot(model6.full, variable = 'I(SWF^2)')
avPlot(model6.full, variable = 'management')
avPlot(model6.full, variable = 'temperature')
par(mfrow = c(1,1))

# Standardized residuals
model6.full.stdres <- stdres(model6.full)
plot(model6.full.stdres, ylim = c(-4,4), ylab = "Standardized residuals")
abline(h = c(-2.5,2.5), col = "red")
label_x <- seq(1,400)
text(subset(model6.full.stdres,model6.full.stdres >2.5), labels=row.names(subset(data.full, model6.full.stdres > 2.5)),
     x = as.character(label_x[model6.full.stdres >2.5]), cex = 0.7, pos = 1)
text(subset(model6.full.stdres,model6.full.stdres < -2.5), labels=row.names(subset(data.full, model6.full.stdres < -2.5)),
     x = as.character(label_x[model6.full.stdres < -2.5]), cex = 0.7, pos = 1)
which(model6.full.stdres > 2.5 | model6.full.stdres < -2.5) 
# Studentized residuals
model6.full.studres <- studres(model6.full)
plot(model6.full.studres, ylim = c(-4,4), ylab = "Studentized residuals")
abline(h = c(-2.5,2.5), col = "red")
text(subset(model6.full.studres,model6.full.studres >2.5), labels=row.names(subset(data.full, model6.full.studres > 2.5)),
     x = as.character(label_x[model6.full.studres >2.5]), cex = 0.7, pos = 1)
text(subset(model6.full.studres,model6.full.studres < -2.5), labels=row.names(subset(data.full, model6.full.studres < -2.5)),
     x = as.character(label_x[model6.full.studres < -2.5]), cex = 0.7, pos = 1)
which(model6.full.studres > 2.5 | model6.full.studres < -2.5)
# Bonferroni Outlier Test
outlierTest(model6.full) # OBS 1 as a outlier with Bonferroni p < 0.05

detach(data.full)

# Q6====
# investigating possible association between duration (outcome) and temperature (predictor).
data.training <- data.full[-d.test, ]
attach(data.training)

# (a) Fit non-parametric models with k=1 and k=2, ====
# for spans 0.25, 0.5, and 0.75 and choose the best-fitting model
# Local linear regression
plot(temperature, duration, main = "Local linear regression")
s <- c(0.25, 0.5, 0.75)
colors <- c("red", "green", "blue")
for (i in 1:length(s)) lines(temperature, predict(loess(duration ~ temperature, span = s[i], degree = 1), 
                                                  data = data.training), col = colors[i])
legend(5, 40, c("span = 0.25", "span = 0.5", "span = 0.75"), lty = 1, col = colors)
# Local quadratic regression
plot(temperature, duration, main = "Local quadratic regression")
for (i in 1:length(s)) lines(temperature, predict(loess(duration ~ temperature, span = s[i], degree = 2), 
                                                  data = data.training), col = colors[i])
legend(5, 40, c("span = 0.25", "span = 0.5", "span = 0.75"), lty = 1, col = colors)

# Check model assumptions
# ========== k=1, span=0.25
fit.loess1<-loess(duration~temperature,degree = 1,span=0.25)
fit.loess1
par(mfrow=c(2,2))
plot(duration~temperature)
lines(loess.smooth(temperature,duration,span=0.25,degree=1),col='red')
plot(residuals(fit.loess1)~temperature)
lines(loess.smooth(temperature,residuals(fit.loess1),span=0.25,degree=1),col='red')
abline(h=0,lty=2)
plot(fitted(fit.loess1),sqrt(abs(residuals(fit.loess1))))
lines(loess.smooth(fitted(fit.loess1),sqrt(abs(residuals(fit.loess1))),span=0.25,degree=1),col='red')
qqnorm(residuals(fit.loess1))
qqline(residuals(fit.loess1),col='red')
# ========== k=1, span=0.5
fit.loess2<-loess(duration~temperature,degree = 1,span=0.5)
fit.loess2
par(mfrow=c(2,2))
plot(duration~temperature)
lines(loess.smooth(temperature,duration,span=0.5,degree=1),col='red')
plot(residuals(fit.loess2)~temperature)
lines(loess.smooth(temperature,residuals(fit.loess2),span=0.5,degree=1),col='red')
abline(h=0,lty=2)
plot(fitted(fit.loess2),sqrt(abs(residuals(fit.loess2))))
lines(loess.smooth(fitted(fit.loess2),sqrt(abs(residuals(fit.loess2))),span=0.5,degree=1),col='red')
qqnorm(residuals(fit.loess2))
qqline(residuals(fit.loess2),col='red')
# ========== k=1, span=0.75
fit.loess3<-loess(duration~temperature,degree = 1,span=0.75)
fit.loess3
par(mfrow=c(2,2))
plot(duration~temperature)
lines(loess.smooth(temperature,duration,span=0.75,degree=1),col='red')
plot(residuals(fit.loess3)~temperature)
lines(loess.smooth(temperature,residuals(fit.loess3),span=0.75,degree=1),col='red')
abline(h=0,lty=2)
plot(fitted(fit.loess3),sqrt(abs(residuals(fit.loess3))))
lines(loess.smooth(fitted(fit.loess3),sqrt(abs(residuals(fit.loess3))),span=0.75,degree=1),col='red')
qqnorm(residuals(fit.loess3))
qqline(residuals(fit.loess3),col='red')
# ========== k=2, span=0.25
fit.loess4<-loess(duration~temperature,degree = 2,span=0.25)
fit.loess4
par(mfrow=c(2,2))
plot(duration~temperature)
lines(loess.smooth(temperature,duration,span=0.25,degree=2),col='red')
plot(residuals(fit.loess4)~temperature)
lines(loess.smooth(temperature,residuals(fit.loess4),span=0.25,degree=2),col='red')
abline(h=0,lty=2)
plot(fitted(fit.loess4),sqrt(abs(residuals(fit.loess4))))
lines(loess.smooth(fitted(fit.loess4),sqrt(abs(residuals(fit.loess4))),span=0.25,degree=2),col='red')
qqnorm(residuals(fit.loess4))
qqline(residuals(fit.loess4),col='red')
# ========== k=2, span=0.5
fit.loess5<-loess(duration~temperature,degree = 2,span=0.5)
fit.loess5
par(mfrow=c(2,2))
plot(duration~temperature)
lines(loess.smooth(temperature,duration,span=0.5,degree=2),col='red')
plot(residuals(fit.loess5)~temperature)
lines(loess.smooth(temperature,residuals(fit.loess5),span=0.5,degree=2),col='red')
abline(h=0,lty=2)
plot(fitted(fit.loess5),sqrt(abs(residuals(fit.loess5))))
lines(loess.smooth(fitted(fit.loess5),sqrt(abs(residuals(fit.loess5))),span=0.5,degree=2),col='red')
qqnorm(residuals(fit.loess5))
qqline(residuals(fit.loess5),col='red')
# ========== k=2, span=0.75
fit.loess6<-loess(duration~temperature,degree = 2,span=0.75)
fit.loess6
par(mfrow=c(2,2))
plot(duration~temperature)
lines(loess.smooth(temperature,duration,span=0.75,degree=2),col='red')
plot(residuals(fit.loess6)~temperature)
lines(loess.smooth(temperature,residuals(fit.loess6),span=0.75,degree=2),col='red')
abline(h=0,lty=2)
plot(fitted(fit.loess6),sqrt(abs(residuals(fit.loess6))))
lines(loess.smooth(fitted(fit.loess6),sqrt(abs(residuals(fit.loess6))),span=0.75,degree=2),col='red')
qqnorm(residuals(fit.loess6))
qqline(residuals(fit.loess6),col='red')
# k =2, span = 0.75 is the best-fitting model

# (b) Fit a quadratic linear model====
par(mfrow = c(1,1))
plot(temperature, duration, main = "Polynomial regression")
# Linear model
fit1 <- lm(duration ~ temperature, data = data.training)
abline(fit1, col = "red")
# Quadratic model
fit2 <- lm(duration ~ temperature + I(temperature^2), data = data.training)
fit2.coef <- fit2$coefficients
curve(fit2.coef[1] + fit2.coef[2]*x + fit2.coef[3]*x^2, 0, 60, add = TRUE, col = "green")
# Cubic model
fit3 <- lm(duration ~ temperature + I(temperature^2) + I(temperature^3), data = data.training)
fit3.coef <- fit3$coefficients
curve(fit3.coef[1] + fit3.coef[2]*x + fit3.coef[3]*x^2 + fit3.coef[4]*x^3, 0, 60, add = TRUE, col = "blue")
# Add legend
legend(5, 40, c("linear", "quadratic", "cubic"), lty = 1, col = c("red", "green", "blue"))

# (c) a plot of the data, the best nonparametric fit, and the linear fit====
par(mfrow = c(1,1))
plot(temperature, duration, main = "Quadratic vs non-parametric regression")
# Local quadratic regression: k = 2, span = 0.75
fit.loess <- loess(duration ~ temperature, span = 0.75, degree = 2)
lines(temperature, predict(fit.loess, data = data.training), col = 'red')
# Linear model
fit1 <- lm(duration ~ temperature, data = data.training)
abline(fit1, col = "black")
# Quadratic model
fit2 <- lm(duration ~ temperature + I(temperature^2), data = data.training)
fit2.coef <- fit2$coefficients
curve(fit2.coef[1] + fit2.coef[2]*x + fit2.coef[3]*x^2, 0, 60, add = TRUE, col = "green")
# Cubic model
fit3 <- lm(duration ~ temperature + I(temperature^2) + I(temperature^3), data = data.training)
fit3.coef <- fit3$coefficients
curve(fit3.coef[1] + fit3.coef[2]*x + fit3.coef[3]*x^2 + fit3.coef[4]*x^3, 0, 60, add = TRUE, col = "blue")

legend(4, 43, c("Non-parametric fit", "linear", "quadratic", "cubic"), lty = 1, col = c("red", "black", "green", "blue"))
# According to the visual interpretation, Non-parametric model fits the data best.

# (d) Test whether the non-parametric model of your choice==== 
# fits the data better than the quadratic model
summary(fit.loess)
summary(fit2) # 69%
# Compare quadratic linear model with non-parametric model
traceS <- fit.loess$trace.hat
SSE0 <- sum(residuals(fit2)^2)
SSE1 <- sum(residuals(fit.loess)^2)
n <- dim(data.training)[1]
Fvalue <- ((SSE0 - SSE1) / (traceS - 3)) / (SSE1 / (n - traceS))
Fvalue
Fcrit <- qf(0.95, traceS - 3, n - traceS)
Fcrit
1 - pf(Fvalue, traceS - 3, n - traceS)
# the difference between the non-parametric model and the quadratic model is s
# ignificant since P-value is zero

# Prediction
attach(data.test)
t.pred <- predict(fit.loess, data.test, se = TRUE)
t.upper <- t.pred$fit + qnorm(0.975) * t.pred$se.fit
t.lower <- t.pred$fit - qnorm(0.975) * t.pred$se.fit
loess <- data.frame("pred" = t.pred$fit, "lower" = t.lower, "upper" = t.upper)
plot(data.test$temperature, data.test$duration)
lines(lowess(data.test$temperature,t.pred$fit))
lines(lowess(data.test$temperature,t.upper))
lines(lowess(data.test$temperature,t.lower))
t.pred <- predict(fit2, data.test, se = TRUE)
t.upper <- t.pred$fit + qnorm(0.975) * t.pred$se.fit
t.lower <- t.pred$fit - qnorm(0.975) * t.pred$se.fit
quadratic <- data.frame("pred" = t.pred$fit, "lower" = t.lower, "upper" = t.upper)
plot(data.test$temperature, data.test$duration)
lines(lowess(data.test$temperature,t.pred$fit))
lines(lowess(data.test$temperature,t.upper))
lines(lowess(data.test$temperature,t.lower))
detach(data.test)
# Assessing goodness of fit
# R-squared
rsq <- function (x, y) cor(x, y) ^ 2
rsq1 <- rsq(loess[,1], duration) # r.squared = 0.8168894
rsq2 <- rsq(quadratic[,1], duration) # r.squared = 0.6893309
RSQ <- c(rsq1, rsq2)
names(RSQ) <- c("Non-parametric", "Linear quadratic")
sort(RSQ)
# Residual sum of squares
RSS1 <- sum(residuals(fit.loess)^2)
RSS2 <- sum(residuals(fit2)^2)
RSS <- c(RSS1, RSS2)
names(RSS) <- c("Non-parametric", "Linear quadratic")
sort(RSS)
# Pearson estimated residual variance
sigma.squared1 <- RSS1 / (n - traceS)
sigma.squared2 <- RSS2 / fit2$df.residual
sigma.squared <- c(sigma.squared1, sigma.squared2)
names(sigma.squared) <- c("Non-parametric", "Linear quadratic")
sort(sigma.squared)
# Mean squared error
MSE1 <- sum(residuals(fit.loess)^2) / (n - traceS)
MSE2 <- sum(residuals(fit2)^2) / (fit2$df.residual)
MSE <- c(MSE1, MSE2)
names(MSE) <- c("Non-parametric", "Linear quadratic")
sort(MSE)
# Root mean squared error
RMSE1 <- sqrt(MSE1)
RMSE2 <- sqrt(MSE2)
RMSE <- c(RMSE1, RMSE2)
names(RMSE) <- c("Non-parametric", "Linear quadratic")
sort(RMSE)
# MSEP
MSEP1 <- mean((loess[,1] - duration)^2)
MSEP2 <- mean((quadratic[,1] - duration)^2)
MSEP <- c(MSEP1, MSEP2)
names(MSEP) <- c("Non-parametric", "Linear quadratic")
sort(MSEP)
compare.results <- data.frame(rbind(RSQ,RSS,sigma.squared, MSE, RMSE, MSEP), 
                              row.names = c('RSQ','RSS','sigma.squared', 'MSE', 'RMSE', 'MSEP'))
names(compare.results) <- c("Non-parametric", "Linear quadratic")
compare.results
# Non-parametric model fits the data better than the quadratic model
caret::postResample(loess[,1], duration)
caret::postResample(quadratic[,1], duration)

detach(data.training)
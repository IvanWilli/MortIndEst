### co-gompertz

# example from the package. pick one year
library(MortCast)
data(mxM, mxF, package = "wpp2017")
country <- "South Africa"
mxm <- subset(mxM, name == country)[,10]
mxf <- subset(mxF, name == country)[,10]
ages <- c(0, 1, seq(5, 100, 5))
est.ages = seq(70, 85, by = 5)
proj.ages = seq(90, 130, by = 5)

# co-kannisto
mxnew <- cokannisto(mxm, mxf, est.ages, proj.ages)
mx_female_cokannisto <- mxnew$female
mx_male_cokannisto <- mxnew$male

# co-gompertz (simple code before turning into a function)
cogompertz <- function(mxF, mxM, est.ages, proj.ages, 
                       method = "coherent", age_conv){
     
     # mxF <- mxf
     # mxM <- mxm
     # age_conv = 120
     
     if(is.null(rownames(mxF))){
          ages_all <- c(0, 1, seq(5, 150, 5))
          ages <- as.integer(ages_all[1:length(mxF)])
     }else{
          ages <- as.integer(names(mxM))     
     }
     mx <- c(mxM, mxF)
     y <- log(mx)
     x <- c(ages, ages)
     g <- c(rep(1, length(mxM)), rep(0, length(mxM)))
     df <- data.frame(y = y, g = g, x = x)
     if(method == "coherent"){
          coefs <- coefficients(lm(y ~ g + x, data = df[df$x %in% est.ages,]))
          female.coefs <- c(c = exp(coefs[[1]]), d = coefs[["x"]])
          male.coefs <- c(c = exp(coefs[[1]] + coefs[["g"]]), d = coefs[["x"]])
          female_mx <- female.coefs["c"] * exp(female.coefs["d"] * proj.ages)
          male_mx <- male.coefs["c"] * exp(male.coefs["d"] * proj.ages)
     }
     if(method == "convergent"){
          female.coefs <- coefficients(lm(y ~ x, data = df[df$x %in% est.ages & df$g == 0,]))
          female_mx <- female.coefs[1] + female.coefs[2] * proj.ages
          init_age <- proj.ages[1] - 5
          df_male_star <- data.frame(x = c(init_age, age_conv),
                                     y = c(df$y[df$x %in% init_age & df$g == 1], 
                                           female_mx[proj.ages == age_conv]))
          lm_male_star <- lm(y ~ x, data = df_male_star)
          male_mx <- predict(lm_male_star, newdata = data.frame(x = proj.ages))
          male_mx[proj.ages>age_conv] <- female_mx[proj.ages>age_conv]
          # plot(ages, log(mxF), xlim = c(0, 130), ylim = c(-8,2)); lines(proj.ages, female_mx)
          # points(ages, log(mxM), col = 2); lines(proj.ages, male_mx, col = 2)
          female_mx <- exp(female_mx)
          male_mx <- exp(male_mx)
     }
     mx_female_cogompertz <- c(mxF[!ages %in% proj.ages], female_mx)
     mx_male_cogompertz <- c(mxM[!ages %in% proj.ages], male_mx)
     mxnew <- data.frame(x = unique(c(ages, proj.ages)), 
                         female = mx_female_cogompertz, male = mx_male_cogompertz)
     return(mxnew)
}

# co-kannisto
mxnew <- cogompertz(mxf, mxm, est.ages, proj.ages, method = "convergent", age_conv = 130)
mx_female_cogompertz <- mxnew$female
mx_male_cogompertz <- mxnew$male

# plot
plot(ages, mxm, log = "y", xlim = c(0, 130), ylim = c(1e-3,3), col = 4, ylab = "5Mx", xlab = "Age")
points(ages, mxf, col = 2)
ages_full <- c(0, 1, seq(5, max(proj.ages), by=5))
lines(ages_full, mx_male_cokannisto, col = 4, lty = 2)
lines(ages_full, mx_female_cokannisto, col = 2, lty = 2)
lines(ages_full, mx_male_cogompertz, col = 4)
lines(ages_full, mx_female_cogompertz, col = 2)
abline(v = 90, lty = 2, col = "grey")
abline(h = 1, lty = 2, col = "grey")
legend("topleft", legend=c("male obs", "female obs",
                               "male cogompertz", "female cogompertz",
                               "male cokannisto", "female cokannisto"), 
       bty = "n",
       pch=c(1,1,rep(NA,4)),
       col=rep(c("blue", "red"),3), 
       lty=c(NA,NA,1,1,2,2))

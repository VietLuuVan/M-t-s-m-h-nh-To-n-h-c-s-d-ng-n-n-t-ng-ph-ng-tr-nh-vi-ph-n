library(deSolve)
sir_1 <- function(beta,alpha, gamma, S0,E0, I0, R0, times) {
  require(deSolve) # for the "ode" function
  
  # the differential equations:
  sir_equations <- function(time, variables, parameters) {
    with(as.list(c(variables, parameters)), {
      dS <- -(beta/8330083) * I * S
      dE <- (beta/8330083) *I * S - alpha * E
      dI <- alpha * E - gamma * I
      dR <-  gamma * I
      return(list(c(dS,dE, dI, dR)))
    })
  }
  
  # the parameters values:
  parameters_values <- c(beta  = beta, alpha = alpha, gamma = gamma)
  
  # the initial values of variables:
  initial_values <- c(S = S0,E = E0, I = I0, R = R0)
  
  # solving
  out <- ode(initial_values, times, sir_equations, parameters_values)
  
  # returning the output:
  as.data.frame(out)
}
sir_values_1 = sir_1(beta = 0.28,alpha = 0.03,gamma = 0.045,S0 = 8310127, E0 = 0, I0 = 11087, R0 = 8869, times = 0:365)
View(sir_values_1)
with(sir_values_1, {
  # plotting the time series of susceptibles:
  plot(time, S, type = "l", col = "blue",
       xlab = "time (days)", ylab = "number of people", ylim = c(0, 8000000), lwd = 3)
  # adding the time series of infectious:
  lines(time, E, col = "yellow", lwd = 3)
  # adding the time series of infectious:
  lines(time, I, col = "red", lwd = 3)
  # adding the time series of recovered:
  lines(time, R, col = "green", lwd = 3)
})

# adding a legend:
legend(200,6000000, c("susceptibles","exposed", "infectious", "recovered"),
       col = c("blue","yellow", "red", "green"),lwd = 3, lty = 1, bty = "o")
#Comparing a model's predictions with data
dat = read.table("C:\\Users\\thoidaipc\\Documents\\statistics\\sir.txt", header = T)
View(dat)
attach(dat)
model_fit <- function(beta,alpha, gamma, dat, N = 8330083, ...) {
  I0 <- cases[1] # initial number of infected (from data)
  times <- day  # time points (from data)
  # model's predictions:
  predictions <- sir_1(beta = beta,alpha = alpha, gamma = gamma,   # parameters
                       S0 = 8310127, E0 = 3326, I0 = 7760, R0 = 8869, # variables' intial values
                       times = times)                # time points
  # plotting the observed prevalences:
  print(predictions$I[50]-dat$cases[50])
  print(predictions$I[60]-dat$cases[60])
  print(predictions$I[80]-dat$cases[80])
  print(predictions$I[90]-dat$cases[90])
  plot(times, dat$cases, col = "red",type = "l", xlab = "days",ylab = "cases", lwd = 4
       , ylim = c(0,500000))
  # adding the model-predicted prevalence:
  lines(times, predictions$I, col = "blue", lwd = 4)
  legend(0,500000,c("data","prediction"), col = c("red","blue"),lwd = c(3,3), lty = 1, bty = "p")
}
model_fit(beta = 0.28,alpha = 0.028, gamma = 0.04, dat, pch = 19, col = "red")

#Estimating model's parameters
dat2 = read.table("C:\\Users\\thoidaipc\\Documents\\statistics\\est.txt", header = T)
dat2
ss <- function(beta, gamma, data = dat2, N = 8000000) {
  I0 <- data$cases[1]
  times <- data$day
  predictions <- sir_1(beta = beta, gamma = gamma,   # parameters
                       S0 = 7980044, I0 = 11087, R0 = 8869, # variables' intial values
                       times = times)                # time points
  sum((predictions$I[-1] - data$cases[-1])^2)
}
ss(beta = 0.132, gamma = 0.09091)
beta_val <- seq(from = 0.132, to = 0.135, le = 10)
gamma_val <- seq(from = 0.01, to = 0.09, le = 10)
min = ss(beta = 0.13, gamma = 0.09)
par = c(beta,gamma)
for (i in beta_val) {
  for (j in gamma_val) {
    if (min < ss(beta =i, gamma = j)) {
      min = ss(beta =i, gamma = j)
      par = c(i,j)
    }
  }
}
par
min_ss_val <- min(ss_val)
min_ss_val
beta_hat <- beta_val[ss_val == min_ss_val]
beta_hat
plot(beta_val, ss_val, type = "l", lwd = 2,
     xlab = expression(paste("infectious contact rate ", beta)),
     ylab = "sum of squares")
# adding the minimal value of the sum of squares:
abline(h = min_ss_val, lty = 2, col = "grey")
# adding the estimate of beta:
abline(v = beta_hat, lty = 2, col = "grey")
gamma_val <- seq(from = 0.01, to = 0.1, le = 100)
ss_val <- sapply(gamma_val, function(x) ss(beta_hat, x))
(min_ss_val <- min(ss_val))
(gamma_hat <- gamma_val[ss_val == min_ss_val])
plot(gamma_val, ss_val, type = "l", lwd = 2,
     xlab = expression(paste("recovery rate ", gamma)),
     ylab = "sum of squares")
abline(h = min_ss_val, lty = 2, col = "grey")
abline(v = gamma_hat, lty = 2, col = "grey")

n <- 10 # number of parameter values to try
beta_val <- seq(from = 0.002, to = 0.0035, le = n)
gamma_val <- seq(from = 0.3, to = 0.65, le = n)
param_val <- expand.grid(beta_val, gamma_val)
ss_val <- with(param_val, Map(ss, Var1, Var2))
ss_val <- matrix(unlist(ss_val), n)
persp(beta_val, gamma_val, -ss_val, theta = 40, phi = 30,
      xlab = "beta", ylab = "gamma", zlab = "-sum of squares")

n <- 30 # number of parameters values
beta_val <- seq(from = 0.002, to = 0.0035, le = n)
gamma_val <- seq(from = 0.3, to = 0.65, le = n)
# calculating the sum of squares:
param_val <- expand.grid(beta_val, gamma_val)
ss_val <- with(param_val, Map(ss, Var1, Var2))


ss_val <- unlist(ss_val)

# minimum sum of squares and parameters values:
(ss_min <- min(ss_val))
ind <- ss_val == ss_min
(beta_hat <- param_val$Var1[ind])
(gamma_hat <- param_val$Var2[ind])

ss_val <- matrix(ss_val, n)
image(beta_val, gamma_val, ss_val,
      xlab = expression(paste("infectious contact rate ", beta, " (/person/day)")),
      ylab = expression(paste("recovery rate ", gamma, " (/day)")))
contour(beta_val, gamma_val,ss_val, add = TRUE)
points(beta_hat, gamma_hat, pch = 3)
box(bty = "o")

#Estimate beta and gamma
ss2 <- function(x) {
  ss(beta = x[1], gamma = x[2])
}
starting_param_val <- c(0.004, 0.5)
ss_optim <- optim(starting_param_val, ss2)
ss_optim

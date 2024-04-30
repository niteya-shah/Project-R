library('invgamma')
library('truncnorm')

T = 1000

phi_0 = 0.96
sigma2_0= 0.002
mu_0 = -0.03
burnin_int = 10
c = 0.001

q_j = c(0.00730, 0.10556, 0.00002, 0.04395, 0.34001, 0.24566, 0.25750)
b_j = c(-10.12999, -3.397281, -8.56686, 2.77786, 0.61942, 1.79518, -1.08819) - 1.2704
w_j = c(5.79596, 2.61369, 5.17950, 0.16735, 0.64009, 0.34023, 1.26261)
J = length(q_j)


sigma2_n = sigma2_0
phi = phi_0
mu = mu_0

data = read.csv('DEXUSUK_2.csv')
length(data[,1])
stock_price = ts(data[1:(T+1),2])
plot(stock_price)
y_t = log(stock_price[2:(T+1)]) - log(stock_price[1:T]) - mean(log(stock_price[2:(T+1)]) - log(stock_price[1:T]))
y_t = y_t * 100

plot(y_t, type='l')

sigma2_correct = 0.0018
mu_correct = -0.01
phi_correct = 0.97

y_t = array(0, T + 1)
h_t_correct = array(0, T + 1)
h_t_correct = rnorm(1, mu_correct, sqrt(sigma2_correct/(1 - phi_correct**2)))
y_t[1] = exp(h_t_correct[1]/2) * rnorm(1, 0, 1)
for(i in 2:(T + 1))
{
  h_t_correct[i] = mu_correct + phi_correct * (h_t_correct[i - 1] - mu_correct) + rnorm(1, 0, sqrt(sigma2_correct))
  y_t[i] = exp(h_t_correct[i]/2) * rnorm(1, 0, 1)
}

mu_calc = mean(h_t_correct)

par(mfrow=c(2,1))
plot(y_t, type='l')
y_star_t = log(y_t**2 + c)
plot(y_star_t, type='l')

h_t = array(0.0, T)

sample_s = function(y_star_t, h_t)
{
  s = array(0, T)
  q_j_h = array(0, J)
  for(i in 1:T)
  {
    for(j in 1:J)
    {
      q_j_h[j] = q_j[j] * dnorm(y_star_t[i], h_t[i] + b_j[j], sqrt(w_j[j]))
    }
    
    s[i] = sample(1:J, 1, prob = q_j_h, replace=TRUE)
  }
  return(s)
}

gen_h = function(h_t, mu, sigma2_n, phi)
{
  h_t[1] = rnorm(1, mu, sqrt(sigma2_n/(1-phi**2)))
  for(i in 2:T)
  {
    h_t[i] = mu + (h_t[i-1] - mu) * phi + rnorm(1, 0, sqrt(sigma2_n))
  }
  return(h_t)
}

phi_log_pdf = function(h_t, mu, sigma2_n, phi)
{
  c = 1.0
  C = 0.5

  if(phi <= 0 || phi >= 1.0)
    return -Inf
  phi_est = 0.5 * log(1 - phi**2) + ((phi**2) * (h_t[1] -mu)**2)/(2 * sigma2_n)
  phi_est = phi_est + dnorm(phi, c, sqrt(C), log = TRUE)
  for(i in 2:T)
  {
    phi_est = phi_est + dnorm(h_t[i], mean = mu + phi * (h_t[i-1] - mu), sd = sqrt(sigma2_n), log = TRUE)
  }
  phi_est
}


mu_log_pdf = function(h_t, mu, sigma2_n, phi)
{
  mu_est = 0
  mu_est = mu_est + dnorm(h_t[1], mu, sqrt(sigma2_n/(1.0 - phi**2)), log = TRUE)
  for(i in 2:T)
  {
    mu_est = mu_est + dnorm(h_t[i], mu + phi * (h_t[i - 1] - mu), sqrt(sigma2_n), log = TRUE)
  }
  mu_est
}

v_log_pdf = function(h_t, mu, sigma2_n, phi)
{
  a = 1000
  v0 = 0.002

  if(sigma2_n <= 0.0)
    return -Inf
  v_est = dinvgamma(sigma2_n, a/2, a * v0/2.0, log = TRUE)
  v_est = v_est + dnorm(h_t[1], mu, sqrt(sigma2_n/(1 - phi**2)), log = TRUE)
  for(i in 2:T)
  {
    v_est = v_est + dnorm(h_t[i], mu + phi * (h_t[i-1] - mu), sqrt(sigma2_n), log = TRUE)
  }
  v_est
}

sample_mu = function(h_t, mu, sigma2_n, phi)
{
  sigma2_mu = sigma2_n/((T-1) * (1-phi)**2 + (1-phi**2))
  mu_cap = sigma2_mu * ((1-phi**2)/sigma2_n * h_t[1] + (1-phi)/sigma2_n * sum(h_t[2:T] - phi * h_t[1:(T-1)]))
  return(rnorm(1, mu_cap, sqrt(sigma2_mu)))
}

sample_phi = function(h_t, mu, sigma2_n, phi)
{
  phi_mcmc = phi
  for(i in 1:burnin_int)
  {
    prop = rtruncnorm(1, 0.9, 1.0, 0.7, 1.0)
    current_prob = phi_log_pdf(h_t, mu, sigma2_n, phi_mcmc)
    prop_prob = phi_log_pdf(h_t, mu, sigma2_n, prop)
    if(!is.finite(prop_prob))
    {
      next 
    }
    if(prop_prob > current_prob || (runif(1) < exp(prop_prob - current_prob)))
    {
      phi_mcmc = prop
    }
  }
  phi_mcmc
}

sample_phi = function(h_t, mu, sigma2_n, phi)
{
  a = 1.0
  b = 0.1
  
  b_2 = 1/(1/b**2 + 1/sigma2_n * sum(h_t[1:(T-1)]**2))
  a_ = b_2 * (a/b**2 + 1/sigma2_n * sum(h_t[1:(T-1)] * (h_t[2:T] - mu)))
  rtruncnorm(1, 0, 1, a_, sqrt(b_2))
}

sample_sigma2_n = function(h_t, mu, sigma2_n, phi)
{
  a = 1000
  v0 = 0.002
  sigma2_mcmc = sigma2_n
  for(i in 1:burnin_int)
  {
    prop = rinvgamma(1, a/2, a * v0/2.0)
    current_prob = v_log_pdf(h_t, mu, sigma2_mcmc, phi)
    prop_prob = v_log_pdf(h_t, mu, prop, phi)

    if(!is.finite(prop_prob))
     {
        next 
     }
    if(prop_prob > current_prob || (runif(1) < exp(prop_prob - current_prob)))
    {
      sigma2_mcmc = prop
    }
  }
  sigma2_mcmc
}

sample_s = function(y_star_t, h_t)
{
  s = array(0, T)
  q_j_h = array(0, J)
  for(i in 1:T)
  {
    for(j in 1:J)
    {
      q_j_h[j] = q_j[j] * dnorm(y_star_t[i], h_t[i] + b_j[j] - 1.2704, sqrt(w_j[j]))
    }
    s[i] = sample(1:J, 1, prob = q_j_h, replace=TRUE)
  }
  return(s)
}

sample_h_t <- function(T, y_star_t, s_t, mu, sigma2_n, phi)
{ 
  x_t_t = array(0, T + 1)
  P_t_t = array(0, T + 1)
  x_t_t_1 = array(0, T)
  P_t_t_1 = array(0, T)

  x_t_t[1] = mu
  P_t_t[1] = sigma2_n/(1 - phi**2)
  for(i in 1:T)
  {
    x_t_t_1[i] = mu + phi * (x_t_t[i] - mu)
    P_t_t_1[i] = phi**2 * P_t_t[i] + sigma2_n
    
    u = s_t[i] 
    K = P_t_t_1[i]/(P_t_t_1[i] + w_j[u])
    x_t_t[i + 1] = x_t_t_1[i] + K * (y_star_t[i] - (b_j[u] + x_t_t_1[i]))
    P_t_t[i + 1] = (1 - K**2) * P_t_t_1[i] + K**2 * w_j[u]
  }
  h_t = array(0, T + 1)
  h_t[T + 1] <- rnorm(1, x_t_t[T + 1], sqrt(P_t_t[T + 1]))
  for(t in T:1)
  {
    
    H <- 1.0 / ( 1.0/(sigma2_n) + 1.0/P_t_t[t] )
    h <- H * ( h_t[t+1]/(sigma2_n) + x_t_t[t]/P_t_t[t] )
    h_t[t] <- rnorm(1,h,sqrt(H))    
  }
  list(h_t, x_t_t, P_t_t)
}

sweep = function(T, y_star_t, h_t, mu, sigma2_n, phi, iters)
{
  phi_t = array(0, iters)
  mu_t = array(0, iters)
  sigma2_n_t = array(0, iters)
  
  phi_t[1] = phi
  mu_t[1] = mu
  sigma2_n_t[1] = sigma2_n
  for(i in 2:(iters + 1))
  {
    phi_t[i] = sample_phi(h_t, mu_t[i-1], sigma2_n_t[i-1], phi_t[i-1])
    mu_t[i] = mu #sample_mu(h_t, mu_t[i-1], sigma2_n_t[i-1], phi_t[i-1])
    sigma2_n_t[i] = sample_sigma2_n(h_t, mu_t[i-1], sigma2_n_t[i-1], phi_t[i-1])
    s_t = sample_s(y_star_t, h_t)
    out = sample_h_t(T, y_star_t, s_t,  mu_t[i], sigma2_n_t[i], phi_t[i])
    h_t = out[[1]]
    x_t_t = out[[2]]
    P_t_t = out[[3]]
    cat(i, "\r")
  }

  list("h_t"<-h_t, "mu"<-mu_t, "sigma2_n"<-sigma2_n_t,"phi"<-phi_t,"x_t_t"<-x_t_t, "P_t_t"<-P_t_t)
}

sigma2_n = sigma2_0
phi = phi_0
mu = mu_0

h_t = gen_h(h_t, mu_correct, sigma2_n, phi)
iters = 1000
out = sweep(T, y_star_t, h_t, mu_correct, sigma2_n, phi, iters)

par(mfrow=c(2,1))
plot(h_t_correct, type='l')
plot(out[[1]], type='l')

par(mfrow=c(3,1))
plot(out[[2]], type='l')
plot(out[[3]], type='l')
plot(out[[4]], type='l')

mean(out[[2]][500:iters])
mean(out[[3]][500:iters])
mean(out[[4]][500:iters])

par(mfrow=c(1,1))
plot(out[[5]], type='l',col='red')
lines(out[[6]] + out[[5]],lty=2,lwd=2)
lines(-out[[6]] + out[[5]],lty=2,lwd=2)

par(mfrow=c(2,1))
s_t = sample_s(y_star_t, out[[1]])
h_t_predicted = sample_h_t(T, y_star_t, s_t, mu_correct, sigma2_correct, 0.99)[[1]]
plot(h_t_correct, type='l')
plot(h_t_predicted, type='l')

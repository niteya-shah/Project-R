library('gmgm')
library('invgamma')
library('truncnorm')

T = 1000

phi_0 = 0.95
sigma2_0= 0.002
mu_0 = 0

c = 0.001

q_j = c(0.00730, 0.10556, 0.00002, 0.04395, 0.34001, 0.24566, 0.25750)
b_j = c(-10.12999, -3.397281, -8.56686, 2.77786, 0.61942, 1.79518, -1.08819)
w_j = c(5.79596, 2.61369, 5.17950, 0.16735, 0.64009, 0.34023, 1.26261)
J = length(q_j)

burnin_int = 10

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

y_star_t = log(y_t**2 + c)
plot(y_star_t, type='l')

h_t = array(0.0, T)
s_t = sample(1:J, T, prob = q_j, replace=TRUE)

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

gen_h = function(h_t, mu, sigma2_n, phi)
{
  h_t[1] = rnorm(1, mu, sqrt(sigma2_n/(1-phi**2)))
  for(i in 2:T)
  {
    h_t[i] = mu + (h_t[i-1] - mu) * phi + rnorm(1, 0, sqrt(sigma2_n))
  }
  return(h_t)
}

prob_phi_prior = function(phi)
{
  phi_star = (phi + 1)/2
  phi_1 = 20
  phi_2 = 1.5
  
  return(dbeta(phi_star, phi_1, phi_2, log=TRUE))
}

g_phi = function(h_t, mu, sigma2_n, phi)
{
  prob_phi_prior(phi) - ((h_t[1] - mu)^2 * (1 - phi**2))/(2 * sigma2_n) + 0.5 * log(1 - phi**2)
}

sample_phi = function(h_t, mu, sigma2_n, phi)
{
  temp = (h_t[1:(T-1)] - mu)
  phi_cap = sum((h_t[2:T] - mu) * temp)/sum(temp)
  V_phi = sigma2_n/sum((h_t[1:T - 1] - mu)^2)
  
  phi_star = rnorm(1, phi_cap, V_phi)
  
  if(phi_star > 0 && phi_star < 1)
  {
    prop_prob = exp(g_phi(phi_star) - g(phi))
    acc_prob  = runif(1)
    if(prob_prob > acc_prob)
    {
      return(phi_star)
    }
  }
  return(phi)
}

sample_sigma2_n = function(h_t, mu, sigma2_n, phi)
{
  sigma_r = 5
  S_phi = 0.01 * sigma_r
  temp = (h_t[1] - mu)**2 * (1 - phi**2) + sum(((h_t[2:T] - mu) - phi * (h_t[1:(T-1)] - mu))**2)
  return(1/rgamma(1, (T + sigma_r)/2, (S_phi + temp)/2))
}

sample_mu = function(h_t, mu, sigma2_n, phi)
{
  sigma2_mu = sigma2_n/((T-1) * (1-phi)**2 + (1-phi**2))
  mu_cap = sigma2_mu * ((1-phi**2)/sigma2_n + (1-phi)/sigma2_n * sum(h_t[2:T] - phi * h_t[1:(T-1)]))
  return(rnorm(1, mu_cap, sqrt(sigma2_mu)))
}

sample_h_t <- function(T, y_star_t, s_t, mu, sigma2_n, phi)
{ 
  c_t = b_j[s_t]
  d = mu * (1 - phi)

  gamma_t = array(0, T)
  P_t = array(0, T)
  v_t = array(0, T)
  F_t = array(0, T)
  K_t = array(0, T)
  L_t = array(0, T)

  P_t_0 = sigma2_n/(1 - phi**2)
  gamma_t_0 = rnorm(1, mu, sqrt(P_t_0))  
  v_t[1] = y_star_t[1] - gamma_t_0 - c_t[1]

  F_t[1] = P_t_0 + sigma2_n
  K_t[1] = P_t_0/F_t[1]
  L_t[1] = phi - K_t[1]
  gamma_t[1] = d + phi * gamma_t_0 + K_t[1] * v_t[1]
  P_t[1] = phi * P_t_0 * L_t[1] + sigma2_n

  for(i in 2:T)
  {
    v_t[i] = y_star_t[i] - gamma_t[i] - c_t[i]
    F_t[i] = P_t[i] + sigma2_n
    K_t[i] = P_t[i]/F_t[i]
    L_t[i] = phi - K_t[i]
    gamma_t[i] = d + phi * gamma_t[i-1] + K_t[i] * v_t[i]
    P_t[i] = phi * P_t[i-1] * L_t[1] + sigma2_n
  }
  h_t = array(0, T)
  h_t[T] <- rnorm(1,gamma_t[T],sqrt(P_t[T]))
  for(t in (T-1):1)
  {
    H <- 1.0 / ( 1.0/(sigma2_n) + 1.0/P_t[t] )
    h <- H * ( h_t[t+1]/(sigma2_n) + gamma_t[t]/P_t[t] )
    h_t[t] <- rnorm(1,h,sqrt(H))
  }
  list(h_t, gamma_t, P_t)
}

sigma2_n = sigma2_0
phi = phi_0
mu = mu_0

h_t = gen_h(h_t, mu, sigma2_n, phi)
plot(h_t, type='l')
for(i in 1:10)
{
  phi = sample_phi(h_t, mu, sigma2_n, phi)
  mu = sample_mu(h_t, mu, sigma2_n, phi)
  sigma2_n = sample_sigma2_n(h_t, mu, sigma2_n, phi)
  s_t = sample_s(y_star_t, h_t)
  out = sample_h_t(T, y_star_t, s_t,  mu, sigma2_n, phi)
  h_t = out[[1]]
  gamma_t = out[[2]]
  P_t = out[[3]]
  cat(phi, mu, sigma2_n, "\n")
}
plot(h_t,type='l')

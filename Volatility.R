library('gmgm')
library('invgamma')
library('truncnorm')

T = 1000
c = 1
C = 0.01

phi_0 = 0.95
v_0 = 0.002
a = 1000

burnin_int = 10

v = v_0
phi = phi_0

sample_mixture = function(y, z)
{
  q_j = c(0.0073, 0.0001, 0.10549, 0.2575, 0.3400, 0.2457, 0.0440)
  b_j = c(-5.7002, -4.9186, -2.6216, -1.1793, -0.3255, 0.2624, 0.7537)
  w_j = c(1.4490, 1.2949, 0.6534, 0.3157, 0.1600, 0.0851, 0.0418)
  chosen = sample(x=1:7, 1, prob=q_j, replace=TRUE)
  
  c(rnorm(n=1, mean=b_j[chosen] + y - z, sd=w_j[chosen]), b_j[chosen] + y - z, w_j[chosen])
}

data = read.csv('DEXUSUK_2.csv')
length(data[,1])
stock_price = ts(data[1:(T+1),2])
plot(stock_price)

r_t = stock_price[1:T]/stock_price[2:(T + 1)] - 1
plot(abs(r_t * 100),type='l')

y_t = log(r_t**2 + 1e-8)/2
plot(y_t, type='l')

mu_0 = mean(y_t) - mean(q_j * b_j)
mu = mu_0


z_t = array(NA, dim=T)
z_t[1] = rnorm(1, mu, v/(1-phi)**2)
for(t in 2:T) z_t[t] = phi * (z_t[t-1] - mu) + mu + rnorm(1,mean=0,sd=sqrt(v))

plot(z_t, type='l')

gamma_t = array(NA, T)
beta_t = array(NA, T)
w_t = array(NA, T)

for(i in 1:T)
{
  temp = sample_mixture(y_t[i], z_t[i])  
  gamma_t[i] = temp[1]
  beta_t[i] = temp[2]
  w_t[i] = temp[3]
}

plot(gamma_t, type='l')

phi_log_pdf = function(phi_v, z_t, mu, v)
{
  if(phi_v <= 0 || phi_v >= 1.0)
    return -Inf
  phi_est = 0.5 * log(1 - phi_v**2) + ((phi_v**2) * (z_t[1] -mu)**2)/(2 * v)
  phi_est = phi_est + log(dnorm(phi_v, c, sqrt(C)))
  for(i in 2:T)
  {
    phi_est = phi_est + log(dnorm(z_t[i], mean = mu + phi_v * (z_t[i-1] - mu), sd = sqrt(v)))
  }
  phi_est
}

mu_log_pdf = function(mu_v, z_t, a, phi, v)
{
  mu_est = mu_0
  mu_est = mu_est + log(dnorm(z_t[1], mu_v, sqrt(v/(1.0 - phi**2))))
  for(i in 2:T)
  {
    mu_est = mu_est + log(dnorm(z_t[i], mu_v + phi * (z_t[i - 1] - mu_v), sqrt(v)))
  }
  mu_est
}

v_log_pdf = function(v_v, z_t, phi, mu)
{
  if(v <= 0.0)
    return -Inf
  v_est = log(dinvgamma(v_v, a/2, a * v_v/2.0))
  v_est = v_est + log(dnorm(z_t[1], mu, v_v/(1 - phi**2)))
  for(i in 2:T)
  {
    v_est = v_est + log(dnorm(z_t[i], mu + phi * (z_t[i-1] - mu), sqrt(v_v)))
  }
  v_est
}

phi_mcmc = phi
for(i in 1:burnin_int)
{
  prop = rtruncnorm(1, 0, 1.0, 0, 1.0)
  current_prob = phi_log_pdf(phi_mcmc, z_t, mu, v)
  prop_prob = phi_log_pdf(prop, z_t, mu, v)
  if(is.finite(prop_prob) && current_prob > prop_prob || (runif(1) < exp(prop_prob - current_prob)))
  {
    phi_mcmc = prop
  }
}

mu_mcmc = mu
for(i in 1:burnin_int)
{
  prop = rnorm(1, 0, 1.0)
  current_prob = mu_log_pdf(mu_mcmc, z_t, a, phi, v)
  prop_prob = mu_log_pdf(prop, z_t, a, phi, v)
  if(is.finite(prop_prob) && current_prob > prop_prob || (runif(1) < exp(prop_prob - current_prob)))
  {
    mu_mcmc = prop
  }
}

v_mcmc = v
for(i in 1:burnin_int)
{
  prop = rinvgamma(1, a/2, a * v/2.0)
  current_prob = v_log_pdf(v_mcmc, z_t, phi, mu)
  prop_prob = v_log_pdf(prop, z_t, phi, mu)
  if(is.finite(prop_prob) && current_prob > prop_prob || (runif(1) < exp(prop_prob - current_prob)))
  {
    v_mcmc = prop
  }
}

phi = phi_mcmc
mu = mu_mcmc
v = v_mcmc

m_0 = mu
C_0 = v/(1-phi**2)


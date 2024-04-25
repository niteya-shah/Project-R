library('gmgm')

T = 1000
c = 1
C = 0.01

phi_0 = 0.95
v_0 = 0.002
a = 1000

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

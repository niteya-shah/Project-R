import numpy as np
import pandas as pd
import scipy
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture
import emcee

df = pd.read_csv('DEXUSUK_1.csv', index_col=0,on_bad_lines='skip')

T = 1000
c = 1
C = 0.01

phi_0 = 0.95
v_0 = 0.002
a = 1000
mu_0 = 0.0

mu = mu_0
v = v_0
phi = phi_0
walkers = 2

def sample_mixture(y, z):
    q_j = np.array([0.0073, 0.0000, 0.1055, 0.2575, 0.3400, 0.2457, 0.0440])
    b_j = np.array([-5.7002, -4.9186, -2.6216, -1.1793, -0.3255, 0.2624, 0.7537]).reshape(-1, 1) + y - z
    w_j = np.array([1.4490, 1.2949, 0.6534, 0.3157, 0.1600, 0.0851, 0.0418]).reshape(-1, 1, 1)

    v_t = GaussianMixture(n_components=7)
    v_t.weights_ = q_j
    v_t.means_ = b_j
    v_t.covariances_ = w_j**2

    return v_t.sample(1)[0].squeeze()
df
P_T = pd.read_csv('DEXUSUK_2.csv',index_col=0).iloc[:T + 1].DEXUSUK.to_numpy()
r_t = P_T[1:]/P_T[:-1] - 1
plt.plot(np.abs(r_t * 100))

y_t = np.log(r_t**2 + 1e-5)/2
plt.plot(y_t,'--')

z_t = np.empty_like(y_t)
z_t[0] = np.random.normal(mu, v/(1-phi)**2)
for i in range(1, T):
    z_t[i] = mu + (z_t[i - 1] - mu) + np.random.normal(0, v)

gamma_t = np.empty_like(z_t)
for i in range(T):
    gamma_t[i] = sample_mixture(y_t[i], z_t[i])
plt.plot(gamma_t)

## PHI
phi_v = phi
scipy.stats.norm.pdf(0)
def phi_log_pdf(phi_v, z_t, mu, v):
    if phi_v <= 0 or phi_v >= 1.0:
        return -np.inf, None
    phi_est = 0.5 * np.log(1 - phi_v**2) + ((phi_v**2) * (z_t[0] - mu)**2)/(2 * v)
    phi_est += scipy.stats.norm.logpdf(phi_v, loc=c, scale=C)
    for i in range(1, T):
        phi_est += scipy.stats.norm.logpdf(z_t[i], loc=mu + phi_v * (z_t[i-1] - mu), scale=v)
    return phi_est, None

def mu_log_pdf(mu_v, z_t, a, phi, v):
    mu_est = 0.0 # TODO
    mu_est += scipy.stats.norm.logpdf(z_t[0], loc=mu_v, scale=v/(1.0-phi**2))
    for i in range(1, T):
        mu_est += scipy.stats.norm.logpdf(z_t[i], loc=mu_v + phi * (z_t[i-1] - mu_v), scale=v)
    return mu_est, None

def v_log_pdf(v_v, z_t, phi, mu):
    if v_v <= 0.0:
        return -np.inf, None
    v_est = scipy.stats.invgamma.logpdf(v_v, a/2, scale = a * v_v/2.0)
    v_est += scipy.stats.norm.logpdf(z_t[0], loc=mu, scale=v_v/(1.0-phi**2))
    for i in range(1, T):
        v_est += scipy.stats.norm.logpdf(z_t[i], loc=mu + phi * (z_t[i-1] - mu), scale=v_v)
    return v_est, None


phi_sampler = emcee.EnsembleSampler(walkers, 1, phi_log_pdf, args=[z_t, mu, v])

phi_state = phi_sampler.run_mcmc(np.random.normal(phi, 0.01,size=walkers)[:, None], 50,skip_initial_state_check=True)
phi_sampler.reset()
out_phi = phi_sampler.run_mcmc(phi_state, 1)
out_phi.coords

mu_sampler = emcee.EnsembleSampler(walkers, 1, mu_log_pdf, args=[z_t, a, phi, v])
mu_state = mu_sampler.run_mcmc(np.random.normal(mu, 0.01,size=walkers)[:, None], 50,skip_initial_state_check=True)
mu_sampler.reset()
out_mu = mu_sampler.run_mcmc(mu_state, 1)
out_mu.coords

v_sampler = emcee.EnsembleSampler(walkers, 1, v_log_pdf, args=[z_t, phi, mu])
v_state = v_sampler.run_mcmc(np.random.normal(v, 0.01,size=walkers)[:, None], 50,skip_initial_state_check=True)
mu_sampler.reset()
out_v = mu_sampler.run_mcmc(v_state, 1)
out_v.coords

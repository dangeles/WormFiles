import collections
import warnings

# Our numerical workhorses
import numpy as np
import pandas as pd
import scipy.optimize
import scipy.stats as st

# Numba to make things faster
import numba

# The MCMC Hammer
import emcee

# Numerical differentiation package
import numdifftools as ndt

# Import plotting tools
import matplotlib.pyplot as plt
import seaborn as sns
import corner


# JB's favorite Seaborn settings for notebooks
rc = {'lines.linewidth': 2,
      'axes.labelsize': 18,
      'axes.titlesize': 18,
      'axes.facecolor': 'DFDFE5'}
sns.set_context('notebook', rc=rc)
sns.set_style('darkgrid', rc=rc)

# Suppress future warnings
warnings.simplefilter(action='ignore', category=FutureWarning)


input_path = '../input/rnai_screen_results/'
output_path = '../output/rnai_screen_results/'

df = pd.read_csv(input_path + 'rnai_heat_shock_data.txt', sep='\t')
names = pd.read_csv(input_path + 'rnai_genes_dict.csv')
names.head()

# rename the columns to something handy
df.columns = ['rnai', 'alive', 'dead', 'date']

# make all codes upper or lower, not both
# first make sure each column is a str
names.code = names.code.apply(str)
df.rnai = df.rnai.apply(str)
# now apply lower
names.code = names.code.apply(str.lower)
df.rnai = df.rnai.apply(str.lower)


# extract the names that have been assayed so far
def translate(x):
    """a function to go between rnai code and gene (human-readable) name."""
    return names[names.code == x].gene_name.values[0]

df['gene'] = df.rnai.apply(translate)

df.sort_values('gene', inplace=True)

# calculate fraction dead
df['fraction_dead'] = df.dead/(df.alive + df.dead)


# standardize dead values to the mean and stdev of the gfp data
gfp_fr_dead_mu = df.fraction_dead.mean()
gfp_fr_dead_sig = df.fraction_dead.std()
total = (df.dead+df.alive)

df['z_dead'] = (df.dead - total*gfp_fr_dead_mu)/(gfp_fr_dead_sig*total)

# plot in either a swarmplot or boxplot
plot = sns.swarmplot(x='gene', y='fraction_dead', data=df)
plt.xticks(rotation=30)
plt.title('Heat Survival After 24 hours')
plt.savefig(output_path + 'swarmplot_gene_heat_shock_assays.pdf')

plot = sns.swarmplot(x='gene', y='z_dead', data=df)
plt.xticks(rotation=30)
plt.title('Heat Survival After 24 hours')
plt.show()

plot = sns.boxplot(x='gene', y='z_dead', data=df)
plt.xticks(rotation=30)
plt.title('Heat Survival After 24 hours')
plt.ylim(-3, 3)
plt.show()

sns.boxplot(x='gene', y='fraction_dead', data=df)
plt.xticks(rotation=30)
plt.title('Heat Survival After 24 hours')
plt.savefig(output_path + 'boxplot_gene_heat_shock_assays.pdf')
plt.show()

sns.boxplot(x='date', y='fraction_dead', data=df[df.gene == 'gfp'])
plt.xticks(rotation=30)
plt.title('Day-to-Day Variation in Heat Survival After 24 hours')
plt.savefig(output_path + 'boxplot_gene_heat_shock_controls_by_date.pdf')
plt.show()


# first, identify outliers in the data
# in theory, it should be binomial distributed, and since
# we are using a normal approximation
def log_posterior_good_bad(p, x):
    """The log posterior for good/bad data model for repeated measurements."""
    # Pull out parameters
    mu, sigma, sigma_bad = p[:3]
    g = p[3:]

    if type(sigma) is list:
        raise ValueError('sigma is list')

    if type(mu) is not np.float64:
        raise ValueError('mu is not float')
    if type(sigma) is not np.float64:
        raise ValueError('mu is not float')
    if type(sigma_bad) is not np.float64:
        raise ValueError('mu is not float')
    if type(g) is not list:
        raise ValueError('g is not list')
    # Check to make sure the prior conditions are ok
    # if any(i < 0.0 for i in g):
    #     return -np.inf
    # if any(i > 1.0 for i in g):
    #     return -np.inf
    #
    # if sigma >= 0:
    #     return -np.inf
    # if sigma_bad < sigma:
    #     return -np.inf
    # log prior
    log_prior = -np.log(sigma) - np.log(sigma_bad)

    # Add in likelihood
    log_like_good = np.log(g / sigma) - ((x - mu) / sigma)**2 / 2.0

    log_like_bad = np.log((1.0 - g) / sigma_bad) \
        - ((x - mu) / sigma_bad)**2 / 2.0

    log_like = np.logaddexp(log_like_good, log_like_bad).sum()

    # Return the whole posterior
    return log_prior + log_like

# Set up MCMC parameters
n_dim = 3 + len(df[df.gene == 'gfp'])  # number of parameters in the model
n_walkers = 200                  # number of MCMC walkers
n_burn = 2000                    # "burn-in" period to let chains stabilize
n_steps = 5000                 # number of MCMC steps to take after burn-in

# Seed random number generator for reproducibility
np.random.seed(42)

# Generate random starting points for walkers.
# p0[i,j] is the starting point for walk i along variable j.
p0 = np.empty((n_walkers, n_dim))
p0[:, 0] = np.random.uniform(10, 30, n_walkers)                # mu
p0[:, 1] = np.random.exponential(5.0, n_walkers)               # sigma
p0[:, 2] = np.random.exponential(20.0, n_walkers)              # sigma_bad
p0[:, 3:] = np.random.uniform(0.0, 1.0, (n_walkers, n_dim-3))  # g_i

# Set up the EnsembleSampler instance
sampler = emcee.EnsembleSampler(n_walkers, n_dim, log_posterior_good_bad,
                                args=(df[df.gene == 'gfp'].z_dead,),
                                threads=6)

# Do the burn-in
pos, prob, state = sampler.run_mcmc(p0, n_burn, storechain=False)

# Reset sampler and run from the burn-in state we got to
_ = sampler.run_mcmc(pos, n_steps)


# Get most probable parameter value
max_ind = np.argmax(sampler.flatlnprobability)
mean_goodbad = sampler.flatchain[max_ind, 0]
sigma = sampler.flatchain[max_ind, 1]
# Get the error bar
sem_goodbad = sampler.flatchain[:, 0].std()

# Report results
print("""
Good/bad data model: {0:.2f} ± {1:.2f} sec/10 min
""".format(mean_goodbad, sem_goodbad))


sns.stripplot(y='z_dead', data=df[df.gene == 'gfp'], jitter=True)
plt.plot(plt.gca().get_xlim(), [mean_goodbad, mean_goodbad], '-',
         color=sns.color_palette()[1], label='Good/bad')
plt.ylabel('mean activity (sec/10 min)')
plt.legend()
plt.show()


# Compute mean goodness of data
g = sampler.flatchain[:, 3:].mean(axis=0)

# Identify outliers
outliers = (g < g.mean() - 1.7*g.std())

# Make strip plot with outliers in red
sns.stripplot(y='z_dead', data=df[df.gene == 'gfp'][~outliers],
              jitter=True)
sns.stripplot(y='z_dead', data=df[df.gene == 'gfp'][outliers],
              jitter=True, color=sns.color_palette()[2])
plt.ylabel('mean activity (sec/10 min)')

plt.show()


sigma
np.sqrt(sigma)

corner.corner(sampler.flatchain[:, :3],
              labels=[r'$\mu$', r'$\sigma$', r'$\sigma_\mathrm{bad}$'],
              bins=100)
plt.show()


# get the number of groups
n_genes_tested = df.gene.unique().shape[0]

# number of tests per group
n = np.array([len(x) for g in df.gene.unique() for x in df[df.gene == g]])

# get the indices of each entry
n = np.array([len(df[df.gene == g]) for g in df.gene.unique()])
n = n.astype(int)

inds = np.concatenate(((0,), n.cumsum()))

z_dead_vals = df.z_dead.values

df.shape
len(inds)

def log_post_hierarchical(params, data, inds, total=df.shape[0]):

    if len(data) < 2:
        raise ValueError('Too few datapoints to run simulation')

    L = len(inds)-1
    mu = params[0: L]
    sigma = params[L: 2*L]
    sigma_bad = params[2*L: 3*L]
    G = params[3*L: 4*L]
    g = params[4*L:4*L + total]
    M = params[-3]
    S = params[-2]
    S_bad = params[-1]

    # cropping
    if any(i <= 0 for i in sigma):
        return -np.inf

    if any(i <= 0 for i in sigma_bad):
        return -np.inf
    if (sigma_bad < sigma).any():
        return -np.inf

    if S <= 0:
        return -np.inf
    if (S_bad <= 0):
        return -np.inf
    if (S_bad < S):
        return - np.inf

    if any(i < 0 for i in g):
        return -np.inf
    if any(i < 0 for i in G):
        return -np.inf
    if any(i > 1 for i in g):
        return -np.inf
    if any(i > 1 for i in G):
        return -np.inf

    # log likelihood calculation
    lp = 0
    for i, m in enumerate(mu):
        # extract data
        # fit each parameter
        # data for log_posterior_good_bad
        x = data[inds[i], inds[i+1]]
        gs = g[inds[i]: inds[i+1]]

        p = [m, sigma[i], sigma_bad[i], gs]
        P = [M, S, S_bad, G[i]]
        # print(p[0:3])
        print(P[0:3])
        print(type(m))
        print(type(M))
        lp += log_posterior_good_bad(p, x)  # outlier detection for bad data
        lp += log_posterior_good_bad(P, m)  # outlier detection, effects!

    return lp


def mcmc(data, inds, n_burn=5000, n_steps=5000, a=2, total=df.shape[0]):
    """
    A parser that sets up and executes an emcee-based MCMC for a
    hierarchical binomial beta distribution.
    Takes in two lists of data (natt, nc) and returns an
    EnsembleSampler object.
    """
    n_dim = 4*len(inds) - 4 + total + 3    # no. of dimensions
    n_walkers = n_dim*50   # number of MCMC walkers

    # p0[i,j] is the starting point for walk i along variable j.
    p0 = np.empty((n_walkers, n_dim))

    # initialize parameters
    p0 = np.empty((n_walkers, n_dim))

    L = len(inds) - 1
    for i in range(0, len(inds)):
        p0[:, i] = np.random.uniform(-5, 5, n_walkers)             # mu
        p0[:, i + L] = np.random.exponential(5.0, n_walkers)       # sigma
        p0[:, i + 2*L] = np.random.exponential(5.0, n_walkers)     # sigma_bad
        p0[:, i + 3*L] = np.random.uniform(0.0, 1.0, n_walkers)    # G_i

    for i in range(0, total):
        p0[:, i + 4*L] = np.random.uniform(0.0, 1.0, n_walkers)    # g_i

    p0[:, -3] = np.random.uniform(-5, 5, n_walkers)                # M
    p0[:, -2] = np.random.exponential(5.0, n_walkers)              # S
    p0[:, -1] = np.random.exponential(5.0, n_walkers)              # S_bad

    # set up the sampler
    sampler = emcee.EnsembleSampler(n_walkers, n_dim, log_post_hierarchical,
                                    args=(data, inds,), threads=6, a=a)

    # Do burn-in
    pos, prob, state = sampler.run_mcmc(p0, n_burn, storechain=False)

    # Sample again, starting from end burn-in state
    _ = sampler.run_mcmc(pos, n_steps)

    return sampler


sampler = mcmc(z_dead_vals, inds, n_burn=5000, n_steps=5000, a=2,
               total=df.shape[0])

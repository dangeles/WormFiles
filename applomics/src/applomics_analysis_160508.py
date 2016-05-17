"""
A script to analyze applomics data.

author: dangeles@caltech.edu
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

df1 = pd.read_csv('../input/apple_inoculation_expt_160508.csv')

# count colonies adjusted for dilution and volume plated
df1['cfu_per_apple'] = df1.colonies * \
                       df1.dilution_factor/df1.volume_plated*df1.apple_mass + 5
df1.dropna(inplace=True)

# plot cfu vs inoculation factor
fig, ax = plt.subplots()
df1[df1.worms == 0].plot('inculation_factor', 'cfu_per_apple', 'scatter',
                         logx=True, logy=True)
df1[df1.worms == 1].plot('inculation_factor', 'cfu_per_apple', 'scatter',
                         logx=True, logy=True)
plt.show()

# since 10**-8 seems like an odd value, so remove it from this experiment.
df2 = df1[df1.inculation_factor > 10**-8].copy()
mean_growth_7 = df1[(df1.worms == 0) &
                    (df1.inculation_factor == 10**-7)].cfu_per_apple.mean()
mean_growth_6 = df1[(df1.worms == 0) &
                    (df1.inculation_factor == 10**-6)].cfu_per_apple.mean()

divider = np.repeat([mean_growth_6, mean_growth_7], [6, 6])
df2['fold_change'] = df2.cfu_per_apple/divider

plt.plot(df2[df2.worms == 0].inculation_factor,
         df2[df2.worms == 0].fold_change, 'bo', ms=10, alpha=0.65,
         label='No Worms/Mean(No Worms)')
plt.plot(df2[df2.worms == 1].inculation_factor,
         df2[df2.worms == 1].fold_change, 'ro', ms=10, alpha=0.65,
         label='Worms/Mean(No Worms)')
plt.xlim(5*10**-8, 2*10**-6)
plt.xscale('log')
plt.ylabel('Fold Change (worms/no worms)')
plt.xlabel('Inoculation Factor (dilution from sat. soln)')
plt.title('Effect of Worms on Bacteria')
plt.legend()
plt.savefig('../output/Fold_Change_Applomics_160508_Expt1.pdf')
plt.show()

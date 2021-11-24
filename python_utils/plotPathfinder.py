#!/usr/bin/env python

# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'
# %%
import pandas as pd
import matplotlib.pyplot as plt


# %%
xtbReactionPath = pd.read_csv('path.csv')

X = xtbReactionPath['opt_step']
y = xtbReactionPath['rel_energy']


# %%
fig, ax = plt.subplots(nrows=1, ncols=1)
ax.scatter(X,y, s=25, c='black')
ax.set_xlabel("Reaction Coordinate")
ax.set_ylabel("Relative Energy (kcal/mol)")

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)


# %%
fig.savefig('reactionpath.png')



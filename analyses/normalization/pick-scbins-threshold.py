import pickle

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from ld import FRAGMENTS

survived = []

with open(FRAGMENTS, 'rb') as stream:
    scbins = pickle.load(stream)

for assembly, df in scbins.items():
    df = df.set_index('source').T

    for threshold in [0, 1, 2, 4, 8, 16, 32, 64]:
        record = (df > threshold).mean(axis=0)
        record.name = f"N={threshold}"
        survived.append(record)

survived = pd.concat(survived, axis=1)
survived = survived.melt(ignore_index=False, var_name="N", value_name="Survived (%)").reset_index()

survived["N"] = survived["N"].str.extract(r"N=(\d+)").astype(int)

fig, ax = plt.subplots()
sns.lineplot(data=survived, x="N", y="Survived (%)", hue="source", legend=False, ax=ax)
sns.despine(fig=fig, ax=ax)
ax.axvline(5, color='red', linestyle='--')
fig.show()
plt.close(fig)

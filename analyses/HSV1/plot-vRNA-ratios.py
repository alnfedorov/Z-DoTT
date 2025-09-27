from urllib.parse import unquote

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from BCBio import GFF
from biobit.toolkit import nfcore

import ld
from assemblies import HSV1
from stories.nextflow.series.internal import B319096

# Fetch all vRNA IDs
vRNAs = set()

contigs = list(GFF.parse(HSV1.gff3))
assert len(contigs) == 1
contig = contigs[0]

types = set()
for record in contig.features:
    if record.type == 'gene':
        for RNA in record.sub_features:
            vRNAs.add(RNA.id)

# Fetch all vRNA TPMs
tables = []
for exp in B319096.experiments:
    # Select only FLAG RIP experiments
    if exp.sample.attributes['treatment'] == 'HSV-1':
        salmon = nfcore.rnaseq.extract.salmon(exp)
        salmon = pd.read_csv(salmon, sep='\t')

        salmon = salmon[~salmon['Name'].str.startswith('ENSMUST')].copy()
        salmon['Name'] = salmon['Name'].apply(unquote)

        # Classify vRNAs and collapse TPMs by category
        salmon['Category'] = salmon['Name'].apply(HSV1.utils.classify)
        salmon['Category'] = salmon['Category'].apply(lambda x: x if x == 'Uncharacterized' else "Canonical")
        salmon = salmon.groupby('Category', as_index=False)[['TPM']].sum()

        # Normalize by the total vTPM
        total = salmon['TPM'].sum()
        salmon['TPM'] = salmon['TPM'] / total * 100

        # Fill the meta
        selection = "RIP input"
        for ip in "FLAG RIP", "Z22 RIP":
            if ip in exp.library.selection:
                selection = ip

        salmon['Selection'] = selection
        salmon['Sample'] = exp.sample.ind

        tables.append(salmon)

# Concatenate all tables
df = pd.concat(tables)

# Use only FLAG RIP and matched input
allowed = df.loc[df['Selection'] == "FLAG RIP", 'Sample'].unique()
df = df.loc[df['Sample'].isin(allowed)].copy()

# Make a bar plot
palette = {"Canonical": "#e7704e", "Uncharacterized": "#9f9fa3"}
order = ['RIP input', 'FLAG RIP']

fig, ax = plt.subplots(figsize=(4.8 * 0.75, 6.4 * 0.75))
sns.barplot(
    data=df, x='Selection', y='TPM', hue='Category', order=order, palette=palette, lw=0.75, edgecolor='black',
    errorbar=("sd", 1), ax=ax, capsize=0.25
)
for container in ax.containers:
    for bar in container:
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_y() + bar.get_height(),
            f"{bar.get_height():.1f}%",
            ha='center',
            va='bottom',
            rotation=15
        )
sns.stripplot(
    data=df, x='Selection', y='TPM', hue='Category', order=order, palette=palette, ax=ax, size=6,
    edgecolor='black', linewidth=0.75, dodge=True, jitter=True, legend=False
)
ax.spines[['top', 'right']].set_visible(False)
ax.set_ylim(0, 100)
ax.set_yticks(range(0, 101, 25))
ax.set_ylabel("vRNA (%)")
ax.set_xlabel("")
ax.set_title("Estimated vRNA concentration")

saveto = ld.PLOTS / f"vRNA-fragments.svg"
saveto.parent.mkdir(parents=True, exist_ok=True)

fig.savefig(saveto, bbox_inches="tight", pad_inches=0, transparent=True)
# fig.show()
plt.close(fig)

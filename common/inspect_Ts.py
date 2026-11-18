"""
inspect_Ts.py  -  explore MOO vs CP couplings and the dimer geometry features.

Reads (from the run folder):
    Ts_compare.csv    from gen_comp_csv.py
    Ts_features.csv   from gen_features_csv.py

Run cell by cell in VS Code (# %% cells), or as a script.
"""

# %% imports
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


from scipy import stats                     # correlations, fits, distributions

# scikit-learn: only needed once you start the ML part
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_absolute_error, r2_score


RUN = 'Y6_20000ps_5A'
run_dir = os.path.expanduser(f'~/MSci/runs/{RUN}')


both  = pd.read_csv(os.path.join(run_dir, 'Ts_compare.csv'))
feats = pd.read_csv(os.path.join(run_dir, 'Ts_features.csv'))
print(f"{len(both)} dimers in compare, {len(feats)} in features")

summary = {}
for carrier, name in [('h', 'Holes'), ('e', 'Electrons')]:
    x = both[f'T{carrier}_CP_meV'].abs()
    y = both[f'T{carrier}_MOO_meV'].abs()

    # log correlation needs both values > 0
    m = (x > 0) & (y > 0)
    r_log = np.corrcoef(np.log10(x[m]), np.log10(y[m]))[0, 1]
    r_lin = np.corrcoef(x, y)[0, 1]

    summary[name] = {
        'N pairs': len(x),
        'N used (log)': int(m.sum()),
        'Pearson r (linear)': r_lin,
        'Pearson r (log10)': r_log,
        'MAE (meV)': (x - y).abs().mean(),
    }

# a bunch of stats


# %% quick look
print(both.columns.tolist())
both.describe().T



#%% Basic comparison between models
#list of possible columns
"""['init', 'fin', 'fin_lattice', 'dist',
            'Th_MOO_meV', 'Th_CP_meV', 'frac_diff_h',
            'Te_MOO_meV', 'Te_CP_meV', 'frac_diff_e',
            'time_ps']"""

min_dist = 1.
colour_feat = 'dist' #change as required
use_log = False #boolean, use log scale or not

fig, axes = plt.subplots(1, 2, figsize=(10, 4.6))
for ax, carrier, title, lab in zip(axes, ['h', 'e'], ['Holes', 'Electrons'], summary):
    x = both[f'T{carrier}_CP_meV'].abs()
    y = both[f'T{carrier}_MOO_meV'].abs()
    c = both[colour_feat] if colour_feat in both else None

    sc = ax.scatter(x, y, c=c, s=18, cmap='viridis')
    pos = pd.concat([x, y])
    pos = pos[pos > 0] #so we only keep non-zero couplings
    lo, hi = pos.min() / 2, pos.max() *1.1
    ax.plot([lo, hi], [lo, hi], 'k--', lw=0.8, label='y = x')


    scale = 'log' if use_log else 'linear'

    ax.set(xscale=scale, yscale=scale, xlim=(lo, hi), ylim=(lo, hi),
            xlabel='|J| counterpoise (meV)', ylabel='|J| MOO (meV)',
            title=f"{title}  (r_log = {summary[lab]['Pearson r (log10)']:.2f})")
    ax.legend(loc='upper left', frameon=False)
    if c is not None:
        fig.colorbar(sc, ax=ax, label='closest C–C distance (Å)')
fig.suptitle(RUN)
fig.tight_layout()
out_png = os.path.join(run_dir, 'Ts_compare.png')
fig.savefig(out_png, dpi=150)
print(f"Wrote {out_png}")

# %%

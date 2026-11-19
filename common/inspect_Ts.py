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

min_dist, max_dist = None, None    # e.g. 3.0, 4.0  (Å); None = no cut
colour_feat = 'dist'
colour_label = {'dist': 'closest C–C distance (Å)', 'time_ps': 'time (ps)'}.get(colour_feat, colour_feat)
use_log = False #boolean, use log scale or not
clip_colour = None   # e.g. (3.1, 3.7) to stretch the colour map over the bulk of points

# dimers identified by (init, fin, fin_lattice)
exclude = set()                    # dropped completely, e.g. {(101, 106, 23)}
flag    = {(116, 162, 0)}                   # drawn as red crosses, but are left out of r_log calc

#instead of listing dimers by hand, could also do so by creating a rule:
flag_rule = {
    'h': [lambda d: d['Th_CP_meV'].abs() > 150],   # one rule, still inside [ ]
    'e': [],                                       # no rules
}

ids = pd.Series(list(zip(both.init, both.fin, both.fin_lattice)), index=both.index)
#array storing our dimer identifying values

keep = ~ids.isin(exclude) #(~ means not, so keep gives a true/false column that False if we should exclude)
if min_dist is not None: keep = keep & both['dist'] >= min_dist
if max_dist is not None: keep = keep & both['dist'] <= max_dist
#our filtering stage
df = both[keep]
ids_keep = ids[keep] #an array of the kept ids

# one shared colour scale, built only from the kept points
cvals = df[colour_feat]
lo_pct, hi_pct = 5, 95          # ignore the lowest and highest 5% when setting colours
vmin, vmax = np.percentile(cvals, [lo_pct, hi_pct])
norm = plt.Normalize(vmin, vmax)


fig, axes = plt.subplots(1, 2, figsize=(10, 4.6))
for ax, carrier, title in zip(axes, ['h', 'e'], ['Holes', 'Electrons']):
    x = df[f'T{carrier}_CP_meV'].abs()
    y = df[f'T{carrier}_MOO_meV'].abs()


    red_cross_set = ids_keep.isin(flag) #ie builds boolean column about which ids have been flagged
    meets_any_rule = pd.Series(False, index=df.index)   # start with nothing flagged
    for rule in flag_rule.get(carrier, []):
        meets_any_rule = meets_any_rule | rule(df)
    #|means or
    is_flagged = red_cross_set | meets_any_rule

    plot_x, plot_y = x[~is_flagged], y[~is_flagged]
    cross_x, cross_y = x[is_flagged], y[is_flagged]

    sc = ax.scatter(plot_x, plot_y, c=cvals[~is_flagged], s=18, cmap='plasma', norm=norm)
    if is_flagged.any():
        ax.scatter(cross_x, cross_y, marker='x', c='red', s=45, lw=1.3, zorder=3,
                   label=f'flagged ({is_flagged.sum()})')

    # stats on the kept, un-flagged points
    m = (~is_flagged) & (x > 0) & (y > 0)
    r_log = np.corrcoef(np.log10(x[m]), np.log10(y[m]))[0, 1]

    pos = pd.concat([x, y])
    pos = pos[pos > 0] #so we only keep non-zero couplings
    if use_log:
        lo, hi = pos.min()/2, pos.max()*2
    else:
        lo, hi = pos.min() / 2, pos.max() *1.1

    ax.plot([lo, hi], [lo, hi], 'k--', lw=0.8, label='y = x')


    scale = 'log' if use_log else 'linear'
    ax.set(xscale=scale, yscale=scale, xlim=(lo, hi), ylim=(lo, hi),
            xlabel='|J| counterpoise (meV)', ylabel='|J| MOO (meV)',
            title=f"{title}  (r_log = {r_log:.2f}, N = {m.sum()})")
    ax.legend(loc='upper left', frameon=False)
    fig.colorbar(sc, ax=ax, label=colour_label,
                 extend='both' if clip_colour else 'neither')


fig.suptitle(RUN)
fig.tight_layout()
out_png = os.path.join(run_dir, 'Ts_compare.png')
fig.savefig(out_png, dpi=150)
print(f"Wrote {out_png}")


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
    x = df[f'T{carrier}_CP_meV'].abs()
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
#%%
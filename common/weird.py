# For inspection of weird molecules
#%%
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
RUN = 'Y6_20000ps_5A'
run_dir = os.path.expanduser(f'~/MSci/runs/{RUN}')

"""['init', 'fin', 'fin_lattice', 'dist',
            'Th_MOO_meV', 'Th_CP_meV', 'frac_diff_h',
            'Te_MOO_meV', 'Te_CP_meV', 'frac_diff_e',
            'time_ps']"""

both  = pd.read_csv(os.path.join(run_dir, 'Ts_compare.csv'))

#%%
extreme = both.loc[both["Th_CP_meV"].abs().idxmax()]
print(extreme)
#%%
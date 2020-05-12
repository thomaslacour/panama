import mst.calculKeff as mst
from material import display as matplt
import matplotlib.pyplot as plt
import numpy as np

from mst.calculKeff import calculKeff

# what result to print from the folowing list :
# alpha_eff, c_eff, n_eff, rho_eff, M_eff, Z_eff
key = 'alpha_eff'


out = calculKeff(
    Model     = 'WT',
    nMode     = 5,
    Mat       = 'panama_PU',
    Inc       = 'billes',
    r_mean_mm = 0.5,
    phivol    = 6.67,
    poly      = 10,
    d         = 2,
    freq      = np.linspace(0, 1, 101), #MHz
    polar     = 'L',
    verbose   = 1,
)



def save_as_csv(list_of_thing):
    try:
        import pandas as pd
    except ModuleNotFoundError:
        return
    csvfile = "data.csv"
    my_df = pd.DataFrame(list_of_thing)
    my_df = my_df.T
    my_df.to_csv(csvfile, index=False, header=['freq','alpha','velocity'])

result = [ out['freq'], out[key] ]
fig = plt.subplot(1,1,1)
fig.axes.plot(*result)
plt.show()
result = [ out['freq'], out['alpha_eff'], out['c_eff'] ]
save_as_csv(result)

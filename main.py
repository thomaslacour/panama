import os, sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from rt import actions as action

# show arguments passed to the script
print(sys.argv)

# COMPUTATION BASED ON A CONFIG FILE
if len(sys.argv)==2 and sys.argv[1][-3:]=='txt':
    config_file = sys.argv[1]
    data = action.read(config_file, v=0) # reading config file (v=1 for verbose)
    do = data["todo"] if data["todo"] is not None else 0 # assign default action
    result = action.process_config(data) # computation
    action.save_as_csv(result, "data.csv")
    plt.show()

else:
    pass

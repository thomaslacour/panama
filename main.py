import os, sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from rt import actions as action

# show arguments passed to the script
print(sys.argv)

# check if a config file is specified
if len(sys.argv)==2: #and sys.argv[1][-3:]!='txt':
    config_file = sys.argv[1]
else:
    config_file = 'conf.txt'

#reading config file (v=1 for verbose)
data = action.read(config_file, v=0)

# assign default action if none is given
# default:0 (rt vs frequency)
do = data["todo"] if data["todo"] is not None else 0

# computation
result = action.process_config(data)

if __name__ == "__main__":
    plt.show(block=True)
    action.save_as_csv(result[0], "data.csv")
    pass

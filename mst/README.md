Name
----

mst.py - compute the effective properties of a heterogeneous material.

Synopsis
--------

   python3 mst.py

Dependence
----------

These python packages are required:

- matplotlib, numpy, math, cmath, re, csv, pandas


### Installation of dependences

If they are not already installed, you should install them with the command

   python3 -m pip install --user matplotlib numpy pandas cmat mat re csv

Description
-----------

The MST module has the following folder architecture:

      src/
      ├ material/
      ├ mst/
      │   ├ calculKeff.py
      │   ├ calculScatteringCoeffSn.py
      │   ├ display.py
      │   └ special.py
      └ mst.py

After its execution, the python script `mst.py` will return a file named `mst.csv` where the output is stored.

See also
--------

> rt.py

Bugs
----

...

Author
------

[Thomas Lacour], [Olivier Poncelet]

  [Thomas Lacour]: mailto:thomas.lacour@u-bordeaux.fr
  [Olivier Poncelet]: mailto:olivier.poncelet@u-bordeaux

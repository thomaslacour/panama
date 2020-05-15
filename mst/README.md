NAME
----

mst.py - compute the effective properties of a heterogeneous material.

SYNOPSIS
--------

> python3 mst.py

DEPENDENCE
----------

These python packages are required:

- matplotlib, numpy, math, cmath, re, csv, pandas

If they are not already installed, you should install them with the command

> python3 -m pip install --user matplotlib numpy pandas cmat mat re csv

DESCRIPTION
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

SEE ALSO
--------

> rt.py

BUGS
----

...

AUTHOR
------

[Thomas Lacour], [Olivier Poncelet]

  [Thomas Lacour]: mailto:thomas.lacour@u-bordeaux.fr
  [Olivier Poncelet]: mailto:olivier.poncelet@u-bordeaux

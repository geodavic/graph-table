# graph-table
Code to render graphs for Roman et. al. juggling project.

# Install

To install, run
```
make install
```
This will install the required python packages.


# Use

To use, ensure the files `matrices.fits`, `subiso.csv`, and `polys.csv` are in the current directory. These are the adjacency matrices, subisomorphism lists, and polynomial lists for the table.

The syntax is demonstrated through the following example:
```
python3 graph.py --scale=0.7 --limit=8 --angle=0
```

This will render the tikz code for the table to `input.tex`. The scale parameter controls how big the graphs are (default is 1), the limit parameter controls how many graphs to draw (default is all of them), and the angle parameter is the angle in degrees at which the graphs are aligned (default is 45).

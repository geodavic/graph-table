# graph-table
Code to render graphs for Berens et. al. juggling project.

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

# Rendering a single graph

If you would like to just render TikZ code for a single graph given an adjacency matrix file, use:

```
python3 single_graph.py --file matrix.fits
```

This will print the TikZ code to your terminal. This has the same `scale` and `angle` options as `graph.py`. It accepts both `.csv` and `.fits` input file types.

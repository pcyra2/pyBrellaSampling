#!/bin/bash
wham  1.3 3.95 53 1e-06 300.0 0 prodmetadata.dat out.pmf 10 60
sed '1d' out.pmf | awk '{print $1,"",$2}' > plot_free_energy.dat


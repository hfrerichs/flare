#!/usr/bin/env gnuplot
set terminal X11 persist

set style line 1 linecolor 3
set style line 2 linecolor 3 linewidth 4

set xlabel "node space: {/Symbol h}"
set ylabel "real space: {/Symbol r}"

set terminal postscript eps enhanced color solid
set size 1.0, 0.5
set output "plot_DELTA_R_SYM.eps"

set multiplot
set size 0.5, 0.5

set origin 0.0, 0.0
set label "{/Symbol D} = 10 %, R = 4" at 0.1,0.7
plot \
"../spacing_DELTA_R_SYM_0.1_4.plt" t "" w l ls 1, \
"../spacing_DELTA_R_SYM_0.1_4.dat" t "" w l ls 2
unset label

set origin 0.5, 0.0
set label "{/Symbol D} = 20 %, R = 8" at 0.1,0.7
plot \
"../spacing_DELTA_R_SYM_0.2_8.plt" t "" w l ls 1, \
"../spacing_DELTA_R_SYM_0.2_8.dat" t "" w l ls 2

unset multiplot

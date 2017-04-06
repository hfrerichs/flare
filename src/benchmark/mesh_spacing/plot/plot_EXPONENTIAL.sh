#!/usr/bin/env gnuplot
set terminal X11 persist

set style line 1 linecolor 3
set style line 2 linecolor 3 linewidth 4

set xlabel "node space: {/Symbol h}"
set ylabel "real space: {/Symbol r}"

set terminal postscript eps enhanced color solid
set size 1.0, 0.5
set output "plot_EXPONENTIAL.eps"

set multiplot
set size 0.5, 0.5

set origin 0.0, 0.0
set label "{/Symbol l} = 0.5" at 0.1,0.9
plot \
"../spacing_EXPONENTIAL_0.5.plt" t "" w l ls 1, \
"../spacing_EXPONENTIAL_0.5.dat" t "" w l ls 2
unset label

set origin 0.5, 0.0
set label "{/Symbol l} = 0.125" at 0.1,0.9
plot \
"../spacing_EXPONENTIAL_0.125.plt" t "" w l ls 1, \
"../spacing_EXPONENTIAL_0.125.dat" t "" w l ls 2

unset multiplot

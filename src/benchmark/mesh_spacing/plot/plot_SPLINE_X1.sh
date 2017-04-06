#!/usr/bin/env gnuplot
set terminal X11 persist

set style line 1 linecolor 3
set style line 2 linecolor 3 linewidth 4
set style line 3 linecolor 1 linewidth 2

set xlabel "node space: {/Symbol h}"
set ylabel "real space: {/Symbol r}"

set terminal postscript eps enhanced color solid
set size 1.0, 0.5
set output "plot_SPLINE_X1.eps"

set multiplot
set size 0.5, 0.5

set origin 0.0, 0.0
set label "{/Symbol h}_1 = 0.8, {/Symbol r}_1 = 0.2, {/Symbol b} = 2" at 0.1,0.85
plot \
"../spacing_SPLINE_X1_0.8_0.2.plt" t "" w l ls 1, \
"../spacing_SPLINE_X1_0.8_0.2.dat" t "" w l ls 2, \
"../spacing_X1_0.8_0.2.dat"        t "" w l ls 3
unset label

set origin 0.5, 0.0
set label "{/Symbol h}_1 = 0.9, {/Symbol r}_1 = 0.7, {/Symbol b} = 2" at 0.1,0.85
plot \
"../spacing_SPLINE_X1_0.9_0.7.plt" t "" w l ls 1, \
"../spacing_SPLINE_X1_0.9_0.7.dat" t "" w l ls 2, \
"../spacing_X1_0.9_0.7.dat"        t "" w l ls 3

unset multiplot

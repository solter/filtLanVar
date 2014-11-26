set term post enh 		  # enhanced PostScript, essentially PostScript
 		 		  # with bounding boxes
set out 'd20.eps'                 # output file

set size 0.8, 0.8

set key right top
set title "d=20"
set xlabel '{/Symbol l}'
set ylabel 'f({/Symbol l})'
set xrange[:3]
set yrange[:1.1]
set parametric
x1=0.5
x2=1.0
y= 0.058081
d=0.02
h=0.01
set arrow from x1+d,y+h to x1-d,y-h nohead
set arrow from x1-d,y+h to x1+d,y-h nohead
set arrow from x2+d,y+h to x2-d,y-h nohead
set arrow from x2-d,y+h to x2+d,y-h nohead

plot 'd20.dat' using 1:3 with line title 'base filter {/Symbol y}({/Symbol l})', \
     'd20.dat' using 1:4 with line title 'polynomial filter {/Symbol r}({/Symbol l})', \
     'd20.dat' using 1:2 with line title '{/Symbol g}'

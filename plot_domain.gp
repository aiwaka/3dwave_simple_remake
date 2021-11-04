reset
set title "domain"
set xrange [-1:1]
set yrange [-1:1]
set zrange [0:2]
set view equal xyz
splot "mesh_gnup.dat" w l

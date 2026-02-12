set datafile columnheaders
set terminal png size 1280, 960 font "Libertinus Sans, 34"

set xlabel "Iteration"
set xtics 1 
set ylabel "Energy"

set key on font ",28"

set title "Harmonic oscillator"
set output "ho_conv.png"
set xrange [0.9:11.1]
plot "../output/ho_convergence.txt" using 1:3:4 \
      with yerrorlines title "VMC energy" pt 1, \
      0.5 title "Exact energy"

set title "Hydrogen atom"
set output "hy_conv.png"
set xrange [0.9:3.1]
plot "../output/hy_convergence.txt" using 1:3:4 \
      with yerrorlines title "VMC energy" pt 1, \
      -0.5 title "Exact energy"

# set title "Helium atom"
# set output "he_conv.png"
# set xrange [0.9:3.1]
# plot "../output/he_convergence.txt" using 1:3:4 \
#       with yerrorlines title "VMC energy" pt 1, \
#       0.5 title "Exact energy"
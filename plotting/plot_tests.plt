set datafile columnheaders
set terminal png size 1280, 960 font "Libertinus Sans, 34"

set xlabel "Alpha"
set ylabel "Energy"

set key on font ",28"

set title "Harmonic oscillator"
set output "ho_fixalpha.png"
set xrange [0.39:0.61]
plot "../output/ho_fixedalpha.txt" using 1:2:3 \
      with yerrorlines title "Results" pt 5, \
      "../output/ho_thijssen.txt" using 1:2:3 \
      with yerrorlines title "Thijssen 2007" pt 5

# set title "Hydrogen atom"
# set output "hy_fixalpha.png"
# plot "../output/hy_convergence.txt" using 1:3:4 \
#       with yerrorlines title "Results" pt 1, \
#       0.5 title "Thijssen 2007"

# set title "Helium atom"
# set output "he_fixalpha.png"
# plot "../output/he_fixedalpha.txt" using 1:2:3 \
#       with yerrorlines title "Results" pt 1, \
#       "../output/he_thijssen.txt" using 1:2:3 \
#       with yerrorlines title "Thijssen 2007" pt 1
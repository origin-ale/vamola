set datafile columnheaders
set terminal png size 1280, 960 font "Libertinus Sans, 34"

set xlabel "Alpha"
set ylabel "Energy"

set key on font ",28"

set output "ho_fixalpha.png"
set xrange [0.39:0.61]
plot "../output/ho_fixedalpha.txt" using 1:2:3 \
      with yerrorlines title "Results (h.o.)" pt 3, \
      "../thijssen/ho_thijssen.txt" using 1:2:3 \
      with yerrorlines title "Thijssen 2007" pt 2

set output "hy_fixalpha.png"
set xrange [0.79:1.21]
plot "../output/hy_fixedalpha.txt" using 1:2:3 \
      with yerrorlines title "Results (H)" pt 1, \
      "../thijssen/hy_thijssen.txt" using 1:2:3 \
      with yerrorlines title "Thijssen 2007" pt 1

set output "he_fixalpha.png"
set xrange [0.04:0.255]
set yrange [-2.9:-2.85]
plot "../output/he_fixedalpha.txt" using 1:2:3 \
      with yerrorlines title "Results (He)" pt 1, \
      "../thijssen/he_thijssen.txt" using 1:2:3 \
      with yerrorlines title "Thijssen 2007" pt 1
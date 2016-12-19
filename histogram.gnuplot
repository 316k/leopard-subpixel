set boxwidth 0.9 relative
set style data histograms
set style fill solid 1.0 border -1

set xtics rotate

plot filename u 1:xtic(2)

pause mouse close

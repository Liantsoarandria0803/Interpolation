set title 'Interpolation Lineaire'
set xlabel 'X'
set ylabel 'Y'
plot 'plot_data.txt' using 1:2 with points title 'Data', 'plot_data.txt' using 1:3 with lines title 'Linear Interpolation'

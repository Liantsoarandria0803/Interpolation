set terminal png size 800,600
set output 'interpolation.png'
set title 'Interpolation'
set xlabel 'X'
set ylabel 'Y'
set grid
plot 'data.txt' using 1:2 with points pt 7 title 'Data', 'linear_data.txt' with lines title 'Linear Interpolation', 'spline_data.txt' with lines title 'Cubic Spline', 'lagrange_data.txt' with lines title 'Lagrange Interpolation', 'data.txt' using 1:2 smooth csplines title 'Gnuplot Cubic Spline'

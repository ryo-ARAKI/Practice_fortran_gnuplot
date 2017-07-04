#Reset all settings
reset
#No regend(凡例なし)
set nokey
#Set gif animation as output
#optimize: use one color map in animation
#delay 10: show one frame for 10/100 second
#size 960,960: size of output
set term gif animate optimize delay 10 size 960,960
#Set the name of output file
set output 'movie_Advection_Nonlinear_1D_Explicit.gif'


#Set range for each axis
set xrange[0:1]
set yrange[0:1]
#Set color bar range

#Set variables

do for [n = 1:250] {
#
plot "Advection_Nonlinear_1D_Explicit.dat" index n using 1:2 with lines # n番目のデータのプロット

}

set out
set terminal wxt enhanced

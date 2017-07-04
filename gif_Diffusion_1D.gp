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
set output 'movie_1D_Euler_Dirichlet.gif'
#set output 'movie_1D_Crank-Nicolson_Dirichlet.gif'


#Set range for each axis
set xrange[0:1]
set yrange[0:0.2]
#Set color bar range

#Set variables

do for [n = 1:10] {
#
plot "Diffusion_1D_Euler_Dirichlet.dat"  index n using 1:2 with lines # n番目のデータのプロット
#plot "Diffusion_1D_Crank-Nicolson_Dirichlet.dat"  index n using 1:2 with lines # n番目のデータのプロット

}

set out
set terminal wxt enhanced

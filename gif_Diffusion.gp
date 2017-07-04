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
set output 'movie_Diffusion_explicit.gif'
#set output 'movie_Diffusion_implicit.gif'


#Set 3-D display using color map
set pm3d
#Set angle
set view 60,120
#Set range for each axis
set xrange[0:1]
set yrange[0:1]
set zrange[0:1]
#Set color bar range
set cbrange[-1:1]
#set isosamples 50   #Number of lines

#Set variables
n0 = 1    # ループ変数の初期値
n1 = 40   # ループ変数の最大値
dn = 1    # ループ変数の増加間隔

do for [n = n0:n1] {
#
splot "Diffusion_explicit.d"  index n using 1:2:3 with lines # n番目のデータのプロット
#splot "Diffusion_implicit.d"  index n using 1:2:3 with lines # n番目のデータのプロット
}

set out
set terminal wxt enhanced

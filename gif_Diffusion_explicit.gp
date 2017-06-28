reset
set nokey
set term gif animate optimize delay 10 size 480,480
set output 'movie_Diffusion_explicit.gif'

set pm3d
set view 60,120
set xrange[0:1]
set yrange[0:1]
set zrange[0:1]
set cbr[-1:1]        #Range of color bar
#set isosamples 50   #Number of lines

#Set variables
n0 = 1    # ループ変数の初期値
n1 = 40   # ループ変数の最大値
dn = 1    # ループ変数の増加間隔

if(exist("n")==0 || n<0) n = n0  # ループ変数の初期化

do for [n = n0:n1] {
#Plot
splot "output_ex6_19.d"  index n using 1:2:3 with lines # n番目のデータのプロット
}

set out
set terminal wxt enhanced

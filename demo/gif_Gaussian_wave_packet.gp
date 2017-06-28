set term gif animate optimize delay 2 size 480,360
set output 'movie_Gaussian_wave_packet.gif'
set xrange[-5:10]
set yrange[0:1]

do for [i = 0:400 ] {
   t=i*0.02
   plot sqrt(1/(1+t*t))*exp(-(x-t)**2/(1+t*t)) lw 2
   }

set out
set terminal wxt enhanced

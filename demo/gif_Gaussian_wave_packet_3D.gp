set term gif animate optimize delay 10 size 480,480
set output 'movie_Gaussian_wave_packet_3d.gif'

set pm3d at b
set xr[-0:1]
set yr[-0:1]
set zr[0:1]
set cbr[0:1]        #Range of color bar
set isosamples 50   #Number of lines

do for [i = 0:50 ] {
   t=i*0.05
   splot sqrt(1/(1+t*t))*exp(-(x-t)**2/(1+t*t))*sqrt(1/(1+t*t))*exp(-(y-2*t)**2/(1+t*t))
   }

set out
set terminal wxt enhanced

do for [i = 0:400 ] {
   t=i*0.02
   plot sqrt(1/(1+t*t))*exp(-(x-t)**2/(1+t*t))
   }

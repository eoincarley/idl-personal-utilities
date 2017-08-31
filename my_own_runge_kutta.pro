function dvdt, T

    ; Define the derivative that is to be numerically integrated 
    ; In this case acceleration
    ; T is time in seconds.

	COMMON share, ar, ad, tr, td
	dvdt = dblarr(n_elements(T))
	dvdt[*] = 1.0/(  1.0d/(ar*exp((T)/tr))  + 1.0d/(ad*exp(-(T)/td) )  )

return, dvdt
END
; These two could be integrated....
function dxdt, xt_integ, T

    ;********Define the derivative that is to be numerically integrated******
    ; In this case velocity
    ; T is time in seconds.

	dxdt = dblarr(n_elements(T))
	index = closest(xt_integ[0,*], T)
	dxdt[*] = xt_integ[1,index[0]]

return, dxdt
END

;-------------------------------------------------------;
;					Main Procedure					    ;
;					
; Construction of my own 4th order Runge-Kutta algorithm to 
; integrate the acceleration expression in Byrne et al. (2010).
								    
pro my_own_runge_kutta

	COMMON share, ar, ad, tr, td
	ar = 1.0       ; See Byrne et al. (2010) for these constants and starting params.
	tr = 138.0     ; e-folding 1 time for acceleration increase (jerk)
	ad = 4950.0
	td = 1249.0    ; e-folding 2 time for acceleration increase


    ;*******************Initial Values an variables***********************
    T = 0.0d        ;6.625  (Shane's Thesis)   ;time units
    V = 40.0e3      ;272.5  (Shane's Thesis) ;km/s
    H = 10.0        ;Step Size for Runge-Kutta
    X = 25.0e6  
    nsteps =  (4.0*3600.0-T)/H  ;time step
    vt_integ = dblarr(2, nsteps)
    tsim = (dindgen(1001)*(4.0*3600.0 - 0.0 )/1000.0 ) +0.0 
    a = 1.0/(  1.0d/(ar*exp(tsim/tr))  + 1.0d/(ad*exp(-tsim/td) )  )
    tstart = anytim(file2time('20020421_004000'), /utim)
         
    ;-----------------------;
    ;   Plot acceleration
    ;
    window,0
    utplot, tstart+tsim, a, charsize=1.5, ytitle='CME Acceleration m/s/s'

    ;-----------------------------;
    ;      Runge-Kutta Loop
    FOR i=0.0, nsteps-2 DO BEGIN
        vt_integ[0,i] = T
        vt_integ[1,i] = V

        ;Compute k1-4 for Runge-Kutta interation
        k1 = dvdt( vt_integ[0,i] );, vt_integ[1,i] )
        k2 = dvdt( vt_integ[0,i]  + H/2.0d);, vt_integ[1,i] + k1/2.0d )
       	k3 = dvdt( vt_integ[0,i]  + H/2.0d);, vt_integ[1,i] + k2/2.0d )
        k4 = dvdt( vt_integ[0,i] ); + H, vt_integ[1,i] + k3)
        T = vt_integ[0,i] + H
        V = vt_integ[1,i] + (H/6.0d)*(k1 + 2.0d*(k2+ k3) + k4)      ;*!dpi ;factor of pi comes in somehwere in 
        															;Vrsnak 2002 JGR 107. Not sure where
        															;exactly
        print, V
    ENDFOR    
    window,1
    utplot, tstart+vt_integ[0,*], vt_integ[1,*]/1.0e3, /ylog, ytitle='CME Velocity (km s!U-1!N)', charsize=1.5,$
    yr=[1,1.0e4], /ys, /xs


    T = 0.0d    
    xt_integ = dblarr(2,nsteps)

    ;-----------------------------;
    ;      Runge-Kutta Loop
    ;      Integrating dxdt (from the last loop). Could be nested in 1st loop.
    FOR i=0.0, nsteps-2 DO BEGIN
        xt_integ[0,i] = T
        xt_integ[1,i] = X

        ;Compute k1-4 for Runge-Kutta interation
        k1 = dxdt( vt_integ, xt_integ[0,i] );, vt_integ[1,i] )
        k2 = dxdt( vt_integ, xt_integ[0,i]  + H/2.0d);, vt_integ[1,i] + k1/2.0d )
       	k3 = dxdt( vt_integ, xt_integ[0,i]  + H/2.0d);, vt_integ[1,i] + k2/2.0d )
        k4 = dxdt( vt_integ, xt_integ[0,i] ); + H, vt_integ[1,i] + k3)
        
        T = xt_integ[0,i] + H
        X = xt_integ[1,i] + (H/6.0d)*(k1 + 2.0d*(k2+ k3) + k4)      ;*!dpi ;factor of pi comes in somehwere in 
        															;Vrsnak 2002 JGR 107. Not sure where
        															;exactly
        print,X
    ENDFOR  

    ;------------------------;
    ;   Plot the results
    ;  
    window,2
    utplot, tstart +xt_integ[0,*], xt_integ[1,*]/1.0e6, /ylog, charsize=1.5, yr=[10,1e5],/xs,/ys,$
    ytitle='CME Height (Mm)'


    window,3
    plot, xt_integ[1,*]/6.955e8 +1.0, vt_integ[1,*]/1.0e3, /ylog, charsize=1.5, yr=[10.0, 1.0e4], /xs, /ys,$
    xtitle='CME Height (R/R!Lsun!N)', ytitle='CME Velocity (km s!U-1!N)', xr=[1,5]


END
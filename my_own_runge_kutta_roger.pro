FUNCTION dvdt,T

    ;********Define the derivative that is to be numerically integrated******

		COMMON share, ar, ad, tr, td

	dvdt = dblarr(n_elements(T))
	
	
	dvdt[*] = 1.0/(  1.0d/(ar*exp((T-1000.0)/tr))  + 1.0d/(ad*exp(-(T-1000.0)/td) )  )

return, dvdt
END

pro my_own_runge_kutta_roger_v2

        ;Following common black param values from Shane's Thesis or Byrne et al. (2010)
		COMMON share, ar, ad, tr, td
    			;ar = 2000.0d             ;Assymptotic value of solar wind velocity 
				;ad = 4950.0d3			   
				;tr = 138.0d				 
				;td = 1249.0
				ar=1.0
				tr=138.0
				ad = 4950.0
				td = 1249.0

     
;********************************Plot***************************************  

tsim = (dindgen(1001)*(4.0*3600.0 - 1.0 )/1000.0 ) +1.0 
a = 1.0/(  1.0d/(ar*exp(tsim/tr))  + 1.0d/(ad*exp(-tsim/td) )  )

plot,tsim,a
;stop



;*******************Initial Values***********************
T = 0.0d     ;6.625  (Shane's Thesis)
V = 40.0e3     ;272.5  (Shane's Thesis)

H = 0.1 	   ;Step Size for Runge-Kutta
nsteps =  (4.0*3600.0-T)/H

vt_integ = dblarr(2,nsteps)

;Runge-Kutta Loop
FOR i=0.0, nsteps-2 DO BEGIN
    vt_integ[0,i] = T
    vt_integ[1,i] = V
    stop
    ;Compute k1-4 for Runge-Kutta interation
    k1 = dvdt( vt_integ[0,i] );, vt_integ[1,i] )
    
    k2 = dvdt( vt_integ[0,i] ); + H/2.0d, vt_integ[1,i] + k1/2.0d )
 
   	k3 = dvdt( vt_integ[0,i] ); + H/2.0d, vt_integ[1,i] + k2/2.0d )

    k4 = dvdt( vt_integ[0,i] ); + H, vt_integ[1,i] + k3)
    
    T = vt_integ[0,i] + H
    V = vt_integ[1,i] + (H/6.0d)*(k1 + 2.0d*(k2+ k3) + k4);*!dpi ;factor of pi comes in somehwere in 
    															;Vrsnak 2002 JGR 107. Not sure where
    															;exactly
    print,V
ENDFOR    


window,1


;tsim = (dindgen(nsteps+1.0)*(4.5*3600.0 +3600.0 )/(nsteps) ) -3600.0  
;a = 1.0/(  1.0d/(ar*exp((tsim-1000.0)/tr))  + 1.0d/(ar*exp(-(tsim-1000.0)/td) )  )
;s = dblarr(n_elements(vt_integ[0,*]))


;s[*] = 40.0e3*vt_integ[0,*] ;+ 0.5*a*(vt_integ[0,*]^2.0)
plot, s/6.95e5, vt_integ[1,*]/1.0e3, /ylog, ytitle='Velocity (km s!U-1!N)', xr=[0,10]
stop

END
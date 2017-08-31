FUNCTION dvdr,R,V

    ;********Define the derivative that is to be numerically integrated******

		COMMON share,vsw_a,alpha,bet,gam 

	dvdr = dblarr(n_elements(V))
	vsw  = dblarr(n_elements(R))
	
	vsw[*]=vsw_a*sqrt( 1.0-exp( -1.0d*(R[*]-2.8d)/8.1d) ) ;Solar wind speed Sheeley et al.
	
	;Following expression from Byrne et al. (2010)
	dvdr[*] = ( 6.95508*(10.0d^(5.0d)) )*( alpha*(R[*]^double(bet) ) )*( 1.0/V[*] )*(abs( V[*]-vsw[*] )^gam)
	

return,dvdr
END

pro my_own_runge_kutta

        ;Following common black param values from Shane's Thesis or Byrne et al. (2010)
		COMMON share,vsw_a,alpha,bet,gam
    			vsw_a = 555.0d             ;Assymptotic value of solar wind velocity 
				alpha = 4.55d-5			   
				bet =  -2.025d				 
				gam =   2.275d				  

cd,'/Users/eoincarley/Data/secchi_l1/20081212/COR1+COR2'
;****************************Read in kinematics data******************************
readcol, 'kins_meanspread.txt', mid_a, mid_r, top_a, top_r, bottom_a, bottom_r, $
midtop_a, midtop_r, mibottom_a, midbottom_r, time, f='D,D,D,D,D,D,D,D,D,D,A'

	   utbasedata = anytim(time[0])
       t = anytim(time) - utbasedata
       rsun = 6.95508e8

       mid_km = mid_r * rsun / 1000.0d
       top_km = top_r * rsun / 1000.0d
       bottom_km = bottom_r * rsun / 1000.0d
       midtop_km = midtop_r * rsun / 1000.0d
       midbottom_km = midbottom_r * rsun / 1000.0d

       mid_vel = deriv(t, mid_km)
       top_vel = deriv(t, top_km)
       bottom_vel = deriv(t, bottom_km)
       midtop_vel = deriv(t, midtop_km)
       midbottom_vel = deriv(t, midbottom_km)

       mid_accel = deriv(t, mid_vel) * 1000.0d ; to be in Metres
       top_accel = deriv(t, top_vel) * 1000.0d
       bottom_accel = deriv(t, bottom_vel) * 1000.0d
       midtop_accel = deriv(t, midtop_vel) * 1000.0d
       midbottom_accel = deriv(t, midbottom_vel) * 1000.0d
       
       
;********************************Plot***************************************  

plot,mid_r,mid_vel,xtitle='Heliocentric Distance [R/R!Lsun!N]',ytitle='Velocity [km/s]',$
psym=1
save,mid_r,mid_vel,filename='CMEvr_20081212_Roger.sav'

;*******************Initial Values***********************
R = 6.625d     ;6.625  (Shane's Thesis)
V = 272.5d     ;272.5  (Shane's Thesis)

H = 0.001 	   ;Step Size for Runge-Kutta
nsteps =  (200.0-R)/H

vr_integ = dblarr(2,nsteps)


;Runge-Kutta Loop
FOR i=0.0, nsteps-1 DO BEGIN
    vr_integ[0,i] = R
    vr_integ[1,i] = V
    
    ;Compute k1-4 for Runge-Kutta interation
    k1 = dvdr( vr_integ[0,i], vr_integ[1,i] )
    
    k2 = dvdr( vr_integ[0,i] + H/2.0d, vr_integ[1,i] + k1/2.0d )
 
   	k3 = dvdr( vr_integ[0,i] + H/2.0d, vr_integ[1,i] + k2/2.0d )

    k4 = dvdr( vr_integ[0,i] + H, vr_integ[1,i] + k3)
    
    R = vr_integ[0,i] + H
    V = vr_integ[1,i] + (H/6.0d)*(k1 + 2.0d*(k2+ k3) + k4)*!dpi ;factor of pi comes in somehwere in 
    															;Vrsnak 2002 JGR 107. Not sure where
    															;exactly
ENDFOR    

oplot,vr_integ[0,*], vr_integ[1,*]

;**************Save params and dvdr integration result******************
params='Initial Height='+string(vr_integ[0,0],format='(f5.2)')+', '+$
	   'Initial Velocity='+string(vr_integ[1,0],format='(f7.3)')+', '+$
	   'Asym Solar Wind Speed='+string(vsw_a,format='(f7.2)')+', '+$
	   'alpha='+string(alpha,format='(f11.8)')+', '+$
	   'beta='+string(bet,format='(f6.3)')+', '+$
	   'Gamma='+string(gam,format='(f6.3)')
	   print,params
save,vr_integ,params,filename='velocity_fit_20081212.sav'


END













PRO common_block_definitions

;
;  Define a number of default parameters
;
  COMMON programconstants, plottype,verbosity,symboltype,dfn,dpn, $
         numdegrees,startingpoint
     plottype=0            ;0: none1: overplots2: separate plots
     verbosity=1           ;0: none1: small data2: most data3: all 360 pb traces
     symboltype=3          ;0:-1:+2:*3:.4:<>5:^6:[]7:x
     dfn=''    			;default filename
dpn=getenv('POLDIR')
;     dpn=''			;default pathname
     numdegrees=360.
     startingpoint=0	;  -131 (andrew's orig value to avoid c3 pylon)

;
;   Inverse Polynomial EXPonentS, min pb CutOff Value
;	Can have multiple trial functional forms for the Inverse Polynomial
;
  COMMON simconstants, ipexponents, weightstrength, cov, nrad, maxrad, polyval ,hgtval

polyval=4	;default value for the order of ipexponents ([[-1,-2,-3,-4]])
hgtval=5.	;default value for the maximum height to convert

	ipexponents = TRANSPOSE([[-1,-2,-3,-4,-5,-6,-7,-8,-9]])
          weightstrength=1.5       ;weights = pb^(-weightstrength)
         cov=1e-14                ;pb values less than cov are ignored in curvefit

  COMMON widgetvars, xss, yss, pbpic, header, windowid, suncenter, widgetu, widgetws, $
                     widgetcov, xs, ys, origline, sunradius, filenames, pltwinid, $
                     logpbpic, dispplt, rhos, pbdata, simne, simpb, clangle, angle, rstart, $
                     nearr, minnef,xarr,yarr,lev,xysize,picwin,pltwin, pbarr, $
                     anglearr, polartype, rhosorig, rmin,widgetres
;
;  sunradius is solar radius in pixels/radius
;  rmin is minimum solar radius to fit polynomials to
;
   xss=512          ;x windowsize of the pic 512
   yss=512          ;y windowsize of the pic 512
   xysize=1024      ; size of the scrollable window

   widgetres=3.		;# of degrees/trace

  COMMON pbconstant, u
   u=.58                                               ;limb darkening coefficient

RETURN
END




PRO PB_INVERTER_MOD,fname,pb,hdr

;Modified by Eoin Carley so coefficients and exponents for
;a number of radial traces can be output as a text file called coeffs_exps.txt
;These coeffs and exponents are used in vdh_inversion to get the radial density profile.

print,'**************'
print,'LOADIMG'
print,'**************'
;
;  Read the pB image, and subtract off the F-coronal contribution
;

COMMON widgetvars
COMMON simconstants
cd,'/Users/eoincarley/Data/22sep2011_event/LASCO_C2/polarized_brightness'
fname='C2-PB-20110922_0257.fts'

pb=LASCO_READFITS(fname,hdr)

tel=''
IF DATATYPE(hdr) ne 'STC' then BEGIN
   tel = strtrim(strupcase(sxpar(hdr, 'TELESCOP')),2)
   if (tel eq 'MK3') or (tel eq 'MK4') then hdr=mlo_fitshdr2struct(hdr) else hdr=lasco_fitshdr2struct(hdr)
ENDIF ELSE tel = hdr.TELESCOP

IF (pb(0) NE -1)  THEN BEGIN
   suncenter=GET_SUN_CENTER(hdr)
   if (tel eq 'MK3') or (tel eq 'MK4') then sunradius = hdr.crradius else $
	   sunradius=GET_SOLAR_RADIUS(hdr,/pixel)
   ELIMINATE_FCORONA, pb, suncenter, sunradius, filenames
   CASE STRUPCASE(hdr.detector) OF
   'C2':    BEGIN
               rmin = 2.35 ;original 2.25 but this introduced some funny artifacts (Eoin Carley)
            END
   'C3':    BEGIN
               rmin = 4.0
            END
   'MK3':   BEGIN
	       rmin = 1.15
	    END
   'MK4':   BEGIN
	       rmin = 1.15
	    END 
  else:
   ENDCASE
;
;  find mininum value for the greatest radius in each of the 4 directions from sun center
;  this will be the number of points in the radius array
;

   sz = SIZE(pb)
  nrad = suncenter.xcen<(sz(1)-suncenter.xcen) < suncenter.ycen<(sz(2)-suncenter.ycen)

   nrad = ROUND(nrad)
  
ENDIF

logpbpic = ALOG(pbpic>widgetcov)
plot_image,logpbpic

;The following produces polarised brightness profiles for the angles
;given by the for loop. These profiles are then used in ne_from_pb
;to produces the coefficients and exponents. These coefficients and
;exponents are then used in vdh_inversion to get the density profile.
;Eoin Carley

FOR i=90.0,100.0,0.01 DO BEGIN    
   j=i+90.0
  eta = j*!DTOR
  angle = eta*!radeg-90		; store angle in common
  ;IF i eq 70.0 THEN BEGIN
  ;	angle_save = angle
  ;ENDIF ELSE BEGIN
  ; angle_save = [angle_save,angle]
  ;ENDELSE
  
  IF (angle LT 0)  THEN angle=angle+360
  clangle = angle
  rhos = FINDGEN(nrad)
  xs = COS(eta) * rhos + suncenter.xcen
  ys = SIN(eta) * rhos + suncenter.ycen
  rhos = rhos/sunradius			;now in solar radii
  hgtval = hgtval<MAX(rhos)
  rhosorig = rhos
  origline = INTERPOLATE(pbpic,xs,ys)	; do bilinear interpolation to pick up the line
  pbdata =origline

;
;	Get rid of datapoints less than rmin
;
stop
  goodrhos = WHERE(rhos GT rmin,nw);	old version had rhos GT 1.
  rhos      = rhos(goodrhos)
  pbdata        = pbdata (goodrhos)
  goodxs   = xs(goodrhos)
  goodys   = ys(goodrhos)


;
;	get rid of data points on the occulting disk (assumes monotonically decreasing pB)
;
  trash    = max(pbdata , msub)			; msub is the index of the max value
  goodrhos = goodrhos(msub:*)
  rhos      = rhos(msub:*)
  pbdata        = pbdata (msub:*)
  goodxs   = goodxs(msub:*)
  goodys   = goodys(msub:*)
  weights = DBLARR(N_ELEMENTS(rhos)) ;same size as rhos & pbdata, but double
;
;	get rid of missing or invalid data points
;	widgetcov is the minimum value of pB considered to be valid
;
  good = WHERE(pbdata  GT widgetcov, ngood)

  IF (ngood NE 0) THEN BEGIN
     goodxs  = goodxs(good)
     goodys  = goodys(good)
     weights(good) = pbdata (good)^(-widgetws)		; widgetws is weight strength
     wht = WHERE (rhos GT hgtval,nwht)
     IF (nwht GT 0)   THEN weights(wht)=0
     minweight = MIN(weights(WHERE(weights GT 0)))
     weights = weights/minweight				; normalize weights to 1
     IF (KEYWORD_SET(doplot)) THEN BEGIN
        WSET, 0
        PLOTS,xs,ys;,color=origline
        PLOTS,goodxs,goodys,color=255
     ENDIF
  ENDIF ELSE BEGIN
     goodxs=[0]
     goodys=[0]
     IF (KEYWORD_SET(doplot)) THEN BEGIN
        WSET,0
        PLOTS,xs,ys;,color=origline
     ENDIF
   ENDELSE
   w = WHERE (rhos(0) EQ rhosorig,nw)
   IF (nw GT 0)  THEN rstart=w(0)  ELSE BEGIN
      PRINT,'ERROR: PB_INVERTER/GETRPW: No good values at PA = ',angle
      rstart=0
   ENDELSE
   PLOTS,xs,ys
   
pbmod = pbdata
save,pbmod,filename='pbmod.sav'
   
   
   
   simne=ne_from_pb(rhos,pbdata,weights,ipexponents(0,0:polyval-1),widgetu,simpb,var, angle=angle)

ENDFOR
RETURN
END


PRO eliminate_fcorona, pic, suncenter, sunradius, filenames

;
;  The VDH inversion requires pB to be only from electron scattering.  Any residual
;  pB from the F-Corona or stray light, for example, must be removed.
;
;  subtracts Koutchmy-Lamy F-corona model, interpolating between FpB data points and
;  Fb polar and ecliptic values

  s = SIZE(pic)
  rhos=distarr(s(1),s(2),suncenter.xcen,suncenter.ycen,dxs,dys);now in pixels
  rhos=rhos/sunradius;now in solar radii
  etas=ATAN(dys,dxs)
  fb = FCOR_KL(rhos, etas)		; koutchmy-lamy model of F-coronal brightness
  fp = FCORPOL_KL (rhos)		; koutchmy-lamy model of F-coronal polarization
  fpb = fb*fp
  pic=temporary(pic)-fpb
RETURN
END


FUNCTION ne_from_pb, rhos, pb, pbweights, exps, u, $
                             simpb, pbvar, pbcoefs, pbpick, angle=ANGLE, doplot=doplot
                    
;	Overall program to control the VDH inversion for a number of
;	trial functions
;	Returns the best fit to Ne
;
;INPUTS
;		Rhos:	Radii where pB is known
;		Pb:		known values of pB
;		Pbweights:	Weights at each point to be used in curve fit
;		Exps:	Exponents of polynomial
;		U:		Limb darkening
;		
;
;KEYWORDS
;		DOPLOT:	If set then plot various parameters, Default is not to plot
;
;OPTIONAL OUTPUTS:
;		Simpb:	Computed values of pB for the best fit
;		Pbvar:	Variance of the best pB fit
;		Pbcoefs:Coefficients of the best pB fit
;		Pbpick:	Pointer to the best pB fit functional form
;
  COMMON programconstants
  COMMON ipexpansionarg, ipexps, goodrhos


  numpoints = N_ELEMENTS(rhos)
  sz = SIZE(exps)
  numtrials = sz(1)
  numterms  = sz(2)

;  tsimpb     = DBLARR(numtrials,numpoints)
  tpbcoefs   = DBLARR(numterms)
  coefs      = DBLARR(numterms)+1

  IF (N_ELEMENTS(WHERE(pbweights GT 0)) LE 50) THEN RETURN,DBLARR(numpoints)+1
;
;	Fit the input data of pB vs Rhos using CURVEFIT for each of the trial polynomials
;   Compute the variance for each of the trials and save all the results
;

      ipexps        = REFORM(exps(0,*))
      tsimpb  = CURVEFIT (rhos,pb,pbweights,coefs,function_name='ipexpansion', $
                                itmax=500,tol=1e-11)
      tpbcoefs = coefs
      tpbvar    = variance(pb,tsimpb)  ;save the variance
;
;	Pick the best trial polynomial fit (minimum variance) and then do the inversion
;

  nsube  = VDH_INVERSION(rhos,u,tpbcoefs, ipexps)
;
IF file_test('coeffs_exp_90-100.txt') eq 1 THEN BEGIN
	readcol,'coeffs_exps_90-100.txt', a, c0, c1, c2, c3, ex0, ex1, ex2, ex3,$
	format = 'D,D,D,D,D,D,D,D,D'
	
	angle = [a,angle]
    coefs0 = [c0,coefs[0]]
	coefs1 = [c1,coefs[1]]
	coefs2 = [c2,coefs[2]]
	coefs3 = [c3,coefs[3]]
	 
    ipexps0 = [ex0,ipexps[0]]
	ipexps1 = [ex1,ipexps[1]]
	ipexps2 = [ex2,ipexps[2]]
	ipexps3 = [ex3,ipexps[3]]
	
	
ENDIF ELSE BEGIN
	angle = angle
	coefs0 = coefs[0]
	coefs1 = coefs[1]
	coefs2 = coefs[2]
	coefs3 = coefs[3]
	
	ipexps0 = ipexps[0]
	ipexps1 = ipexps[1]
	ipexps2 = ipexps[2]
	ipexps3 = ipexps[3]
	
ENDELSE	

writecol,'coeffs_exps_90-100.txt',angle,coefs0,coefs1,coefs2,coefs3,ipexps0,ipexps1,ipexps2,ipexps3,$
fmt = '(D,D,D,D,D,D,D,D,D)'


;	RETURN ONLY THE BEST PB TRIAL VARIANCE AND COEFFICIENTS
;
  IF (n_params() GE 6) THEN simpb   = tsimpb
  IF (n_params() GE 7) THEN pbvar   = tpbvar
  IF (n_params() GE 8) THEN pbcoefs = tpbcoefs
  IF (n_params() GE 9) THEN pbpick  = tpbpick
  RETURN, nsube
  
END
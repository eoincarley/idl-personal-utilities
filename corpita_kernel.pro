; Routine that determines the pixels to be analysed by the
; corpita_pipeline routine. This is unique to each event, but the
; same for each image studied for that event, allowing it to be
; calculated once, speeding up processing time. 
;+ 
; NAME:
;    CORPITA_KERNEL
;
; PURPOSE:
;    Calculate arc sectors for given CBF event to be used
;    with CorPITA algorithm. 
;
; INPUTS:
;    MAP_0 -> Base image map.
;    INDEX_0 -> Header of base image FITS file.
;    X_A, Y_A -> X- & Y-coordinates of source location
;    TILT_DEG -> Array of angles for which the processing must be
;                applied to the images (360 degrees by default).
;    TOT_ANGLE -> The arc size to be studied (10 degrees by default).
;    INCR -> The increment angle (1 degree by default)
;    ST_ANG -> The starting angle (0 degrees by default)
;    SAVEDIR -> Location of directory for save files.
;    T_START -> Flare start time
;
; OPTIONAL INPUTS:
;    VERB -> Set to produce continuous progress updates.
;
; OUTPUTS:
;    IND_LON -> Array of indices where arc reaches limb
;    INDOK -> Pointer array of good pixels for each arc
;    CNTOK -> Number of good pixels for each arc
;    REV_IND -> Pointer array of Reverse Indices from histogram
;               identifying good pixels.
;    BINCOUNT -> Histogram of pixels in each arc
;    ARC_IND -> Pointer of pixels for each arc
;
; VERSION: 1.4 (Updated from v1.3 01-Sep-2014; Updated from v1.2 03-Jun-2013; Updated from v1.1 21-Jan-2013)
;
; HISTORY:
;    Oct-2011 DL Written (first draft)
;    19-Dec-2011 DL Date & passband now passed as strings, not
;                   floats. 
;    10-Jan-2012 DL Header updated and keywords removed to reflect
;                   modifications to calling routine.
;    13-Mar-2012 DL Header and code clean-up.
;    19-Jun-2012 DL Changed date to YYYY/MM/DD format.
;    08-Jan-2013 DL Major revisions:
;                   - Now using histograms instead of masks to
;                     identify pixels for subsequent processing. 
;                   - Save-file name changed for brevity and clarity.
;    21-Jan-2013 DL Version updated to v1.2. Now using delvarx to
;                   remove unused variables from memory.  
;    21-Feb-2013 DL Routine renamed from pipeline_kernel.pro to
;                   corpita_kernel.pro for consistency. Now defining
;                   arc edge lines using pointers. 
;    25-Feb-2013 DL Obsolete keywords removed
;    22-Mar-2013 DL Code cleanup and removal of obsolete variables
;    11-Apr-2013 DL Modified to allow arcs to start at angles other than
;                   0 degrees
;    03-Jun-2013 DL Version update. Header & code cleanup to reflect
;                   changes elsewhere and to remove unused variables.
;    14-Jun-2013 DL Replace delvarx with undefine to ensure pointer
;                   removal.
;    01-Sep-2014 DL Code updated to v1.4 to reflect code cleanup and
;                   modifications elsewhere. 
;
;-

pro corpita_kernel, map_0, index_0, x_A=x_A, y_A=y_A, tilt_deg=tilt_deg, tot_angle=tot_angle, $
                    incr=incr, st_ang=st_ang, savedir=savedir, t_start=t_start, ind_lon=ind_lon, $
                    verb=verb
  
  sz = size(map_0.data)

; Get longitude and latitude of source point.
  map2wcs, map_0, wcs
  wcs_convert_from_coord, wcs, [x_A, y_A], 'hg', lon_srce, lat_srce
  latr_srce = temporary(lat_srce*!dtor) & lonr_srce = temporary(lon_srce*!dtor)
  srce = [latr_srce, lonr_srce]
  
; Get coordinates of image
  crd = wcs_get_coord(wcs)
  wcs_convert_from_coord, wcs, crd, 'hg', lon_all, lat_all
  latr = temporary(lat_all*!dtor) & lonr = temporary(lon_all*!dtor)

; Identify off-limb pixels & set them equal to NaN
  h = float(reform(sqrt((crd[0,*,*] - 0.)^2. + (crd[1,*,*] - 0.)^2.), sz[1], sz[2]))
  missing = !values.f_nan

; Remove unused variables from memory
  ;undefine, lon_all, lat_all, lon_srce, lat_srce, crd
  
; Get Great Circle distance of each pixel in the image from the source 
  d_lambda = (lonr - lonr_srce) ; Longitude Difference
  d_phi = (latr - latr_srce)    ; Latitude Difference
  
  a = (((cos(latr))*(sin(d_lambda)))^2. + (((cos(latr_srce))*(sin(latr))) - ((sin(latr_srce))*(cos(latr))*(cos(d_lambda))))^2.)
  b = ((sin(latr_srce))*(sin(latr))) + ((cos(latr_srce))*(cos(latr))*(cos(d_lambda)))	
  arctan, b, sqrt(a), grt_ang, grt_ang_deg

; Resulting array of great circle distances. 
  grt_ang_deg = reform(grt_ang_deg, sz[1], sz[2])

; Get the azimuth of each pixel in the image relative to the source.
  arctan,(cos(latr_srce)*sin(latr) - sin(latr_srce)*cos(latr)*cos(lonr-lonr_srce)), (cos(latr)*sin(lonr-lonr_srce)), az_r, az_d  
  az_d = reform(az_d, sz[1], sz[2])
  az_d[where(h ge index_0.rsun_obs)] = missing
  
; Delete unused variables
  ;undefine, a, b, grt_ang, az_r, d_lambda, d_phi, lonr, latr, lonr_srce, latr_srce

; Using 1 degree increments along arc and then a comparison across arcs.  
  rad_incr = 181                   ; No. increments in radial direction.
  angles = findgen(rad_incr)*!dtor ; Array of angles from source towards limb. 

; Define array of indices where arc reaches limb
  ind_lon = fltarr(n_elements(tilt_deg)) ; The index where the arc reaches the limb 

; Determine the corresponding arc lines.
  for i_ang = st_ang, n_elements(tilt_deg)-1, incr do begin
     

; Indices of pixels that we're interested in
     lwr_lim = i_ang - (tot_angle/2.)
     upr_lim = i_ang + (tot_angle/2.)  

     if lwr_lim lt 0. then lwr_lim = lwr_lim + 360.
     if upr_lim gt 360. then upr_lim = upr_lim - 360.

     if lwr_lim gt upr_lim then $
        a_indices = [where((az_d ge lwr_lim) and (h lt (index_0.rsun_obs-1.))), where((az_d lt upr_lim) and (h lt (index_0.rsun_obs-1.)))] $
     else $
        a_indices = where((az_d ge lwr_lim) and (az_d lt upr_lim) and (h lt (index_0.rsun_obs-1.)))

; Number of bins we're interested in here. The final array will
; have numbers up to maxbin and then NaNs.
     ind_lon[i_ang/incr] = ceil(max(grt_ang_deg[a_indices]))

; Using code from nrgf.pro by Huw Morgan to histogram out the
; pixels by great circle angle from the source.
     bin_count = histogram(grt_ang_deg[a_indices], nbin = ind_lon[i_ang/incr], min = 0, binsize = 1, reverse_ = ri, /nan) 
     ind_ok = where(bin_count gt 0, cnt_ok, comp = indnotok, ncomp = cntnotok) 
 
     if i_ang eq st_ang then begin
        cntok = cnt_ok 
        rev_ind = [ptr_new(ri, /no_copy)]
        indok = [ptr_new(ind_ok, /no_copy)]
        bincount = [ptr_new(bin_count, /no_copy)]
        arc_ind = [ptr_new(a_indices, /no_copy)]
     endif else begin
        cntok = [cntok, cnt_ok]
        rev_ind = [rev_ind, ptr_new(ri, /no_copy)]
        indok = [indok, ptr_new(ind_ok, /no_copy)]
        bincount = [bincount, ptr_new(bin_count, /no_copy)]
        arc_ind = [arc_ind, ptr_new(a_indices, /no_copy)]
     endelse

    
     
; **************************************
; Define arc sector bounds for plotting
; **************************************
     lwr_arc_ind = where((az_d ge lwr_lim) and (az_d lt lwr_lim+0.25) and (h lt (index_0.rsun_obs-1.)))
     upr_arc_ind = where((az_d le upr_lim) and (az_d gt upr_lim-0.25) and (h lt (index_0.rsun_obs-1.)))

     lwr_arc_ind_sort = lwr_arc_ind[sort(grt_ang_deg[lwr_arc_ind])]
     line_deg = round(grt_ang_deg[lwr_arc_ind_sort]*10.)/10.0
     lwr_arc_ind_sort = lwr_arc_ind_sort[uniq(line_deg)]
     line_deg = line_deg[uniq(line_deg)]

     xypix = array_indices(map_0.data, lwr_arc_ind_sort)
     str_ang = +string(i_ang, format='(I3)')
    if i_ang eq st_ang then begin
        xy_pixel_profs = CREATE_STRUCT('name', 'xy_pixel_profs', 'xypix_'+str_ang, xypix, 'pos_deg'+str_ang, line_deg) 
    endif else begin
        xy_pixel_profs = add_tag(xy_pixel_profs, xypix, 'xypix_'+str_ang)
        xy_pixel_profs = add_tag(xy_pixel_profs, line_deg, 'pos_deg'+str_ang)
    endelse    

     plots, xypix[0, *], xypix[1, *], color=140

     if i_ang gt 250 then stop
     ;prof = interpolate(map_0.data, xypix[0, *], xypix[1, *]) 
     ;window, 1, xs=900, ys=500
     ;plot, line_deg, smooth(prof, 5), /ylog, yr=[1, max(prof)], psym=3
     ;line_deg = line_deg[uniq(line_deg)]

     if i_ang eq st_ang then lwr_arc = [ptr_new(lwr_arc_ind, /no_copy)] else lwr_arc = [lwr_arc, ptr_new(lwr_arc_ind, /no_copy)]
     if i_ang eq st_ang then upr_arc = [ptr_new(upr_arc_ind, /no_copy)] else upr_arc = [upr_arc, ptr_new(upr_arc_ind, /no_copy)]

; Delete variables
     ;undefine, lwr_lim, upr_lim, cnt_ok, ri, ind_ok, bin_count, a_indices, lwr_arc_ind, upr_arc_ind

  endfor

  ;save, filename = savedir+'/CorPITA_arcs_'+time2file(t_start, /sec)+'.sav', $
  ;      lwr_arc, upr_arc, ind_lon, /compress
  ;save, filename = savedir+'/CorPITA_locs_'+time2file(t_start, /sec)+'.sav', $
  ;      missing, ind_lon, indok, cntok, rev_ind, bincount, arc_ind, wcs, /compress
  ;save, filename = savedir+'/CorPITA_pixel_vals_'+time2file(t_start, /sec)+'.sav', $
  ;      az_d, grt_ang_deg, h, /compress

  ;undefine, az_d, h, grt_ang_deg, cntok, rev_ind, bincount, arc_ind, wcs

end

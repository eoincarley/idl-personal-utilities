;+ 
; PURPOSE: 
;	Computes the mean background out ouf a defined
;	set of neighbouring profiles and subtracts it.
; CALLING SEQUENCE:
;	result = ConstBackSub( image, min, max [, direction ] )
; INPUTS:
;	image: 2D array
;	min, max: pixel number of beginning and end of the
;		background region, respectively
;	direction: 'X' or 'Y'
; KEYWORDS:
;	BACKGROUND: to get the background array subtracted
;	AUTOMATIC: the background array to subtract is
;		searched automatically by taking some vectors
;		in the image with lowest mean and lowest
;		variance as background vectors. If a number is
;		passed it is taken as percet of channels to 
;		consider else the default is 5 %.
;	BKGARR: the background field computed by the AUTOMATIC proc.
; PROCEDURE:
;	Between min and max, the values are averaged in
;	X or Y, and the resulting vector is subtracted from the
;	original image.
; SIDE EFFECT:
;	The image is returned as FLOAT type.
; MODIFICATION HISTORY:
;
;       2004-2-24: some adaptations while integrating in the objects
;                  # - stuff
;                  sig_array stuff
;                  csillag@fh-aargau.ch
;       1998-6-1:  Performance optimization, June 98, P.Messmer
;	1993-7-1:  AUTOMATIC in july 93, A.Cs
;	1994-1:    BKGARR in Jan 94, A.Cs.
;	1991:      Created in August 1991 by A.Csillaghy,
;		   Inst. of Astronomy, ETH Zurich
;-

FUNCTION ConstBackSubEdit, image, min, max, direction, $
	BACKGROUND = background, $
	AUTOMATIC = automatic, BKGARR = bkgArr, verbose = verbose, min_stdev_list=min_stdev_list



;  t = systime(1)

checkvar, direction, 'X'
yes_x =  strupcase( direction )  eq 'X'

im = float( yes_x ? transpose( image ) : image )

imsize = size( im, /struct )

nx = imsize.dimensions[0]
ny = imsize.dimensions[1]

IF Keyword_Set( AUTOMATIC ) THEN BEGIN

; First compute average over time  for every frequency channel
; Therefor the average is computed along dimension 0

    average_arr = Avg( im, 1 )

; subtract this average value from every channel

    im = im - average_arr#(fltarr(ny) + 1)

; now compute the standard-deviation for every timestep

    sdev_arr = sig_array( im, 1 )

; build the list of background candidates: take those 
; time steps with lowest standard deviation

    list = sort(sdev_arr)

; keep only a certain amount of the possible background
; candidates

    IF automatic EQ 1 THEN automatic = 0.05
    
    nPart = (ny*automatic) < ( ny-1) > 0
    list = list( 0: nPart )
    
ENDIF ELSE BEGIN
    
    min = min > 0 < (ny-1) 
    max = max > 0 < (ny-1)
    list = LIndGen( max-min ) + min
    
ENDELSE
min_stdev_list = list

RETURN, im

END

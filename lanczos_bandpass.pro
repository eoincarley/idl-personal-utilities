;----------------------------------------------------------------------------
; $Id: lanczos_bandpass.pro,v 1.5 2002/06/08 22:43:44 johnny Exp johnny $
;+
; NAME:
;   LANCZOS_BANDPASS
;
; PURPOSE:
;   Bandpass filters a timeseries using a Lanczos smoothed symmetric
;   non-recursive Fourier method filter.
;
; CATEGORY:
;   Time series analysis.
;
; CALLING SEQUENCE:
;   Result = LANCZOS_BANDPASS(in_data, in_fc1, in_fc2, in_n)
;
; INPUTS:
;   in_data:      Input timeseries to filter.  Floating/double vector
;                 dimensioned (NT).  Timestep (data interval) assumed 
;                 constant between points in vector.  Not changed by 
;                 procedure.
;
;   in_fc1:       Cut-in frequency (i.e. filter has zero response from
;                 f=[0,fc1]).  In units of cycles per data interval.
;                 Floating/double scalar.  Not changed by procedure.
;
;   in_fc2:       Cut-out frequency (i.e. filter has zero response from
;                 f=[fc2,Nyquist]).  In units of cycles per data interval.
;                 Floating/double scalar.  Not changed by procedure.
;
;   in_n:         The parameter n, which defines how many 2n+1 weights to
;                 use in the filter (so weights go from k=-n to n).
;
; KEYWORD PARAMETERS:
;   DETREND:      If keyword is set true, detrending is accomplished by
;                 removing a linear fit based on minimizing the chi-squared 
;                 error, prior to filtering.  After filtering, the trend
;                 is added back in to the filtered timeseries.
;
;   ZERO_ENDPTS:  If keywords is set true, the endpoints of the dataset
;                 are set to 0 prior to filtering, then returned to their
;                 original value in the filtered timeseries.
;
; OUTPUTS:
;   Result:       Filtered timeseries, same units and size as input time-
;                 series.  Type floating.  Created/overwritten by function.
;
; FILE DEVICE I/O:
;   None.
;
; COMMON BLOCKS:
;   None.
;
; EXAMPLE:
;   Establish a noisy dataset:
;      original = SIN((FINDGEN(2000)/400)^2.5)
;      noisy    = original+((RANDOMU(SEED,2000)-.5)/ 2)
;   Set the filter to low pass:
;      fc1 = 0.0
;      fc2 = 0.05
;      n   = 53
;   Lanczos filter the timeseries:
;      Result = LANCZOS_BANDPASS(noisy, fc1, fc2, n, /Detrend)
;   Plot the result to see how the filter removes high frequency noise.
;
; MODIFICATION HISTORY:
; - 6 Jun 2002:  Orig. ver. Johnny Lin, CIRES/University of Colorado.
;   Email:  air_jlin@yahoo.com.  Passed minimally passable tests.
; - 8 Jun 2002:  Add case so if Zero_Endpts is true, the filtered data
;   timeseries endpoints are reset to the "original" value prior to
;   filtering.  Passed minimally passable tests.
;
; NOTES:
; - Written for IDL 5.3.
; - References:
;   * Duchon, C. E. (1979):  "Lanczos filtering in one and two dimensions,"
;     J. Applied Meteorology, Vol. 18, pp. 1016-1022.
;   * Hartmann, D. L. (2002):  Atmospheric Sciences 552 Class Notes:
;     Objective Analysis, Chapter 7:  Filtering.  Seattle, WA, University
;     of Washington.  Online at:  http://www.atmos.washington.edu/~dennis/
;     552_Notes_ftp.html.
;   * Press, W. H., B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling
;     (1989):  Numerical Recipes in Pascal.  Cambridge University Press.
; - Press et al. (1989, Sec. 12.9, p. 479-480) notes when using the FFT 
;   to filter real data, the filter function H should obey H(-f)=H(f)*.  
;   This function meets that criteria, since the Lanczos filter is real 
;   and symmetric in f (even for an odd timeseries, since in that case 
;   the Nyquist frequency is omitted).
; - Press et al. (1989, Sec. 12.9, p. 480) suggests you should detrend 
;   the data prior to filtering by subtracting the line defined by the 
;   endpoints of the timeseries.  This enables the endpoints to equal 
;   zero.  In my tests, this can make a difference in the filtered 
;   values you get at the ends of the timeseries, but away from the ends 
;   this detrending makes little difference.  However, I found that 
;   retrending this line sometimes makes the filtered values a poorer 
;   fit to the dataset, since the endpoints of the data can at times 
;   have values that are very far apart from 0, creating an "artific-
;   ially" steep line for part of the timeseries.  Thus, instead of 
;   accomplishing the zeroing of the endpoints with detrending, a 
;   keyword (Zero_Endpts) is available to allow the user to set only 
;   the endpoints to zero.  Note this detrending I'm referring to here 
;   is entirely separate from the minimized chi-square error linear 
;   fit detrending keyword option (Detrend).
; - For a geophysical timeseries (like rainfall), the difference between 
;   using and not using the Zero_Endpts keyword is small (c. 0.1% of the 
;   filtered value) and limited to the points near the ends of the time-
;   series.  For the timeseries given in the example, the influence is
;   also limited to the points near the ends, but is more on the order
;   of 1-6%.
; - Nyquist frequency is +/- 0.5 cycles/data interval.
; - Keywords are optional unless otherwise noted.
; - No procedures called with _Extra keyword invoked.
; - No user-written procedures called.
;-
; Copyright (c) 2002 Johnny Lin.  For licensing, distribution conditions,
; and contact information, see http://www.johnny-lin.com/lib.html.
;----------------------------------------------------------------------------

FUNCTION LANCZOS_BANDPASS, in_data, in_fc1, in_fc2, in_n  $
                         , DETREND     = detrend  $
                         , ZERO_ENDPTS = zero_endpts  $
                         , _EXTRA      = extra




; ------------------ Initial Declarations and Error Checks ------------------

COMPILE_OPT IDL2
ON_ERROR, 0

data  = in_data               ;- protect input parameters
fc1   = in_fc1[0]
fc2   = in_fc2[0]
flt_n = FLOAT(in_n[0])           ;+ floating value of parameter n
int_n = LONG(in_n[0])            ;+ integer value of parameter n

NT     = N_ELEMENTS(data)     ;- number of elements of data (as integer
flt_NT = FLOAT(NT)            ;  and as floating type
dtype  = SIZE(data, /Type)    ;- type code of input data

if (SIZE(data, /N_Dimensions) ne 1) then  $      ;- check data is a vector
   MESSAGE, 'error--data not vector'

if (fc2 le fc1) then  $                    ;- check fc2 > fc1
   MESSAGE, 'error--fc2 must be larger'

if ( flt_n lt (1.3/(fc2-fc1)) ) then  $    ;- check n large enough to ensure
   MESSAGE, 'error--n too small'           ;  unit response at band center




; ------------ Detrend Input Data and Zero Endpoints (if chosen) ------------
;
; Key outputs of this section:  data, data_0, data_ntm1, trendline1.  
; Regarding zeroing endpoints, see Press et al. (1989, Sec. 12.9, p. 480).

if (KEYWORD_SET(detrend) eq 1) then begin   ;- remove linear trend from data
   xval = FINDGEN(NT)                          ;+ x-axis "values"

   fit_param = LINFIT(xval, data)              ;+ curve-fit parameters
   slope = fit_param[1]
   yint  = fit_param[0]

   trendline1 = (slope*xval) + yint            ;+ remove trend
   data = TEMPORARY(data) - trendline1
endif


if (KEYWORD_SET(zero_endpts) eq 1) then begin   ;- make the first and last
    data_0     = data[0]                        ;  points of the timeseries 
    data_ntm1  = data[NT-1]
    data[0]    = 0.0
    data[NT-1] = 0.0
endif




; -------------- Calculate Weight Function for Smooth Response --------------
;
; Note:  At k=0 (which exists at element int_n), the expression for wbar_k 
; (Duchon 1979, eqn. 11) must be evaluated by l'Hopital's Rule.  This is
; why k[int_n] is set to NaN for the first set of calculations of wbar_k,
; then wbar_k[int_n] is calculated separately.
;
; Key output from this section:  k, wbar_k.

k = FINDGEN((2*int_n)+1) - flt_n     ;- Variable k=[-n,n]
k[int_n] = !VALUES.F_NAN             ;- Inititalize k=0

sigma         = SIN(!PI*k/flt_n) / (!PI*k/flt_n)
wbar_k        = (  ( SIN(2.*!PI*fc2*k) / (!PI*k) )  $
                 - ( SIN(2.*!PI*fc1*k) / (!PI*k) )  )  $
              * sigma
wbar_k[int_n] = 2.0*(fc2-fc1)




; ----------------------- Calculate Frequency Vector ------------------------
;
; The frequency vector (freq) is the frequency values corresponding to the
; Fourier transform returned by the IDL FFT function, is in units cycles per
; data interval, and is the key output of this section.

freq = REPLICATE(!VALUES.F_NAN, NT)        ;- initialize freq vector
NF   = FLOOR(flt_NT/2.)                    ;- no. of pos. non-0 freq. pts.
freq[0:NF] = (FINDGEN(NT))[0:NF]  $        ;- fill f=0 and f>0 freq.
           / flt_NT

if ((NT mod 2) eq 0) then begin            ;- if NT is even, fill f<0 freq.
   freq[NF+1:*] = -REVERSE(freq[1:NF-1])
endif else begin                           ;- if NT is odd, fill f<0 freq.
   freq[NF+1:*] = -REVERSE(freq[1:NF])
endelse

tmp = WHERE(FINITE(freq) eq 0, count)                   ;- error check freq
if (count ne 0) then MESSAGE, 'error--bad freq fill'    ;  fill




; ------------------ Calculate Frequency Response Function ------------------
;
; Selected key variables:
;   tmpsum    An n-element vector which is the value of the argument inside
;             the summation symbol in eqn. 7 in Duchon (1979), for each
;             value of k=[1,n], for a freq[i].
;
; Key output from this section:  Rbar_n.

Rbar_n     = FLTARR(NT)           ;- init. smoothed freq. response fctn.
pos_k      =      k[int_n+1:*]    ;- Values of k=[1,n]
pos_wbar_k = wbar_k[int_n+1:*]    ;- Values of wbar_k corresponding to pos_k

for i=0,NT-1 do begin
    tmpsum    = pos_wbar_k * COS(2.*!PI*freq[i]*pos_k)
    Rbar_n[i] = wbar_k[int_n] + (2. * TOTAL(tmpsum))
endfor




; ------------------ Apply Filter and Retrend (if chosen) -------------------

data_filtered = FLOAT( FFT(FFT(data,1)*Rbar_n, -1) )       ;- FFT filter


if (KEYWORD_SET(zero_endpts) eq 1) then begin   ;- restore the first and last
    data_filtered[0]    = data_0                ;  points of the timeseries 
    data_filtered[NT-1] = data_ntm1
endif


if (KEYWORD_SET(detrend) eq 1) then begin                  ;- retrend
   data_filtered = TEMPORARY(data_filtered) + trendline1
endif




; ----------------------------- Prepare Output ------------------------------

Result = TEMPORARY(data_filtered)          ;- return filtered timeseries

RETURN, Result




END     ;=== end of function ===

; ========== end of file ==========
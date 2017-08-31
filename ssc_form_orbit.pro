;+
; Project     :	STEREO - SSC
;
; Name        :	SSC_FORM_ORBIT()
;
; Purpose     :	Form planetary orbit positions
;
; Category    :	STEREO, Orbit
;
; Explanation :	Calculates the orbital period of a planet or the Moon, and
;               returns the orbital positions over that period.  Used for
;               overplotting the orbit path.
;
;               The orbital positions are first calculated in an inertial
;               coordinate system, and then transformed into the desired
;               system.  This avoids the orbital distortions seen in rotating
;               coordinate systems.
;
; Syntax      :	Orbit = SSC_FORM_ORBIT(DATE, BODY, SYSTEM)
;
; Examples    :	See SSC_PLOT_WHERE
;
; Inputs      :	DATE    = Scalar date to calculate the orbit for.  The
;                         positions are calculated for +/- half the orbital
;                         period around this date.
;
;               BODY    = The name of the body to calculate the orbit for.  Can
;                         either be the name of a planet, or "Moon".
;
;               SYSTEM  = The name of the coordinate system.
;
; Opt. Inputs :	None.
;
; Outputs     :	The result of the function is a 3x361 array of coordinate
;               positions.
;
; Opt. Outputs:	None.
;
; Keywords    :	EDATE  = Date to use for addressing the ephemeris.  Used when
;                        DATE is outside the valid range of the ephemeris.
;                        Mainly used for Ulysses.
;
;               AU     = If set, then the coordinates are returned in
;                        Astronomical Units, instead of the default of
;                        kilometers.
;
; Calls       :	LOAD_STEREO_SPICE_GEN, GET_STEREO_COORD, ANYTIM2UTC, TAI2UTC,
;               UTC2TAI, CONVERT_STEREO_COORD
;
; Common      :	None.
;
; Restrictions:	None.
;
; Side effects:	Loads the generic kernels if not already loaded.
;
; Prev. Hist. :	None.
;
; History     :	Version 1, 23-Jun-2006, William Thompson, GSFC
;               Version 2, 27-Jun-2006, William Thompson, GSFC
;                       Added keyword EDATE
;               Version 3, 13-Sep-2011, WTT, better MU values
;
; Contact     :	WTHOMPSON
;-
;
function ssc_form_orbit, date, body, system, au=au, edate=edate
on_error, 2
;
;  Check the input parameters.
;
if n_params() ne 3 then message, $
  'Syntax: Orbit = SSC_FORM_ORBIT(DATE, BODY, SYSTEM)'
if n_elements(date) ne 1 then message, 'DATE must be a scalar'
;
;  Make sure that the generic kernels are loaded.
;
load_stereo_spice_gen
;
;  Set up the inertial coordinate system, and the gravitational parameter mu.
;
test = strupcase( ntrim(body) )
if (test eq 'MOON') or (test eq '301') then begin
    base = 'GEI'
    mu = 3.986004418D5          ;G * M_earth (km^3/s^2)
end else begin
    base = 'HAE'
    mu = 1.32712440018D11       ;G * M_sun
endelse
;
;  Get the complete state vector, and use it to calculate the conic
;  coefficients.
;
if n_elements(edate) eq 0 then edate = date
state = get_stereo_coord(edate, body, system=base)
utc = anytim2utc(edate,/ccsds)
cspice_utc2et, utc, et
cspice_oscelt, state, et, mu, elts
;
;  Calculate the semi-major axis A, and the orbital period P.
;
a = elts[0] / (1 - elts[1])
p = 2*!dpi * sqrt(a^3 / mu)
;
;  Set up a series of times ranging over +/- half the period, about 1 degree
;  apart.  Extrapolate the coordinates in the inertial coordinate system.
;
;
utc = tai2utc(utc2tai(utc) + p*(dindgen(361)/360 - 0.5d0), /ccsds)
coord = dblarr(3,361)
for i=0,360 do begin
    cspice_utc2et, utc[i], et
    cspice_conics, elts, et, temp
    coord[*,i] = temp[0:2]
endfor
;
;  Convert to the requested coordinate system.
;
convert_stereo_coord, date, coord, base, system
;
;  If requested, convert to Astronomical Units.
;
if keyword_set(au) then coord = coord / 1.4959787D8
;
return, coord
end
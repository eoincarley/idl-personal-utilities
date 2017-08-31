function get_fits_times, folder

current = curdir()

cd, folder

files = findfile('*.fts')
times = dblarr( n_elements(files) )

IF n_elements(files) gt 0.0 THEN BEGIN

	FOR i = 0, n_elements(files)-1 DO BEGIN
		he = LASCO_FITSHDR2STRUCT( HEADFITS(files[i]) )
		times[i]  = anytim(he.date_obs, /utim)
	ENDFOR	
	
ENDIF ELSE BEGIN
print,'Did not find any files with .fts extension'
ENDELSE


cd,current
return,times


END
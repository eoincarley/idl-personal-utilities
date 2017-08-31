pro find_rhessi_discons

	cd,'~/Downloads/'
	file='count_rhessi.txt'
	r=rd_tfile(file,/nocomment,/auto)
	t= r[0,*]
	z= r[3,*]/1.

	times = anytim('2003-01-01T'+ t[*], /utim)

	!p.multi=[0,1,2]

	utplot, times, z, $
		/xs, $
		/ys, $
		/ylog, $
		yr=[1e1, 1e4]


	utplot, times, deriv(times, z), $
		/xs, $
		/ys, $
		yr=[-500, 500];, $
		;/ylog;, $
		;yr=[1e2, 1e4]

	result = deriv(times, z)	

	positions1 = where(result ge 50.0)
	positions2 = where(result le -50.0)
	;positions2 = positions2[n_elements(positions2)-1]

			
STOP

END
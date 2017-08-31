pro normalising_gradient

	loadct, 74
	
	!p.charsize = 1
	cd, '~/Data/22Sep2011_event/herringbones'
	radio_spectro_fits_read, 'BIR_20110922_104459_01.fit', data_raw, times, freq
  
	; Now forward drift bursts. Second.
	; best performance is angle1 = 138, angle2 = 175

	t1_index = closest(times,anytim(file2time('20110922_105110'),/utim))
	t2_index = closest(times,anytim(file2time('20110922_105400'),/utim))
	f1_index = closest(freq, 43.0)
	f2_index = closest(freq, 17.0)
	inten0 = -0.4
	inten1 = 1.0
	outname = 'peak_ft_second_master_forward.sav'
	burst_file = 'bursts_ft_second_master_forward.txt'
	smooth_param = 1
	continu = 1

	data_bs = constbacksub(data_raw, /auto)	
	data_normal = data_bs
	FOR i = 0, n_elements(data_bs[0, *])-1 DO BEGIN
		data_normal[*, i] = data_bs[*, i]/max(data_bs[*, i])
	ENDFOR
	;data_normal = smooth(data_normal, 2)	;filter_image(data_bs, fwhm=2)

	window, 1, xs=2400, ys=700
	spectro_plot, bytscl(data_normal, inten0, inten1), times, freq, $
		/ys, $
		ytitle = '!6Frequency [MHz]', $
		yticks = 5, $
		yminor = 4, $
		yr = [freq[f1_index],freq[f2_index]], $
		xrange = [times[t1_index],times[t2_index]], $
		/xs, $
		xtitle = 'Start time: '+anytim(times[t1_index], /cc, /trun)+' (UT)', $
		charsize = 2.0


	;------------------------------------------------;
	;			Manually select the bursts
	;
	bf = 41
	btimes = 0.0	
	WHILE continu eq 1 DO BEGIN
		loadct, 74, /silent
		
		spectro_plot, bytscl(data_normal, inten0, inten1), times, freq, $
			/ys, $
			ytitle = '!6Frequency [MHz]', $
			yticks = 5, $
			yminor = 4, $
			yr = [freq[f1_index],freq[f2_index]], $
			xrange = [times[t1_index],times[t2_index]], $
			/xs, $
			xtitle = 'Start time: '+anytim(times[t1_index], /cc, /trun)+' (UT)', $
			charsize = 2.0

			;plots, btimes, bf, color=4, symsize=0.4, psym=8

		finished = 0
		WHILE finished ne 1 DO BEGIN 
	  			loadct, 0, /silent
	            point, btimes, bf, /data
	            finished = 1
	    ENDWHILE	
	    tindex = btimes
	    findex = bf
	    FOR i=0, n_elements(btimes)-1 DO BEGIN
            tindex[i] = closest(times, btimes[i])
            findex[i] = closest(freq, bf[i])
        ENDFOR
        inten = interpolate(data_bs, tindex, findex) 

        write_text, btimes, bf, inten, burst_file
        ;manual = 'n'
        ;READ, continu, prompt = 'Conintue? (1/0):'
    ENDWHILE

stop
END

pro write_text, bt, bff, inten, burst_filename

  IF file_test(burst_filename) eq 1 THEN BEGIN
    readcol,burst_filename, btprev, bffprev, intenprev,$
    format = 'A,D,D'
  
    bt = [btprev, '-', anytim(bt, /ccs)]
    bff = [bffprev, !Values.F_NAN, bff]
    inten = [intenprev, !Values.F_NAN, inten]

  ENDIF

  writecol, burst_filename, anytim(bt, /ccs), bff, inten, fmt='(A,D,D)'

END
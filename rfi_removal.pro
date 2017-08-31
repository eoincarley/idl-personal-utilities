function rfi_removal, raw_data

	;raw_data: the raw spectrogram
	
	;return rfi removed, background subtracted data
	
	
	bs_data = constbacksub(raw_data, /auto)

	rfi_array = abs(sobel(bs_data)) > 20.0 < 100.0  ;requires a little bit of tuning here

	indeces = where(rfi_array gt 2.0*stdev(rfi_array))
	rfi_xy = array_indices(rfi_array, indeces)

	;---------------------Do a 10 point average around every RFI point ---------------;
	raw_new = raw_data
	FOR i = 0.0, n_elements(rfi_xy[0,*])-1 DO BEGIN
		x_index = rfi_xy[0,i]
		y_index = rfi_xy[1,i]
		print,x_index, y_index
		surround = 10  ; Size of box around RFI point
		array_size = size(raw_new)
		xsize = array_size[1]
		ysize = array_size[2]
	
		section = raw_new[(x_index-surround)>0.0:(x_index+surround)<(xsize-1), (y_index-surround)>0.0:(y_index+surround)<(ysize-1)]
		raw_new[(x_index-1)>0.0:(x_index+1)<(xsize-1), (y_index-1)>0.0:(y_index+1)<(ysize-1)] = mean(section > 0.0)
	ENDFOR
bs_new = constbacksub(raw_new, /auto)
return,bs_new

END

pro remove_nans, input, output, number_positions, nan_positions, return_img=return_img

	true_false = finite(input) ;0s for NaN, 1s for number
	number_positions = where(true_false eq 1)
	nan_positions = where(true_false eq 0)
	
	if keyword_set(return_img) then begin
		xy_indeces = array_indices(input, nan_positions)
		output = input
		output[xy_indeces[0, *], xy_indeces[1, *]] = 0.0
	endif else begin
		output = input[number_positions]
	endelse	


END
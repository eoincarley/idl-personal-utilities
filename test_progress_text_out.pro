pro test_progress_text_out



	openw, 30, '~/test.txt'
	for i=0, 1000 do begin

		progress_percent, i, 0, 1000
		progress_percent_txt_output, i, 0, 1000
		wait, 0.1
	endfor
	close, 30



END
; Routine to calculate operating speed of computer

pro whosthefastest

	t1 = systime(1)
	for i=0, 9 do time_test3
	t2 = systime(1)
	print, 'Average of 10 runs :'+arr2str((t2-t1)/10)+'seconds'

end

dirsToEnter="1.OLD-RUNS-07-09-16/ ""2.NEW-RUNS-07-09-16/ ""3.OLD100-RUNS-07-09-16/ ""4.NEW100-RUNS-07-09-16/ ""5.NRs-RUNS-07-09-16/ ""6.NRs100-RUNS-07-09-16/ ""7.B0-total-RUNS-07-09-16/ ""8.B0-gg-RUNS-07-09-16/ ""9.B0-str-RUNS-07-09-16/ ""10.B0-unstr-RUNS-07-09-16/"
for i in $dirsToEnter ;
	do cd $i ;
		qsub *.pbs ;
		cd ../; 
	done
## Only run once:
## head -1 all_simulated_phenos.txt | cut -f2- | sed 's/\t/\n/g' - > traits.txt

i=1

while read trait
do
    mkdir -p GAPIT_${trait}
    echo ${trait} > GAPIT_${trait}/params.txt

	echo -e "Starting job $i \n"
    cd GAPIT_${trait}
    nohup Rscript ../run_GAPIT.R > GAPIT_${trait}.log &
    cd ../ 
    
    rem=$(($i%30))
    if [ $rem -eq 0 ]; then
	    echo -e "Waiting for job $i to finish... \n"
    	wait
    fi
    
    i=$(($i + 1))
    
    	
done < traits.txt

  
		 

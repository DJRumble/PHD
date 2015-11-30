#! /bin/csh
alias echo "echo > /dev/null" #cleans output to screen up slightly
convert
kappa
unalias echo

###################################################################
#Start of parameters
###################################################################
set REDUCED_DIR = $1 #location of raw maps, in units 1 (both in the same folder)
set home_dir = $2 #location of top level folder to output files to
set file450 = $3 #name of raw 450 micron map, without extension
set file850 = $4 #name of raw 850 map, without extension
set snr = $5 #number below which to cut low signal to noise
set type = $6 #b or t accepted, for beta or temperature map creation
set constvalue = $7 #temperature or beta to hold constant in calculations for b/T
set scalar450 = $8 #scalar to convert units 1 to Jy/pixel (ideally pW -> Jy/pixel but consistent conversions accepted)
set scalar850 = $9 
set fwhm450 = $10 #fwhm of 450micron telescope
set fwhm850 = $11
set output = `echo $file450 | sed 's/_s450//'` #creates output filename from 450 micron filename
set output = `echo $output | sed 's/_450//'`
set constoutput = `echo $constvalue | sed 's/\./_/'` #deals with non-integer constants for filename usage

###################################################################
#Start of bulk script
###################################################################

############### 1 -- Check raw data file #############

if(-e $home_dir"/map/"$output"/"$output"450850clip.sdf") then #cleans up old maps to be overwritten
    rm $home_dir"/map/"$output"/"$output"450850clip.sdf"
endif    
if("$type" == "t" && -e $home_dir"/map/"$output"450850temp"$constoutput"snr"$snr".sdf") then 
    rm $home_dir"/map/"$output"/"$output"450850temp"$constoutput"snr"$snr".sdf" #deletes only 1 map per run, by snr/type
else if("$type" == "b" && -e $home_dir"/map/"$output"450850beta"$constoutput"snr"$snr".sdf") then
    rm $home_dir"/map/"$output"/"$output"450850beta"$constoutput"snr"$snr".sdf"
endif
if (! -d $home_dir) then
    mkdir $home_dir #creates file structure if needed
endif #uses ~/450/filename; ~/850/filename; ~/map/filename to store outputs and temporary 450/850 reduction maps 
if (! -d $home_dir/450) then
    mkdir $home_dir/450
endif
if(! -d $home_dir/850) then
    mkdir $home_dir/850
endif
if(! -d $home_dir/map) then
    mkdir $home_dir/map
endif
if(! -d $home_dir/450/$file450) then
    mkdir $home_dir/450/$file450
endif
if(! -d $home_dir/850/$file850) then
    mkdir $home_dir/850/$file850
endif
if(! -d $home_dir/map/$output) then
    mkdir $home_dir/map/$output
endif
echo "Reducing "$REDUCED_DIR"/"$file450".sdf and "$REDUCED_DIR"/"$file850".sdf:"

############### 2 -- scalar to convert units 1 to Jy/pixel (ideally pW -> Jy/pixel but consistent conversions accepted)  #############
echo Converting Units: #scales to Jy/pixel
cmult in=$REDUCED_DIR/$file450.sdf scalar=$scalar450 out=$home_dir/450/$file450/s450scaled.sdf
cmult in=$REDUCED_DIR/$file850.sdf scalar=$scalar850 out=$home_dir/850/$file850/s850scaled.sdf

############## 3 -- Thresh maps  #############
echo Thresholding negative values: #removes negatives
stats $home_dir/450/$file450/s450scaled.sdf QUIET
set max450 = `parget maximum stats` #uses maximum value to force maximums to remain correct,
stats $home_dir/850/$file850/s850scaled.sdf QUIET #but force script to not prompt for a maximum value
set max850 = `parget maximum stats`
thresh in=$home_dir/450/$file450/s450scaled.sdf out=$home_dir/450/$file450/s450thresh.sdf thrlo=0 newlo=0 thrhi=$max450 newhi=$max450 QUIET
thresh in=$home_dir/850/$file850/s850scaled.sdf out=$home_dir/850/$file850/s850thresh.sdf thrlo=0 newlo=0 thrhi=$max850 newhi=$max850 QUIET

######### 4 -- Convolve maps with the beam size of the other #########
echo Convolving maps: #4"/pixel for 450 micron pre-alignment; 6"/pixel for 850micron, convolve with other telescope beam
set smooth450 = `echo "$fwhm850/4" | bc -l`
set smooth850 = `echo "$fwhm450/6" | bc -l`
gausmooth in=$home_dir/450/$file450/s450thresh.sdf fwhm=$smooth450 out=$home_dir/450/$file450/s450convolve.sdf
gausmooth in=$home_dir/850/$file850/s850thresh.sdf fwhm=$smooth850 out=$home_dir/850/$file850/s850convolve.sdf

############### 5 -- Align the 450 map with the 850 map ############
echo Aligning maps:
ndftrace $home_dir/450/$file450/s450convolve.sdf QUIET
set ndim450 = `parget ndim ndftrace` #get number of dimensions for maps to check if itermap cut or not
ndftrace $home_dir/850/$file850/s850convolve.sdf QUIET
set ndim850 = `parget ndim ndftrace`
if($ndim450 == 3 && $ndim850 == 3) then #combination of 2/3D maps to align 450 4" pixels onto 850 6" pixels
    wcsframe $home_dir/850/$file850/s850convolve.sdf sky-spectrum #removes 3rd dimension
    collapse in=$home_dir/850/$file850/s850convolve.sdf axis=3 out=$home_dir/850/$file850/s850collapse.sdf QUIET
    wcsframe $home_dir/450/$file450/s450convolve.sdf sky-spectrum
    collapse in=$home_dir/450/$file450/s450convolve.sdf axis=3 out=$home_dir/450/$file450/s450collapse.sdf QUIET
    wcsalign in=$home_dir/450/$file450/s450collapse.sdf out=$home_dir/450/$file450/s450align.sdf ref=$home_dir/850/$file850/s850collapse.sdf method=sincsinc rebin=TRUE accept QUIET
endif
if($ndim450 == 3 && $ndim850 == 2) then
    wcsframe $home_dir/450/$file450/s450convolve.sdf sky-spectrum
    collapse in=$home_dir/450/$file450/s450convolve.sdf axis=3 out=$home_dir/450/$file450/s450collapse.sdf QUIET
    wcsalign in=$home_dir/450/$file450/s450collapse.sdf out=$home_dir/450/$file450/s450align.sdf ref=$home_dir/850/$file850/s850convolve.sdf method=sincsinc rebin=TRUE accept QUIET
endif
if($ndim450 == 2 && $ndim850 == 3) then
    wcsframe $home_dir/850/$file850/s850convolve.sdf sky-spectrum
    collapse in=$home_dir/850/$file850/s850convolve.sdf axis=3 out=$home_dir/850/$file850/s850collapse.sdf QUIET
    wcsalign in=$home_dir/450/$file450/s450convolve.sdf out=$home_dir/450/$file450/s450align.sdf ref=$home_dir/850/$file850/s850collapse.sdf method=sincsinc rebin=TRUE accept QUIET
endif 
if($ndim450 == 2 && $ndim850 == 2) then
    wcsalign in=$home_dir/450/$file450/s450convolve.sdf out=$home_dir/450/$file450/s450align.sdf ref=$home_dir/850/$file850/s850convolve.sdf method=sincsinc rebin=TRUE accept QUIET
endif

############### 6 -- Make Ratio maps ###########################
ndftrace $home_dir/450/$file450/s450align.sdf QUIET
set error450 = `parget variance ndftrace` #checks for variance array for maps to allow snr cuts
if($ndim850 == 2) then
    ndftrace $home_dir/850/$file850/s850convolve.sdf QUIET
    set error850 = `parget variance ndftrace`
    if($error450 == "TRUE" && $error850 == "TRUE") then #cuts snr if possible
	echo "Clipping SNR below "$snr":"
	errclip in=$home_dir/450/$file450/s450align.sdf out=$home_dir/450/$file450/s450clip.sdf mode=snr limit=$snr QUIET
	errclip in=$home_dir/850/$file850/s850convolve.sdf out=$home_dir/850/$file850/s850clip.sdf mode=snr limit=$snr QUIET
	echo Creating ratio map:
	div in1=$home_dir/450/$file450/s450clip.sdf in2=$home_dir/850/$file850/s850clip.sdf out=$home_dir/map/$output/$output"450850.sdf"
	#source ./linplot_ratio $home_dir/map/$output $output"450850" $home_dir/450/$file450 $home_dir/850/$file850 s450clip s850clip
    else if($error450 == "FALSE" || $error850 == "FALSE") then #if either map doesn't have variance, no cuts done
	echo Creating ratio map:
        div in1=$home_dir/450/$file450/s450align.sdf in2=$home_dir/850/$file850/s850convolve.sdf out=$home_dir/map/$output/$output"450850.sdf" #uses correct 850 map based on whether collapse was used or not
	#source ./linplot_ratio $home_dir/map/$output $output"450850" $home_dir/450/$file450 $home_dir/850/$file850 s450align s850convolve
    endif
else if($ndim850 == 3) then #same code, just uses collapse rather than convolve due to 850 map needing to be collapsed first
    ndftrace $home_dir/850/$file850/s850collapse.sdf QUIET
    set error850 = `parget variance ndftrace`
    if($error450 == "TRUE" && $error850 == "TRUE") then
	echo "Clipping SNR below "$snr":"
	errclip in=$home_dir/450/$file450/s450align.sdf out=$home_dir/450/$file450/s450clip.sdf mode=snr limit=$snr QUIET
	errclip in=$home_dir/850/$file850/s850collapse.sdf out=$home_dir/850/$file850/s850clip.sdf mode=snr limit=$snr QUIET
	echo Creating ratio map:
	div in1=$home_dir/450/$file450/s450clip.sdf in2=$home_dir/850/$file850/s850clip.sdf out=$home_dir/map/$output/$output"450850.sdf"
	#source ./linplot_ratio $home_dir/map/$output $output"450850" $home_dir/450/$file450 $home_dir/850/$file850 s450clip s850clip
    else if($error450 == "FALSE" || $error850 == "FALSE") then
	echo Creating ratio map:
	div in1=$home_dir/450/$file450/s450align.sdf in2=$home_dir/850/$file850/s850collapse.sdf out=$home_dir/map/$output/$output"450850.sdf"
	#source ./linplot_ratio $home_dir/map/$output $output"450850" $home_dir/450/$file450 $home_dir/850/$file850 s450align s850collapse
    endif
endif

######### 8 -- Create FITS files and pass to 'temperature.py' ###########

if($error450 == "TRUE" && $error850 == "TRUE") then
    echo Cutting high variance ratio pixels: #cuts top 3% of data if possible, to remove outliers, rough edges etc
    stats $home_dir/map/$output/$output"450850.sdf" order percentiles=97 comp=variance QUIET
    set var99 = `parget perval stats`
    errclip in=$home_dir/map/$output/$output"450850.sdf" out=$home_dir"/map/"$output"/"$output"450850clip.sdf" mode=variance limit=$var99 QUIET
    rm $home_dir/map/$output/$output"450850.sdf" #deletes old map to save storage
    #source ./linplot_ratio $home_dir/map/$output $output"450850clip" $home_dir/450/$file450 $home_dir/850/$file850 s450clip s850clip
endif
if($error450 == "TRUE" && $error850 == "TRUE") then
    if ("$type" == "b") then #creates beta map with constant temperature as specified
	echo "Creating beta map with temperature = "$constvalue"K:"
	maths "(log(IA*(exp(31.97/"$constvalue")-1)/(exp(16.93/"$constvalue")-1))/log(1.888857))-3" ia=$home_dir"/map/"$output"/"$output"450850clip.sdf" out=$home_dir"/map/"$output"/"$output"450850beta"$constoutput"snr"$snr".sdf"
	echo "Created "$home_dir"/map/"$output"/"$output"450850beta"$constoutput"snr"$snr".sdf"
    else if ("$type" == "t") then
	echo "Creating temperature map with ß = "$constvalue":" #creates FITS file to output to python
	ndf2fits in=$home_dir"/map/"$output"/"$output"450850clip.sdf" out=$home_dir"/map/"$output"/*" QUIET
	ndftrace $home_dir/map/$output/$output"450850clip.sdf" QUIET
	set dim = `parget dims ndftrace` #gets dimensions for FOR loops
	python temperature.py $home_dir/map/$output $output"450850clip" $dim $constvalue yesvar #passes dimensions, constant and a flag to use variance since it has it
	fits2ndf in=$home_dir"/map/"$output"/"$output"450850clip.fit" out=$home_dir"/map/"$output"/"$output"450850temp"$constoutput"snr"$snr".sdf" QUIET #makes NDF from resulting temperature FITS
	rm $home_dir"/map/"$output"/"$output"450850clip.fit"
	echo "Cutting high variance temperature pixels:" #cuts anything with deltaT > 30K
	errclip in=$home_dir/map/$output/$output"450850temp"$constoutput"snr"$snr".sdf" out=$home_dir/map/$output/$output"450850temp"$constoutput"clipsnr"$snr".sdf" mode=variance limit=30 QUIET
	thresh in=$home_dir/map/$output/$output"450850temp"$constoutput"clipsnr"$snr".sdf" out=$home_dir/map/$output/$output"450850temp"$constoutput"snr"$snr".sdf" thrlo=0 newlo=0 thrhi=999.98 newhi=bad QUIET #cuts anything calculated > 1000K, as unphysical and breaks variance calculations
	rm $home_dir/map/$output/$output"450850temp"$constoutput"clipsnr"$snr".sdf"
	echo "Created "$home_dir"/map/"$output"/"$output"450850temp"$constoutput"snr"$snr".sdf"
    endif
endif
if($error450 == "FALSE" || $error850 == "FALSE") then
    if ("$type" == "b") then #same code as before
	echo "Creating beta map with temperature = "$constvalue"K:"
	maths "(log(IA*(exp(31.97/"$constvalue")-1)/(exp(16.93/"$constvalue")-1))/log(1.888857))-3" ia=$home_dir"/map/"$output"/"$output"450850.sdf" out=$home_dir"/map/"$output"/"$output"450850beta"$constoutput"snr"$snr".sdf"
	echo "Created "$home_dir"/map/"$output"/"$output"450850beta"$constoutput"snr"$snr".sdf"
    else if ("$type" == "t") then
	echo "Creating temperature map:" #same as before, with flag to show no variance array
	ndf2fits in=$home_dir"/map/"$output"/"$output"450850.sdf" out=$home_dir"/map/"$output"/*" QUIET
	ndftrace $home_dir/map/$output/$output"450850.sdf" QUIET
	set dim = `parget dims ndftrace`
	python temperature.py $home_dir/map/$output $output"450850" $dim $constvalue novar
	fits2ndf in=$home_dir"/map/"$output"/"$output"450850.fit" out=$home_dir"/map/"$output"/"$output"450850temp"$constoutput"snr"$snr".sdf" QUIET 
	rm $home_dir"/map/"$output"/"$output"450850.fit" #no error analysis without variance array
	echo "Created "$home_dir"/map/"$output"/"$output"450850temp"$constoutput"snr"$snr".sdf"
    endif
endif


rm -r $home_dir/450 #tidies up temporary files to free up storage space
rm -r $home_dir/850

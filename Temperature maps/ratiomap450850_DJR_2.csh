#! /bin/csh
alias echo "echo > /dev/null" #cleans output to screen up slightly
convert
kappa
unalias echo

###################################################################
#Start of parameters
###################################################################

set REDUCED_DIR = $1 #/scratch/damian/run7 #location of raw maps, in units 1 (both in the same folder)
set home_dir = $2 #location of top level folder to output files to
set file450 = $3 #serpensmwc297_s450_ext_mask2_cal491_mos_th #name of raw 450 micron map, without extension
set file850 = $4 #serpensmwc297_s850_ext_mask2_cal537_mos_th #name of raw 850 map, without extension
set snr = $5 #number below which to cut low signal to noise
set type = $6 #b or t accepted, for beta or temperature map creation
set constvalue = $7 #temperature or beta to hold constant in calculations for b/T
set map = $8 #name of map being used to select appropriate STDV values
set output = `echo $file450 | sed 's/_s450//'` #creates output filename from 450 micron filename
set output = `echo $output | sed 's/_450//'`
set constoutput = `echo $constvalue | sed 's/\./_/'` #deals with non-integer constants for filename usage

###################################################################
#Start of bulk script
###################################################################


############### 1 -- Check raw data file #############
#cleans up old maps to be overwritten
echo 1

if(-e $home_dir"/map/"$output"/"$output"450850clip.sdf") then 
    rm $home_dir"/map/"$output"/"$output"450850clip.sdf"
endif    

#deletes only 1 map per run, by snr/type
if("$type" == "t" && -e $home_dir"/map/"$output"450850temp"$constoutput"snr"$snr".sdf") then 
    rm $home_dir"/map/"$output"/"$output"450850temp"$constoutput"snr"$snr".sdf"
else if("$type" == "b" && -e $home_dir"/map/"$output"450850beta"$constoutput"snr"$snr".sdf") then
    rm $home_dir"/map/"$output"/"$output"450850beta"$constoutput"snr"$snr".sdf"
endif

#creates file structure if needed
if (! -d $home_dir) then
    mkdir $home_dir 
endif  

#uses ~/450/filename; ~/850/filename; ~/map/filename to store outputs and temporary 450/850 reduction maps
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

############### 2 -- Make SNR cut to data based on STDV #############
echo 2
#Check the map and identify the correct STDV to use in S/N cut.

#scalar to convert units 1 to Jy/pixel (ideally pW -> Jy/pixel but consistent conversions accepted)
#set scalar450 = 491
#set scalar850 = 537

#or if i have been using the wrong factors all this time and realy should have been using FCFarcsec/FCFbeam x pixel area
set scalar450 = 0.1535
set scalar850 = 0.1569

#STDV of maps
set STDV_south_450 = 0.0724853408915 #mJy
set STDV_south_850 = 0.0202817776445 #mJy
set STDV_main_450 = 0.0524739027056 #mJy
set STDV_main_850 = 0.0123515747858 #mJy
set STDV_mwc297_450 = 0.102735175629 #mJy
set STDV_mwc297_850 =  0.0163324719417 #mJy
set STDV_E1_450 = 0.107306686051 #mJy
set STDV_E1_850 = 0.0139176325298 #mJy

#Prints to terminal
echo "Clipping SNR below "$snr":"

if($map == south) then
    set sigma450 = `calc $STDV_south_450\*$scalar450\*$snr`
    set sigma850 = `calc $STDV_south_850\*$scalar850\*$snr`
else if ($map == main) then
    set sigma450 = `calc $STDV_main_450\*$scalar450\*$snr`
    set sigma850 = `calc $STDV_main_850\*$scalar850\*$snr`
else if ($map == mwc297) then
    set sigma450 = `calc $STDV_mwc297_450\*$scalar450\*$snr`
    set sigma850 = `calc $STDV_mwc297_850\*$scalar850\*$snr`
else if ($map == E1) then
    set sigma450 = `calc $STDV_E1_450\*$scalar450\*$snr`
    set sigma850 = `calc $STDV_E1_850\*$scalar850\*$snr`
else
    set sigma450 = 1
    set sigma850 = 1
endif

#Set N (the number of convolving pixels in the function)
#set N_450 = 6.809
#set N_850 = 6.663
#and becasue i can't find out how to do powers n cshell here are the SQRT of the above numbers
set N_450 = 2.609
set N_850 = 2.581

#Calculate the mean STDV by sigma/SQRT(N)
set Msigma450 = `calc $sigma450/$N_450`
set Msigma850 = `calc $sigma850/$N_850`

echo $Msigma450
echo $Msigma850
echo $map


############### 3 -- Convert to map units wth FCF  #############
echo 3
#scales maps with parameter conversion factors
#scales to Jy/pixel
echo Converting Units: 
cmult in=$REDUCED_DIR/$file450.sdf scalar=$scalar450 out=$home_dir/450/$file450/s450scaled.sdf
cmult in=$REDUCED_DIR/$file850.sdf scalar=$scalar850 out=$home_dir/850/$file850/s850scaled.sdf

######### 4 -- Convolve maps with the beam size of the other #########
echo 4
#Convolution of threshed maps
#4"/pixel for 450 micron pre-alignment; 6"/pixel for 850micron, convolve with other telescope beam

#Set FWHM for SCUBA-2 maps - MB (Main Beam) S (Secondary Beam)
set fwhm450_MB = 7.9
set fwhm850_MB = 13.0
set fwhm450_S = 25.0
set fwhm850_S = 48.0

#set normalisation constants
set alpha450 = 0.94
set alpha850 = 0.98
set beta450 = 0.06
set beta850 = 0.02

#Calculate 
#NOTE: MB only values = 8.5, 13.5, combined = 8.3, 14.6, Dempsey13 = 9.8, 14.6

set fwhm450 = 8.2949 # `calc $alpha450\*$fwhm450_MB+$beta450\*$fwhm450_S`
set fwhm850 = 14.5499 # `calc $alpha850\*$fwhm850_MB+$beta850\*$fwhm850_S`

echo Convolving maps: 
#set smooth450 = `echo "$fwhm850/4" | bc -l`
#set smooth850 = `echo "$fwhm450/6" | bc -l`

#gausmooth in=$home_dir/450/$file450/s450thresh.sdf fwhm=$smooth450 out=$home_dir/450/$file450/s450convolve.sdf
#gausmooth in=$home_dir/850/$file850/s850thresh.sdf fwhm=$smooth850 out=$home_dir/850/$file850/s850convolve.sdf

#Main Beam#
echo main beam convolution:

set fwhm450_MB = 7.9
set fwhm850_MB = 13.0

set smooth450_MB = `echo "$fwhm850_MB/4" | bc -l`
set smooth850_MB = `echo "$fwhm450_MB/6" | bc -l`

gausmooth in=$home_dir/450/$file450/s450scaled.sdf fwhm=$smooth450_MB out=$home_dir/450/$file450/s450convolve_MB.sdf
gausmooth in=$home_dir/850/$file850/s850scaled.sdf fwhm=$smooth850_MB out=$home_dir/850/$file850/s850convolve_MB.sdf

#Secondary Beam#
echo secondary beam convolution:

set fwhm450_S = 25.0
set fwhm850_S = 48.0

set smooth450_S = `echo "$fwhm850_S/4" | bc -l`
set smooth850_S = `echo "$fwhm450_S/6" | bc -l`

gausmooth in=$home_dir/450/$file450/s450scaled.sdf fwhm=$smooth450_S out=$home_dir/450/$file450/s450convolve_S.sdf
gausmooth in=$home_dir/850/$file850/s850scaled.sdf fwhm=$smooth850_S out=$home_dir/850/$file850/s850convolve_S.sdf

#Normilisation#
echo normalisation:

cmult in=$home_dir/450/$file450/s450convolve_MB.sdf scalar=$alpha450 out=$home_dir/450/$file450/s450norm_MB.sdf
cmult in=$home_dir/450/$file450/s450convolve_S.sdf scalar=$beta450 out=$home_dir/450/$file450/s450norm_S.sdf
cmult in=$home_dir/850/$file850/s850convolve_MB.sdf scalar=$alpha850 out=$home_dir/850/$file850/s850norm_MB.sdf
cmult in=$home_dir/850/$file850/s850convolve_S.sdf scalar=$beta850 out=$home_dir/850/$file850/s850norm_S.sdf

add in1=$home_dir/450/$file450/s450norm_MB.sdf in2=$home_dir/450/$file450/s450norm_S.sdf out=$home_dir/450/$file450/s450convolve.sdf
add in1=$home_dir/850/$file850/s850norm_MB.sdf in2=$home_dir/850/$file850/s850norm_S.sdf out=$home_dir/850/$file850/s850convolve.sdf


############### 5 -- Thresh maps  #############
echo 5
#Threshes maps to remove negative components from scaled maps
#Also clips maps based on SNR 
echo Thresholding negative values: #removes negatives
stats $home_dir/450/$file450/s450convolve.sdf QUIET
set max450 = `parget maximum stats` #uses maximum value to force maximums to remain correct,
stats $home_dir/850/$file850/s850convolve.sdf QUIET #but force script to not prompt for a maximum value
set max850 = `parget maximum stats`
thresh in=$home_dir/450/$file450/s450convolve.sdf out=$home_dir/450/$file450/s450thresh.sdf thrlo=$Msigma450 newlo=bad thrhi=$max450 newhi=$max450 QUIET
thresh in=$home_dir/850/$file850/s850convolve.sdf out=$home_dir/850/$file850/s850thresh.sdf thrlo=$Msigma850 newlo=bad thrhi=$max850 newhi=$max850 QUIET


############### 6 -- Align the 450 map with the 850 map ############
echo 6
#Align 450 onto 850 convolved maps (in 2 or 3D where necessary)
echo Aligning maps:
ndftrace $home_dir/450/$file450/s450thresh.sdf QUIET
set ndim450 = `parget ndim ndftrace` #get number of dimensions for maps to check if itermap cut or not
ndftrace $home_dir/850/$file850/s850thresh.sdf QUIET
set ndim850 = `parget ndim ndftrace`
if($ndim450 == 3 && $ndim850 == 3) then #combination of 2/3D maps to align 450 4" pixels onto 850 6" pixels
    wcsframe $home_dir/850/$file850/s850thresh.sdf sky-spectrum #removes 3rd dimension
    collapse in=$home_dir/850/$file850/s850thresh.sdf axis=3 out=$home_dir/850/$file850/s850collapse.sdf QUIET
    wcsframe $home_dir/450/$file450/s450thresh.sdf sky-spectrum
    collapse in=$home_dir/450/$file450/s450thresh.sdf axis=3 out=$home_dir/450/$file450/s450collapse.sdf QUIET
    wcsalign in=$home_dir/450/$file450/s450collapse.sdf out=$home_dir/450/$file450/s450align.sdf ref=$home_dir/850/$file850/s850collapse.sdf method=sincsinc rebin=TRUE accept QUIET
endif
if($ndim450 == 3 && $ndim850 == 2) then
    wcsframe $home_dir/450/$file450/s450thresh.sdf sky-spectrum
    collapse in=$home_dir/450/$file450/s450thresh.sdf axis=3 out=$home_dir/450/$file450/s450collapse.sdf QUIET
    wcsalign in=$home_dir/450/$file450/s450collapse.sdf out=$home_dir/450/$file450/s450align.sdf ref=$home_dir/850/$file850/s850thresh.sdf method=sincsinc rebin=TRUE accept QUIET
endif
if($ndim450 == 2 && $ndim850 == 3) then
    wcsframe $home_dir/850/$file850/s850thresh.sdf sky-spectrum
    collapse in=$home_dir/850/$file850/s850thresh.sdf axis=3 out=$home_dir/850/$file850/s850collapse.sdf QUIET
    wcsalign in=$home_dir/450/$file450/s450thresh.sdf out=$home_dir/450/$file450/s450align.sdf ref=$home_dir/850/$file850/s850collapse.sdf method=sincsinc rebin=TRUE accept QUIET
endif 
if($ndim450 == 2 && $ndim850 == 2) then
    wcsalign in=$home_dir/450/$file450/s450thresh.sdf out=$home_dir/450/$file450/s450align.sdf ref=$home_dir/850/$file850/s850thresh.sdf method=sincsinc rebin=TRUE accept QUIET
endif


############### 7 -- Make Ratio maps ###########################
echo 7

#checks for variance array for maps to allow snr cuts

### IF statement on when to make cuts (no. of dimensions and wether there is varience) ###
ndftrace $home_dir/450/$file450/s450align.sdf QUIET
set error450 = `parget variance ndftrace` 
if($ndim850 == 2) then
    ndftrace $home_dir/850/$file850/s850thresh.sdf QUIET
    set error850 = `parget variance ndftrace`
    if($error450 == "TRUE" && $error850 == "TRUE") then #cuts snr if possible

	echo Creating ratio map:
	div in1=$home_dir/450/$file450/s450align.sdf in2=$home_dir/850/$file850/s850thresh.sdf out=$home_dir/map/$output/$output"450850.sdf"
    #if either map doesn't have variance, no cuts done
    else if($error450 == "FALSE" || $error850 == "FALSE") then 
	echo Creating ratio map:
    #uses correct 850 map based on whether collapse was used or not
        div in1=$home_dir/450/$file450/s450align.sdf in2=$home_dir/850/$file850/s850thresh.sdf out=$home_dir/map/$output/$output"450850.sdf"
 endif

    #same code, just uses collapse rather than convolve due to 850 map needing to be collapsed first
else if($ndim850 == 3) then
    ndftrace $home_dir/850/$file850/s850collapse.sdf QUIET
    set error850 = `parget variance ndftrace`
    if($error450 == "TRUE" && $error850 == "TRUE") then
	echo "Clipping SNR below "$snr":"

	echo Creating ratio map:
	div in1=$home_dir/450/$file450/s450align.sdf in2=$home_dir/850/$file850/s850collapse.sdf out=$home_dir/map/$output/$output"450850.sdf"
    else if($error450 == "FALSE" || $error850 == "FALSE") then
	echo Creating ratio map:
	div in1=$home_dir/450/$file450/s450align.sdf in2=$home_dir/850/$file850/s850collapse.sdf out=$home_dir/map/$output/$output"450850.sdf"
    endif
endif


if($error450 == "TRUE" && $error850 == "TRUE") then
    #cuts top 3% of data if possible, to remove outliers, rough edges etc
    echo Cutting high variance ratio pixels: 
    stats $home_dir/map/$output/$output"450850.sdf" order percentiles=97 comp=variance QUIET
    set var99 = `parget perval stats`
    errclip in=$home_dir/map/$output/$output"450850.sdf" out=$home_dir"/map/"$output"/"$output"450850clip.sdf" mode=variance limit=$var99 QUIET
    #deletes old map to save storage
    #rm $home_dir/map/$output/$output"450850.sdf"
endif

######### 8 -- Create FITS files and pass to 'temperature.py' ###########
echo 8
if($error450 == "TRUE" && $error850 == "TRUE") then
    #creates beta map with constant temperature as specified
    echo path TRUE TRUE
    if ("$type" == "b") then 
	echo "Creating beta map with temperature = "$constvalue"K:"
	maths "(log(IA*(exp(31.97/"$constvalue")-1)/(exp(16.93/"$constvalue")-1))/log(1.888857))-3" ia=$home_dir"/map/"$output"/"$output"450850clip.sdf" out=$home_dir"/map/"$output"/"$output"450850beta"$constoutput"snr"$snr"_2.sdf"
	echo "Created "$home_dir"/map/"$output"/"$output"450850beta"$constoutput"snr"$snr"_2.sdf"

    else if ("$type" == "t") then
	#creates FITS file to output to python
	echo "Creating temperature map with ß = "$constvalue":"
	ndf2fits in=$home_dir"/map/"$output"/"$output"450850clip.sdf" out=$home_dir"/map/"$output"/*" QUIET
	ndftrace $home_dir/map/$output/$output"450850clip.sdf" QUIET
	set dim = `parget dims ndftrace` #gets dimensions for FOR loops

	###### SEND TO temperature.py ######
	#passes dimensions, constant and a flag to use variance since it has it
	python temperature.py $home_dir/map/$output $output"450850clip" $dim $constvalue yesvar 
	#makes NDF from resulting temperature FITS
	fits2ndf in=$home_dir"/map/"$output"/"$output"450850clip.fit" out=$home_dir"/map/"$output"/"$output"450850temp"$constoutput"snr"$snr"_2.sdf" QUIET 
	rm $home_dir"/map/"$output"/"$output"450850clip.fit"
	#cuts anything with deltaT > 30K
	echo "Cutting high variance temperature pixels:"
	errclip in=$home_dir/map/$output/$output"450850temp"$constoutput"snr"$snr"_2.sdf" out=$home_dir/map/$output/$output"450850temp"$constoutput"clipsnr"$snr"_2.sdf" mode=variance limit=30 QUIET
	#cuts anything calculated > 1000K, as unphysical and breaks variance calculations
	thresh in=$home_dir/map/$output/$output"450850temp"$constoutput"clipsnr"$snr"_2.sdf" out=$home_dir/map/$output/$output"450850temp"$constoutput"snr"$snr"_2.sdf" thrlo=0 newlo=0 thrhi=999.98 newhi=bad QUIET 
	rm $home_dir/map/$output/$output"450850temp"$constoutput"clipsnr"$snr"_2.sdf"
	echo "Created "$home_dir"/map/"$output"/"$output"450850temp"$constoutput"snr"$snr"_2.sdf"
    endif
endif


if($error450 == "FALSE" || $error850 == "FALSE") then
    echo path FALSE FALSE
    if ("$type" == "b") then #same code as before
	echo "Creating beta map with temperature = "$constvalue"K:"
	maths "(log(IA*(exp(31.97/"$constvalue")-1)/(exp(16.93/"$constvalue")-1))/log(1.888857))-3" ia=$home_dir"/map/"$output"/"$output"450850.sdf" out=$home_dir"/map/"$output"/"$output"450850beta"$constoutput"snr"$snr"_2.sdf"
	echo "Created "$home_dir"/map/"$output"/"$output"450850beta"$constoutput"snr"$snr"_2.sdf"

    else if ("$type" == "t") then
	#same as before, with flag to show no variance array
	echo "Creating temperature map:" 
	ndf2fits in=$home_dir"/map/"$output"/"$output"450850.sdf" out=$home_dir"/map/"$output"/*" QUIET
	ndftrace $home_dir/map/$output/$output"450850.sdf" QUIET
	set dim = `parget dims ndftrace`
	###### SEND TO temperature.py ######
	python temperature.py $home_dir/map/$output $output"450850" $dim $constvalue novar
	fits2ndf in=$home_dir"/map/"$output"/"$output"450850.fit" out=$home_dir"/map/"$output"/"$output"450850temp"$constoutput"snr"$snr"_2.sdf" QUIET 
	#no error analysis without variance array
	rm $home_dir"/map/"$output"/"$output"450850.fit" 
	echo "Created "$home_dir"/map/"$output"/"$output"450850temp"$constoutput"snr"$snr"_2.sdf"
    endif
endif

#tidies up temporary files to free up storage space#
#rm -r $home_dir/450
#rm -r $home_dir/850

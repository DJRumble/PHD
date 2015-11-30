#DJR UoE
#27/01/2014 - ratio450850_8.csh
#
#This code produces temperature maps
#
#it is based on TJW original code but contians the following improvements:
#1) Applying SNR cut of data
#2) Main and Secondary beam convolution
#3) swapping thresh and convolution processes around (and adjusting for these ch#anges)
#4) moves wcsalign forward infront of thresh (but not convolve)
#5) uses Flux per Pixel maps (remove process 2) and updates RMS numbers
#6) normalisation corrections to convolution and corrections to RMS numbers
#7) Now produces error and percentage error maps in both error and temp. maps in parrallel
#8) Implements Error array into the clipping process instead of the original Variance array
#9) Allows user to input what proportion of high error pixels to remove  
#10) Hard code in fixed parameters, mathethmatics and pixel extraction. 

#Output files go to directory output_7
##################################################################

#! /bin/csh
alias echo "echo > /dev/null" #cleans output to screen up slightly
convert
kappa
unalias echo

###################################################################
#Start of free parameters
###################################################################

set REDUCED_DIR = $1 #location of raw maps, in units 1 (both in the same folder)
set home_dir = $2 #location of top level folder to output files to
set file450 = $3 #name of raw 450 micron map, without extension
set file850 = $4  #name of raw 850 map, without extension
#set snr = $5 #number below which to cut low signal to noise 
set snr = 5  ######## FIXED AT 5 ########
#set type = $6 #b or t accepted, for beta or temperature map creation 
set type = t ######## FIXED AT t ########
#set constvalue = $7 #temperature or beta to hold constant in calculations for b/T
set constvalue = 1.8 ######## FIXED AT 1.8 ########
set percent = $5 #Pecentage error in temp. maps above which are cut (e.g. 3)
set map = $6 #name of map being used to select appropriate STDV values
set output = `echo $file450 | sed 's/_s450//'` #creates output filename from 450 micron filename
set output = `echo $output | sed 's/_450//'`
set constoutput = `echo $constvalue | sed 's/\./_/'` #deals with non-integer constants for filename usage

###################################################################
#Start of fixed parameters
###################################################################

alias acalc '       awk "BEGIN{ print \!* }" '

set pi = 3.14159265359

#Set FWHM for SCUBA-2 maps - MB (Main Beam) S (Secondary Beam)
set fwhm450_MB = 7.9
set fwhm850_MB = 13.0 
set fwhm450_S = 25
set fwhm850_S = 48

#set normalisation constants
set alpha450 = 0.94
set alpha850 = 0.98
set beta450 = 0.06
set beta850 = 0.02

# A (and B) = 1/(2*pi*(a4*(FM4**2)+b4*(FS4**2))) or 1/(2*pi*(a8*(FM8**2)+b8*(FS8**2)))
set A = `acalc 1/(2*$pi*($alpha450*($fwhm450_MB^2)+$beta450*($fwhm450_S^2)))`
set B = `acalc 1/(2*$pi*($alpha850*($fwhm850_MB^2)+$beta850*($fwhm850_S^2)))`

#set numbers for multiplication with convovle maps - in general they are of the form 2pi*sig^2*A*alpha
set A_a = `acalc ($fwhm450_MB^2)*$A*$alpha450*2*$pi`
set A_b = `acalc ($fwhm450_S^2)*$A*$beta450*2*$pi`
set B_a = `acalc($fwhm850_MB^2)*$B*$alpha850*2*$pi`
set B_b = `acalc($fwhm850_S^2)*$B*$beta850*2*$pi`

#Prints to terminal
echo $map

#STDV of maps - updated 17 june 2013 - http://wiki.astro.ex.ac.uk/bin/view/JCMTGouldBelt/SCUBA-2IR1NoiseEstimates

if($map == south) then
    #set STDV_south_450 = 0.00727832102989  #Jy/pixel
    set STDV_south_450 = 0.00847368513332   #Jy/pixel skyloop
    #set STDV_south_850 = 0.00141253836219  #Jy/pixel
    set STDV_south_850 = 0.00151656703904 #Jy/pixel skyloop
    set sigma450 = `calc $STDV_south_450\*$snr`
    set sigma850 = `calc $STDV_south_850\*$snr`
    set var450 = `acalc ($STDV_south_450 ^ 2.0)*(($alpha450 ^ 2.0) + ($beta450 ^ 2.0))*(0.66666666666666) `
    set var850 = `acalc ($STDV_south_850 ^ 2.0)*(($alpha850 ^ 2.0) + ($beta850 ^ 2.0))*(0.66666666666666) `
else if ($map == main) then
    set STDV_main_450 = 0.0481830581482        #0.00837840813063  #Jy/pixel
    set STDV_main_850 = 0.0101566273883       #0.00167943638908  #Jy/pixel
    set STDV_main_850noco = 0.00168595226123  #Jy/pixel
    set sigma450 = `calc $STDV_main_450\*$snr`
    set sigma850 = `calc $STDV_main_850\*$snr`
    set var450 = `acalc ($STDV_main_450 ^ 2.0)*(($alpha450 ^ 2.0) + ($beta450 ^ 2.0))*(0.66666666666666) `
    set var850 = `acalc ($STDV_main_850 ^ 2.0)*(($alpha850 ^ 2.0) + ($beta850 ^ 2.0))*(0.66666666666666) `
else if ($map == mwc297) then
    #set STDV_mwc297_450 = 0.0209323405711  #Jy/pixel --  DJR
    set STDV_mwc297_450 = 0.0164665781382 #Jy/pixel -- IR1
    #set STDV_mwc297_850 = 0.00261840219553 #Jy/pixel -- DJR
    set STDV_mwc297_850 = 0.00218410360456 #Jy/pixel -- IR1
    set sigma450 = `calc $STDV_mwc297_450\*$snr`
    set sigma850 = `calc $STDV_mwc297_850\*$snr`
    set var450 = `acalc ($STDV_mwc297_450 ^ 2.0)*(($alpha450 ^ 2.0) + ($beta450 ^ 2.0))*(0.66666666666666) `
    set var850 = `acalc ($STDV_mwc297_850 ^ 2.0)*(($alpha850 ^ 2.0) + ($beta850 ^ 2.0))*(0.66666666666666) `
else if ($map == E1) then
    set STDV_E1_450 = 0.0185238342182  #Jy/pixel
    set STDV_E1_850 = 0.00202344639093  #Jy/pixel
    set sigma450 = `calc $STDV_E1_450\*$snr`
    set sigma850 = `calc $STDV_E1_850\*$snr`
    set var450 = `acalc ($STDV_E1_450 ^ 2.0)*(($alpha450 ^ 2.0) + ($beta450 ^ 2.0))*(0.66666666666666) `
    set var850 = `acalc ($STDV_E1_850 ^ 2.0)*(($alpha850 ^ 2.0) + ($beta850 ^ 2.0))*(0.66666666666666) `
else if ($map == ophA) then ###### not working #####
    set STDV_ophA_450 = 0.0131259511768  #Jy/pixel
    set STDV_ophA_850 = 0.0187261 #Jy/pixel
    set sigma450 = `calc $STDV_ophA_450\*$snr`
    set sigma850 = `calc $STDV_ophA_850\*$snr`
    set var450 = `acalc ($STDV_ophA_450 ^ 2.0)*(($alpha450 ^ 2.0) + ($beta450 ^ 2.0))*(0.66666666666666) `
    set var850 = `acalc ($STDV_ophA_850 ^ 2.0)*(($alpha850 ^ 2.0) + ($beta850 ^ 2.0))*(0.66666666666666) `
else if ($map == skyam) then
    set STDV_mwc297_450_skyam = 0.0191233795195 #Jy/pixel
    set STDV_mwc297_850_skyam = 0.0016914313601 #Jy/pixel
    set sigma450 = `calc $STDV_mwc297_450_skyam\*$snr`
    set sigma850 = `calc $STDV_mwc297_850_skyam\*$snr`
    set var450 = `acalc ($STDV_mwc297_450_skyam ^ 2.0)*(($alpha450 ^ 2.0) + ($beta450 ^ 2.0))*(0.66666666666666) `
    set var850 = `acalc ($STDV_mwc297_850_skyam ^ 2.0)*(($alpha850 ^ 2.0) + ($beta850 ^ 2.0))*(0.66666666666666) `
else if ($map == skytest) then
    set STDV_mwc297_450_skytest = 0.0196827514374 #Jy/pixel
    set STDV_mwc297_850_skytest = 0.00173805640792 #Jy/pixel
    set sigma450 = `calc $STDV_mwc297_450_skytest\*$snr`
    set sigma850 = `calc $STDV_mwc297_850_skytest\*$snr`
    set var450 = `acalc ($STDV_mwc297_450_skytest ^ 2.0)*(($alpha450 ^ 2.0) + ($beta450 ^ 2.0))*(0.66666666666666) `
    set var850 = `acalc ($STDV_mwc297_850_skytest ^ 2.0)*(($alpha850 ^ 2.0) + ($beta850 ^ 2.0))*(0.66666666666666) `
else if ($map == skyspam) then
    set STDV_mwc297_450_skyspam = 0.00960455540972 #Jy/pixel
    set STDV_mwc297_850_skyspam = 0.000740904059697   #Jy/pixel
    set sigma450 = `calc $STDV_mwc297_450_skyspam\*$snr`
    set sigma850 = `calc $STDV_mwc297_850_skyspam\*$snr`
    set var450 = `acalc ($STDV_mwc297_450_skyspam ^ 2.0)*(($alpha450 ^ 2.0) + ($beta450 ^ 2.0))*(0.66666666666666) `
    set var850 = `acalc ($STDV_mwc297_850_skyspam ^ 2.0)*(($alpha850 ^ 2.0) + ($beta850 ^ 2.0))*(0.66666666666666) `
else if ($map == skyspext) then
    set STDV_mwc297_450_skyspext = 0.00962163452763  #Jy/pixel
    set STDV_mwc297_850_skyspext = 0.000738955670671   #Jy/pixel
    set sigma450 = `calc $STDV_mwc297_450_skyspext\*$snr`
    set sigma850 = `calc $STDV_mwc297_850_skyspext\*$snr`
    set var450 = `acalc ($STDV_mwc297_450_skyspext ^ 2.0)*(($alpha450 ^ 2.0) + ($beta450 ^ 2.0))*(0.66666666666666) `
    set var850 = `acalc ($STDV_mwc297_850_skyspext ^ 2.0)*(($alpha850 ^ 2.0) + ($beta850 ^ 2.0))*(0.66666666666666) `
else
    set sigma450 = 1
    set sigma850 = 1
endif

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
if(! -d $home_dir/math) then
    mkdir $home_dir/math
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
if(! -d $home_dir/math/$output) then
    mkdir $home_dir/math/$output
endif

echo "Reducing "$REDUCED_DIR"/"$file450".sdf and "$REDUCED_DIR"/"$file850".sdf:"

######### 2 -- Convolve maps with the beam size of the other #########
echo 2
#Convolution of maps
#4"/pixel for 450 micron pre-alignment; 6"/pixel for 850micron, convolve with other telescope beam
##### to implement ####### code that reads pixel size from the header

ndftrace $REDUCED_DIR/$file450.sdf QUIET
set a = `parget FPIXSCALE ndftrace` #get pixel scale from the header - currently does all three axes. 
set p450 = `echo $a | awk '{print $1}'`
ndftrace $REDUCED_DIR/$file850.sdf QUIET
set b = `parget FPIXSCALE ndftrace` #get pixel scale from the header - currently does all three axes. 
set p850 = `echo $b | awk '{print $1}'`

echo "Pixel sizes are "$p450" and "$p850

#Secondary Beam# these FWHM come from Dempsey13
echo primary beam convolution:

set fwhm450_MB = 7.9
set fwhm850_MB = 13.0

set smooth450_MB = `echo "$fwhm850_MB/$p450" | bc -l`
set smooth850_MB = `echo "$fwhm450_MB/$p850" | bc -l`

gausmooth in=$REDUCED_DIR/$file450.sdf fwhm=$smooth450_MB out=$home_dir/450/$file450/s450convolve_MB.sdf
gausmooth in=$REDUCED_DIR/$file850.sdf fwhm=$smooth850_MB out=$home_dir/850/$file850/s850convolve_MB.sdf

#Secondary Beam#
echo secondary beam convolution:

set fwhm450_S = 25.0
set fwhm850_S = 48.0

set smooth450_S = `echo "$fwhm850_S/$p450" | bc -l`
set smooth850_S = `echo "$fwhm450_S/$p850" | bc -l`

gausmooth in=$REDUCED_DIR/$file450.sdf fwhm=$smooth450_S out=$home_dir/450/$file450/s450convolve_S.sdf
gausmooth in=$REDUCED_DIR/$file850.sdf fwhm=$smooth850_S out=$home_dir/850/$file850/s850convolve_S.sdf

#Normilisation#
echo normalisation:

cmult in=$home_dir/450/$file450/s450convolve_MB.sdf scalar=$A_a out=$home_dir/450/$file450/s450norm_MB.sdf
cmult in=$home_dir/450/$file450/s450convolve_S.sdf scalar=$A_b out=$home_dir/450/$file450/s450norm_S.sdf
cmult in=$home_dir/850/$file850/s850convolve_MB.sdf scalar=$B_a out=$home_dir/850/$file850/s850norm_MB.sdf
cmult in=$home_dir/850/$file850/s850convolve_S.sdf scalar=$B_b out=$home_dir/850/$file850/s850norm_S.sdf

add in1=$home_dir/450/$file450/s450norm_MB.sdf in2=$home_dir/450/$file450/s450norm_S.sdf out=$home_dir/450/$file450/s450convolve.sdf
add in1=$home_dir/850/$file850/s850norm_MB.sdf in2=$home_dir/850/$file850/s850norm_S.sdf out=$home_dir/850/$file850/s850convolve.sdf

############### 3 -- Align the 450 map with the 850 map ############
echo 3
#Align 450 onto 850 convolved maps (in 2 or 3D where necessary)

echo Aligning maps:
ndftrace $home_dir/450/$file450/s450convolve.sdf QUIET
set ndim450 = `parget ndim ndftrace` #get number of dimensions for maps to check if itermap cut or not
ndftrace $home_dir/850/$file850/s850convolve.sdf QUIET
set ndim850 = `parget ndim ndftrace`
if($ndim450 == 3 && $ndim850 == 3) then #combination of 2/3D maps to align 450 4" pixels onto 850 6" pixels
    echo path: 3 3
    wcsframe $home_dir/850/$file850/s850convolve.sdf sky-spectrum #removes 3rd dimension
    collapse in=$home_dir/850/$file850/s850convolve.sdf axis=3 out=$home_dir/850/$file850/s850collapse.sdf QUIET
    wcsframe $home_dir/450/$file450/s450convolve.sdf sky-spectrum
    collapse in=$home_dir/450/$file450/s450convolve.sdf axis=3 out=$home_dir/450/$file450/s450collapse.sdf QUIET
    wcsalign in=$home_dir/450/$file450/s450collapse.sdf out=$home_dir/450/$file450/s450align.sdf ref=$home_dir/850/$file850/s850collapse.sdf method=nearest conserve=TRUE accept QUIET
endif
if($ndim450 == 3 && $ndim850 == 2) then
    echo path: 3 2 
    wcsframe $home_dir/450/$file450/s450convolve.sdf sky-spectrum
    collapse in=$home_dir/450/$file450/s450convolve.sdf axis=3 out=$home_dir/450/$file450/s450collapse.sdf QUIET
    wcsalign in=$home_dir/450/$file450/s450collapse.sdf out=$home_dir/450/$file450/s450align.sdf ref=$home_dir/850/$file850/s850convolve.sdf method=nearest conserve=TRUE accept QUIET
endif
if($ndim450 == 2 && $ndim850 == 3) then
    echo path: 2 3
    wcsframe $home_dir/850/$file850/s850convolve.sdf sky-spectrum
    collapse in=$home_dir/850/$file850/s850convolve.sdf axis=3 out=$home_dir/850/$file850/s850collapse.sdf QUIET
    wcsalign in=$home_dir/450/$file450/s450convolve.sdf out=$home_dir/450/$file450/s450align.sdf ref=$home_dir/850/$file850/s850collapse.sdf method=nearest conserve=TRUE accept QUIET
endif 
if($ndim450 == 2 && $ndim850 == 2) then
    echo path: 2 2 
    wcsalign in=$home_dir/450/$file450/s450convolve.sdf out=$home_dir/450/$file450/s450align.sdf ref=$home_dir/850/$file850/s850convolve.sdf method=nearest conserve=TRUE accept QUIET
endif

############### 4 -- Thresh maps  #############
echo 4
#Threshes maps to remove negative components from scaled maps, also clips maps based on SNR 

echo "Clipping SNR below "$snr":"

#Set SQRT(N) (the number of convolving pixels in the function)
# This refers to N = N_MB*alpha^2 + N_SB*beta^2 where N_comp = 1.133 (MB^2) / (pixel area^2)

set N_MB_450 = `acalc ((1.133*($fwhm450_MB^2.0))/($p450^2.0))`
set N_MB_850 = `acalc ((1.133*($fwhm850_MB^2.0))/($p850^2.0))`
set N_SB_450 = `acalc ((1.133*($fwhm450_S^2.0))/($p450^2.0))`
set N_SB_850 = `acalc ((1.133*($fwhm850_S^2.0))/($p850^2.0))`

set N_450 = `acalc ((($N_MB_450*($alpha450^2.0))+(($N_SB_450*($beta450^2.0))))^0.5)`
set N_850 = `acalc ((($N_MB_850*($alpha850^2.0))+(($N_SB_850*($beta850^2.0))))^0.5)`

set N_450_aligned = `calc $N_450\*0.6666666666667`

#Calculate the mean STDV by sigma/SQRT(N)
set Msigma450 = `calc $sigma450/$N_450_aligned`
set Msigma850 = `calc $sigma850/$N_850`

if($ndim850 == 3) then
    set s850 = 's850collapse.sdf'
endif
if($ndim850 == 2) then
    set s850 = 's850convolve.sdf'
endif 

echo Thresholding negative values: #removes negatives
stats $home_dir/450/$file450/s450align.sdf QUIET
set max450 = `parget maximum stats` #uses maximum value to force maximums to remain correct,
stats $home_dir/850/$file850/$s850 QUIET #but force script to not prompt for a maximum value
set max850 = `parget maximum stats`

#aligned - only aligned has an altered RMS due to the 'align' process.
thresh in=$home_dir/450/$file450/s450align.sdf out=$home_dir/450/$file450/s450align_th.sdf thrlo=$Msigma450 newlo=bad thrhi=$max450 newhi=$max450 QUIET

#collapsed
thresh in=$home_dir/850/$file850/$s850 out=$home_dir/850/$file850/s850collapse_th.sdf thrlo=$Msigma850 newlo=bad thrhi=$max850 newhi=$max850 QUIET

#residual
thresh in=$home_dir/850/$file850/s850convolve.sdf out=$home_dir/850/$file850/s850convolve_th.sdf thrlo=$Msigma850 newlo=bad thrhi=$max850 newhi=$max850 QUIET

############### 5 -- Make Ratio maps ###########################
echo 5

#checks for variance array for maps to allow snr cuts

### IF statement on when to make cuts (no. of dimensions and wether there is varience) ###
ndftrace $home_dir/450/$file450/s450align_th.sdf QUIET
set error450 = `parget variance ndftrace` 
if($ndim850 == 2) then
    echo $ndim850
    ndftrace $home_dir/850/$file850/s850convolve_th.sdf QUIET
    set error850 = `parget variance ndftrace`
    if($error450 == "TRUE" && $error850 == "TRUE") then #cuts snr if possible
	echo TRUE
	echo Creating ratio map:
	div in1=$home_dir/450/$file450/s450align_th.sdf in2=$home_dir/850/$file850/s850convolve_th.sdf out=$home_dir/map/$output/$output"450850.sdf"
    #if either map doesn't have variance, no cuts done
    else if($error450 == "FALSE" || $error850 == "FALSE") then 
	echo FALSE
	echo Creating ratio map:
    #uses correct 850 map based on whether collapse was used or not
        div in1=$home_dir/450/$file450/s450align_th.sdf in2=$home_dir/850/$file850/s850convolve_th.sdf out=$home_dir/map/$output/$output"450850.sdf"
 endif

    #same code, just uses collapse rather than convolve due to 850 map needing to be collapsed first
else if($ndim850 == 3) then
    echo $ndim850
    ndftrace $home_dir/850/$file850/s850collapse_th.sdf QUIET
    set error850 = `parget variance ndftrace`
    if($error450 == "TRUE" && $error850 == "TRUE") then
	echo TRUE
	echo "Clipping SNR below "$snr":"
	echo Creating ratio map:
	div in1=$home_dir/450/$file450/s450align_th.sdf in2=$home_dir/850/$file850/s850collapse_th.sdf out=$home_dir/map/$output/$output"450850.sdf"
    else if($error450 == "FALSE" || $error850 == "FALSE") then
	echo FALSE
	echo Creating ratio map:
	div in1=$home_dir/450/$file450/s450align_th.sdf in2=$home_dir/850/$file850/s850collapse_th.sdf out=$home_dir/map/$output/$output"450850.sdf"
    endif
endif

######### 6 -- creates new Varience arrays #########
echo 6
# creates new Varience array as existing one is not correct (exact reason unknown)
echo 'Creating new Variance array of ratio  maps:'

#Arithmatic 
thresh in=$home_dir/450/$file450/s450align.sdf out=$home_dir/450/$file450/s450var.sdf thrlo=0 thrhi=0 newlo=$var450 newhi=$var450 QUIET
thresh in=$home_dir/850/$file850/$s850 out=$home_dir/850/$file850/s850var.sdf thrlo=0 thrhi=0 newlo=$var850 newhi=$var850 QUIET

div in1=$home_dir/450/$file450/s450var.sdf in2=$home_dir/450/$file450/s450align.sdf out=$home_dir/math/$output/s450frcerr.sdf
div in1=$home_dir/850/$file850/s850var.sdf in2=$home_dir/850/$file850/$s850 out=$home_dir/math/$output/s850frcerr.sdf
mult in1=$home_dir/math/$output/s450frcerr.sdf in2=$home_dir/math/$output/s450frcerr.sdf out=$home_dir/math/$output/s450frcerr2.sdf
mult in1=$home_dir/math/$output/s850frcerr.sdf in2=$home_dir/math/$output/s850frcerr.sdf out=$home_dir/math/$output/s850frcerr2.sdf
add in1=$home_dir/math/$output/s450frcerr2.sdf in2=$home_dir/math/$output/s850frcerr2.sdf out=$home_dir/math/$output/ratio_err.sdf
pow in=$home_dir/math/$output/ratio_err.sdf  power=0.5 out=$home_dir/map/$output/$output'450850_fracerror.sdf'
cmult in=$home_dir/map/$output/$output'450850_fracerror.sdf' scalar=100 out=$home_dir/map/$output/$output'450850_errorpercent.sdf'
mult in1=$home_dir/map/$output/$output'450850_fracerror.sdf' in2=$home_dir/map/$output/$output'450850.sdf' out=$home_dir/map/$output/$output'450850_error.sdf'

#cuts top 3% of data if possible, to remove outliers, rough edges etc
echo Cutting high variance ratio pixels 3%:
#stats $home_dir/map/$output/$output"450850.sdf" order percentiles=97 comp=variance QUIET 
stats $home_dir/map/$output/$output"450850_error.sdf" order percentiles=97 QUIET 
set var99 = `parget perval stats`
thresh in=$home_dir/map/$output/$output"450850_error.sdf" out=$home_dir/math/$output/97_mask.sdf thrlo=$var99 newlo=1 thrhi=$var99 newhi=bad QUIET

mult in1=$home_dir/map/$output/$output"450850.sdf" in2=$home_dir/math/$output/97_mask.sdf out=$home_dir"/map/"$output"/"$output"450850clip.sdf" QUIET

######### 7 -- Create FITS files and pass to 'temperature.py' ###########
echo 7
#Produces RAW temperature map with no tidying or corrections due to variance.

if ("$type" == "b") then 
    echo "Creating beta map with temperature = "$constvalue"K:"
    maths "(log(IA*(exp(31.97/"$constvalue")-1)/(exp(16.93/"$constvalue")-1))/log(1.888857))-3" ia=$home_dir"/map/"$output"/"$output"450850clip.sdf" out=$home_dir"/map/"$output"/"$output"450850beta"$constoutput"snr"$snr".sdf"
    echo "Created "$home_dir"/map/"$output"/"$output"450850beta"$constoutput"snr"$snr".sdf"

else if ("$type" == "t") then
    echo "Creating temperature map with ß = "$constvalue":"
    #creates FITS file to output to python
    ndf2fits in=$home_dir"/map/"$output"/"$output"450850clip.sdf" out=$home_dir"/map/"$output"/"$output"450850clip.fit" QUIET
    ndftrace $home_dir/map/$output/$output"450850clip.sdf" QUIET
    set dim = `parget dims ndftrace` #gets dimensions for FOR loops

    echo 'Sending to Temperature.py'
    ###### SEND TO temperature.py ######
    #passes dimensions, constant and a flag to use variance since it has it
    python temperature.py $home_dir/map/$output $output"450850clip" $dim $constvalue yesvar 
    #makes NDF from resulting temperature FITS

    echo 'Returning from Temeprature.py'
    fits2ndf in=$home_dir"/map/"$output"/"$output"450850clip.fit" out=$home_dir"/map/"$output"/"$output"450850temp"$constoutput"snr"$snr"_RAW.sdf" QUIET 
    rm $home_dir"/map/"$output"/"$output"450850clip.fit"
    
endif

######### 8 -- Makes Variance arrays for Temp #########
echo 8
# creates new Varience array for temp. maps from Analytical solution 
echo 'Creating new Variance array of temp  maps:'

#Temp. is a variable of 4 different functions;
# Y = exp(16.93/T)
# Z = exp(31.97/T)
# V = exp(-15.04/T)
# U = exp(-31.97/T)

# I need to create number maps for the division for each constant involved here.
echo thresh
thresh in=$home_dir/map/$output/$output'450850.sdf' out=$home_dir/math/$output/16_93.sdf thrlo=0 thrhi=0 newlo=16.93 newhi=16.93 QUIET
thresh in=$home_dir/map/$output/$output'450850.sdf' out=$home_dir/math/$output/31_97.sdf thrlo=0 thrhi=0 newlo=31.97 newhi=31.97 QUIET
thresh in=$home_dir/map/$output/$output'450850.sdf' out=$home_dir/math/$output/-15_04.sdf thrlo=0 thrhi=0 newlo=-15.04 newhi=-15.04 QUIET
thresh in=$home_dir/map/$output/$output'450850.sdf' out=$home_dir/math/$output/-31_97.sdf thrlo=0 thrhi=0 newlo=-31.97 newhi=-31.97 QUIET

# now for each function
div in1=$home_dir/math/$output/16_93.sdf in2=$home_dir"/map/"$output"/"$output"450850temp"$constoutput"snr"$snr"_RAW.sdf" out=$home_dir/math/$output/y.sdf
expe in=$home_dir/math/$output/y.sdf out=$home_dir/math/$output/Y.sdf
div in1=$home_dir/math/$output/31_97.sdf in2=$home_dir"/map/"$output"/"$output"450850temp"$constoutput"snr"$snr"_RAW.sdf" out=$home_dir/math/$output/z.sdf
expe in=$home_dir/math/$output/z.sdf out=$home_dir/math/$output/Z.sdf
div in1=$home_dir/math/$output/-15_04.sdf in2=$home_dir"/map/"$output"/"$output"450850temp"$constoutput"snr"$snr"_RAW.sdf" out=$home_dir/math/$output/v.sdf
expe in=$home_dir/math/$output/v.sdf out=$home_dir/math/$output/V.sdf
div in1=$home_dir/math/$output/-31_97.sdf in2=$home_dir"/map/"$output"/"$output"450850temp"$constoutput"snr"$snr"_RAW.sdf" out=$home_dir/math/$output/u.sdf
expe in=$home_dir/math/$output/u.sdf out=$home_dir/math/$output/U.sdf

# the next section of the equation is very lengthy and will take many steps. 

#numerator (top) - Arithmatic
echo Numerator
mult in1=$home_dir"/map/"$output"/"$output"450850temp"$constoutput"snr"$snr"_RAW.sdf" in2=$home_dir"/map/"$output"/"$output"450850temp"$constoutput"snr"$snr"_RAW.sdf" out=$home_dir/math/$output/sq_T.sdf
csub in=$home_dir/math/$output/Z.sdf scalar=1 out=$home_dir/math/$output/Z_1.sdf
csub in=$home_dir/math/$output/Y.sdf scalar=1 out=$home_dir/math/$output/Y_1.sdf
div in1=$home_dir/math/$output/Z_1.sdf in2=$home_dir/math/$output/Y_1.sdf out=$home_dir/math/$output/Z_Y.sdf
add in1=$home_dir/math/$output/Z.sdf in2=$home_dir/math/$output/U.sdf out=$home_dir/math/$output/ZU.sdf
csub in=$home_dir/math/$output/ZU.sdf scalar=2 out=$home_dir/math/$output/ZU_2.sdf
mult in1=$home_dir/math/$output/ZU_2.sdf in2=$home_dir/math/$output/sq_T.sdf out=$home_dir/math/$output/sq_T_ZU2.sdf
mult in1=$home_dir/math/$output/Z_Y.sdf in2=$home_dir/math/$output/sq_T_ZU2.sdf out=$home_dir/math/$output/numerator.sdf

#dinominator (bottom) - Arithmatic
echo dinominator

cmult in=$home_dir/math/$output/Y.sdf scalar=15.04 out=$home_dir/math/$output/15_04Y.sdf
cmult in=$home_dir/math/$output/V.sdf scalar=16.93 out=$home_dir/math/$output/16_93V.sdf
add in1=$home_dir/math/$output/15_04Y.sdf in2=$home_dir/math/$output/16_93V.sdf out=$home_dir/math/$output/YV.sdf
csub in=$home_dir/math/$output/YV.sdf scalar=31.97 out=$home_dir/math/$output/dinominator.sdf

#Producing the final error and %error maps - Arithmatic
div in1=$home_dir/math/$output/numerator.sdf in2=$home_dir/math/$output/dinominator.sdf out=$home_dir/math/$output/frac.sdf
div in1=$home_dir/map/$output/$output'450850_error.sdf' in2=$home_dir/map/$output/$output'450850.sdf' out=$home_dir/math/$output/Fratio.sdf
mult in1=$home_dir/math/$output/frac.sdf in2=$home_dir/math/$output/Fratio.sdf out=$home_dir/map/$output/$output'450850temp_error.sdf' #Temp Error Array
div in1=$home_dir/map/$output/$output'450850temp_error.sdf' in2=$home_dir"/map/"$output"/"$output"450850temp"$constoutput"snr"$snr"_RAW.sdf" out=$home_dir/math/$output/frac_errT.sdf
cmult in=$home_dir/math/$output/frac_errT.sdf scalar=100 out=$home_dir/map/$output/$output'450850temp_errorpercent.sdf' #Percentage temp error array

###### Variance corrections ######

#cuts anything with %del T > input
echo "Cutting >" $percent "% error pixels in temperature map:"
thresh in=$home_dir/map/$output/$output'450850temp_errorpercent.sdf' out=$percent'_mask.sdf' thrlo=$percent newlo=1 thrhi=$percent newhi=bad QUIET
mult in1=$home_dir/map/$output/$output"450850temp"$constoutput"snr"$snr"_RAW.sdf" in2=$percent'_mask.sdf' out=$home_dir/map/$output/$output"450850temp"$constoutput"clipsnr"$snr".sdf" QUIET

#cuts anything calculated > 1000K, as unphysical and breaks variance calculations
echo "Cutting high (>999.98K) pixels:"
thresh in=$home_dir/map/$output/$output"450850temp"$constoutput"clipsnr"$snr".sdf" out=$home_dir/map/$output/$output"450850temp%"$percent".sdf" thrlo=0 newlo=0 thrhi=999.98 newhi=bad QUIET

#!!!!!!!!!#
echo "Created "$home_dir"/map/"$output"/"$output"450850temp"$constoutput"snr"$snr".sdf" 
#!!!!!!!!!#

##########tidies up temporary files to free up storage space##########
echo 'Remove files'

#rm -r $home_dir/450
#rm -r $home_dir/850
rm -r $home_dir/math

#remove excess files - Process 8 
rm $home_dir/map/$output/$output"450850temp"$constoutput"clipsnr"$snr".sdf"
rm $home_dir/map/$output/$output"450850temp"$constoutput"snr"$snr"_RAW.sdf"
rm $home_dir"/map/"$output"/"$output"450850clip.sdf"

rm $percent'_mask.sdf'

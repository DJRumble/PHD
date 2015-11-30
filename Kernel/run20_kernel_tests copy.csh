#DJR UoE
#25/10/2014 - run20.csh
#
#This code produces maps that are convolved with the MB and SB components of the JCMT beam (at that wavelength) - it then compartes those with a beam produced with Kate Pattle's convolution Kernel. Runs a few tests and plots a profile of the beams for inspection.

#! /bin/csh
alias echo "echo > /dev/null" #cleans output to screen up slightly
alias acalc '       awk "BEGIN{ print \!* }" '
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
set LOGIC = $5 #TRUE or FALSE depending on whether maps use the same pixel size - NOTE FOR KATE [when using SCUBA-2 maps set this to TRUE].

###################################################################
#Start of fixed parameters
###################################################################

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

set C_2SQRT2ln2 = 2.3548

# A (and B) = 1/(2*pi*(a4*(FM4**2)+b4*(FS4**2))) or 1/(2*pi*(a8*(FM8**2)+b8*(FS8**2)))
set A = `acalc 1/(2*$pi*($alpha450*($fwhm450_MB^2)+$beta450*($fwhm450_S^2)))`
set B = `acalc 1/(2*$pi*($alpha850*($fwhm850_MB^2)+$beta850*($fwhm850_S^2)))`

#set numbers for multiplication with convovle maps - in general they are of the form 2pi*sig^2*A*alpha
set A_a = `acalc ($fwhm450_MB^2)*$A*$alpha450*2*$pi`
set A_b = `acalc ($fwhm450_S^2)*$A*$beta450*2*$pi`
set B_a = `acalc ($fwhm850_MB^2)*$B*$alpha850*2*$pi`
set B_b = `acalc ($fwhm850_S^2)*$B*$beta850*2*$pi`

### make files ###

#uses ~/450/filename; ~/850/filename; ~/map/filename to store outputs and temporary 450/850 reduction maps
if(! -d $home_dir/$file450) then
    mkdir $home_dir/$file450
endif

### bulk code ###
ndftrace $REDUCED_DIR/$file450.sdf QUIET
set a = `parget FPIXSCALE ndftrace` #get pixel scale from the header - currently does all three axes. 
set p450 = `echo $a | awk '{print $1}'`
ndftrace $REDUCED_DIR/$file850.sdf QUIET
set b = `parget FPIXSCALE ndftrace` #get pixel scale from the header - currently does all three axes. 
set p850 = `echo $b | awk '{print $1}'`

echo  pixel sizes are: 450 =  $p450 and 850 = $p850

#Secondary Beam# these FWHM come from Dempsey13
echo primary beam convolution:

set fwhm450_MB = 7.9
set fwhm850_MB = 13.0

set smooth450_MB = `echo "$fwhm450_MB/$p450" | bc -l`
set smooth850_MB = `echo "$fwhm850_MB/$p850" | bc -l`

gausmooth in=$REDUCED_DIR/$file450.sdf fwhm=$smooth450_MB out=$home_dir/$file450/s450convolve_MB.sdf
gausmooth in=$REDUCED_DIR/$file850.sdf fwhm=$smooth850_MB out=$home_dir/$file450/s850convolve_MB.sdf
#Secondary Beam#
echo secondary beam convolution:

set fwhm450_S = 25.0
set fwhm850_S = 48.0

set smooth450_S = `echo "$fwhm450_S/$p450" | bc -l`
set smooth850_S = `echo "$fwhm850_S/$p850" | bc -l`

gausmooth in=$REDUCED_DIR/$file450.sdf fwhm=$smooth450_S out=$home_dir/$file450/s450convolve_S.sdf
gausmooth in=$REDUCED_DIR/$file850.sdf fwhm=$smooth850_S out=$home_dir/$file450/s850convolve_S.sdf

#Normilisation#
echo normalisation:

cmult in=$home_dir/$file450/s450convolve_MB.sdf scalar=$A_a out=$home_dir/$file450/s450norm_MB.sdf
cmult in=$home_dir/$file450/s450convolve_S.sdf scalar=$A_b out=$home_dir/$file450/s450norm_S.sdf
cmult in=$home_dir/$file450/s850convolve_MB.sdf scalar=$B_a out=$home_dir/$file450/s850norm_MB.sdf
cmult in=$home_dir/$file450/s850convolve_S.sdf scalar=$B_b out=$home_dir/$file450/s850norm_S.sdf

add in1=$home_dir/$file450/s450norm_MB.sdf in2=$home_dir/$file450/s450norm_S.sdf out=$home_dir/$file450/s450convolve_pix.sdf
add in1=$home_dir/$file450/s850norm_MB.sdf in2=$home_dir/$file450/s850norm_S.sdf out=$home_dir/$file450/s850convolve_pix.sdf

ndf2fits in=$home_dir/$file450/s450convolve_pix.sdf out=$home_dir/$file450/s450convolve_pix.fits

#running Kernels
echo running Kernel
idl run_aniano_convol.pro   #NOT FOR KATE [you need to manually hard code in the Kernel wrapper the directory in which 's450convolve_pix.fits' is stored, likewise the output needs to be in the same directory. 

fits2ndf in=$home_dir/$file450/K850convolve.fits out=$home_dir/$file450/K850convolve.sdf
fits2ndf in=$home_dir/$file450/kernel_result.fits out=$home_dir/$file450/kernel_result.sdf

if ($LOGIC == TRUE) then
    #align K850 onto a 850 grid
    echo running alignment
    wcsframe $home_dir/$file450/s850convolve_pix.sdf sky-spectrum
    collapse in=$home_dir/$file450/s850convolve_pix.sdf axis=3 out=$home_dir/$file450/s850collapse.sdf QUIET
    wcsalign in=$home_dir/$file450/K850convolve.sdf out=$home_dir/$file450/K850align.sdf ref=$home_dir/$file450/s850collapse.sdf method=nearest conserve=TRUE accept QUIET
    set K850 = $home_dir/$file450/K850align.sdf
else if ($LOGIC == FALSE) then
    echo no alignment

    set K850 = $home_dir/$file450/K850convolve.sdf
endif

#Calculate beam properties
stats $home_dir/$file450/s450convolve_pix.sdf QUIET
set peak450 = `parget maximum stats`
set int450 = `parget total stats`
stats $home_dir/$file450/s850convolve_pix.sdf QUIET
set peak850 = `parget maximum stats`
set int850 = `parget total stats`
#stats $home_dir/$file450/kernel_result.sdf QUIET
#set Kpeak850 = `parget maximum stats`
#set Kint850 = `parget total stats`
stats $K850 QUIET
set peakk850 = `parget maximum stats`
set intk850 = `parget total stats`

echo 450: Peak is $peak450 and intergrated flux if $int450
echo 850: Peak is $peak850 and intergrated flux if $int850
echo K850: Peak is $peakk850 and intergrated flux if $intk850
#echo Kernel: Peak is $Kpeak850 and intergrated flux if $Kint850

set effsig450 = `acalc ($int450/($peak450*2*$pi))^0.5`
set effsig850 = `acalc ($int850/($peak850*2*$pi))^0.5`
set effsigk850 = `acalc ($intk850/($peakk850*2*$pi))^0.5`
#set Keffsig850 = `acalc ($Kint850/($Kpeak850*2*$pi))^0.5`

set efffwhm450 = `acalc (($effsig450*$C_2SQRT2ln2)*$p450)`
set efffwhm850 = `acalc (($effsig850*$C_2SQRT2ln2)*$p850)`
set efffwhmk850 = `acalc (($effsigk850*$C_2SQRT2ln2)*$p850)`
#set Kefffwhm850 = `acalc (($Keffsig850*$C_2SQRT2ln2)*$p850)`

#echo Does align [ $efffwhmk850 ] equal non-align [ $efffwhmkc850 ]
echo $K850

echo effective beam at 450 is $efffwhm450
echo effective beam at 850 is $efffwhm850
echo effective beam at K850 is $efffwhmk850
#echo effective beam of kernel is $Kefffwhm850

rm $home_dir/$file450/s450convolve_pix.fits
rm $home_dir/$file450/K850convolve.fits

ndfcopy $K850'(-15:15,0:0)' $home_dir/$file450/k850line.sdf
ndfcopy $home_dir/$file450/s850convolve_pix.sdf'(-15:15,0:0)' $home_dir/$file450/s850line.sdf

#if ($efffwhm850 -gt $efffwhmk850) then
#linplot $home_dir/$file450/k850line.sdf device='xwi' 
#linplot $home_dir/$file450/s850line.sdf device='xwi' clear='false'
#else if ($efffwhm850 -lt $efffwhmk850) then
linplot $home_dir/$file450/s850line.sdf device='xwi' 
linplot $home_dir/$file450/k850line.sdf device='xwi' clear='false'
#endif

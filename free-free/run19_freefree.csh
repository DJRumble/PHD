#DJR UoE
#25/03/2014 - run19.csh

#This is complete code that takes VLA radio data and SCUBA-2 data and processes them, producing maps with contamination subtracted. This requires: unit conversion in VLA data, cross convolution of VLA with SCUBA-2 (450/850), Alignment of VLA maps to SCUBA-2 grids. Scaling of VLA data into smm via Free-Free emission, convolution of each wavelength resemble the JCMT beam. Fianlly subtraction of contaminant maps from IR1 maps to produce new maps.

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
set file850 = $4  #name of raw directory/850 map, without extension
set fileVLA = $5  #name of raw directory/VLA map (2cm,3.6cm,6cm), without extension

###################################################################
#Start of fixed parameters
###################################################################

set pi = 3.14159265359

#Free-Free power index
set alpha_freefree = 1.03

#Set FWHM for SCUBA-2 maps - MB (Main Beam) S (Secondary Beam)
set fwhm450 = 9.6
set fwhm850 = 14.1

set fwhm450_MB = 7.9
set fwhm850_MB = 13.0 
set fwhm450_S = 25
set fwhm850_S = 48

#Set FWHM for VLA maps (single component)
set fwhm2cm = 4.1
set fwhm3_6cm = 7.2
set fwhm6cm = 18.43    #[5.1,21.2]     #14.4
set pa_6cm = 30.7

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
set B_a = `acalc ($fwhm850_MB^2)*$B*$alpha850*2*$pi`
set B_b = `acalc ($fwhm850_S^2)*$B*$beta850*2*$pi`

### make files ###

#uses ~/450/filename; ~/850/filename; ~/map/filename to store outputs and temporary 450/850 reduction maps
if (! -d $home_dir/450) then
    mkdir $home_dir/450
endif
if(! -d $home_dir/850) then
    mkdir $home_dir/850
endif
if(! -d $home_dir/VLA) then
    mkdir $home_dir/VLA
endif
if(! -d $home_dir/450/$file450) then
    mkdir $home_dir/450/$file450
endif
if(! -d $home_dir/850/$file850) then
    mkdir $home_dir/850/$file850
endif
if(! -d $home_dir/VLA/$fileVLA) then
    mkdir $home_dir/VLA/$fileVLA
endif

### BULK CODE  ###
#############################################################################################

### 1) Change map units based on FWHM of VLA beam

echo 1: Units
#Converts units of VLA data from Jy per beam to Jy per Pixel. FWHM currently uses a typical size of VLA beam at various wavelenghts as a circle when infact it is described by a elipsodial shape - UPDATE

ndftrace $REDUCED_DIR/sandell_radio/$fileVLA.sdf QUIET
set a = `parget FPIXSCALE ndftrace` #get pixel scale from the header - currently does all three axes. 
set pvla = `echo $a | awk '{print $1}'`

echo Pixel size of VLA: $pvla

set unit = `acalc ($pvla^2)/($fwhm6cm^2)`

cmult in=$REDUCED_DIR/sandell_radio/$fileVLA.sdf scalar=$unit out=$home_dir/VLA/$fileVLA/$fileVLA'_pix.sdf'

### 2) Cross convolve VLA data with SCUBA-2 data to achieve like resolution

echo 2: Resolution with JCMT beam
#Convolve the VLA beam with a single component JCMT beam for 450 and 850. Take from Dempsey13.

ndftrace $REDUCED_DIR/IR1/$file450'.sdf' QUIET
set a = `parget FPIXSCALE ndftrace` #get pixel scale from the header - currently does all three axes. 
set p450 = `echo $a | awk '{print $1}'`
ndftrace $REDUCED_DIR/IR1/$file850'.sdf' QUIET
set b = `parget FPIXSCALE ndftrace` #get pixel scale from the header - currently does all three axes. 
set p850 = `echo $b | awk '{print $1}'`

echo  pixel sizes are: 450 =  $p450 and 850 = $p850

#set smooth450 = `echo "$fwhm450/$p450" | bc -l`
#set smooth850 = `echo "$fwhm850/$p850" | bc -l`

#gausmooth in=$home_dir/VLA/$fileVLA/$fileVLA'_pix.sdf' fwhm=$smooth450 out=$home_dir/VLA/$fileVLA/$fileVLA'_smooth450.sdf'
#gausmooth in=$home_dir/VLA/$fileVLA/$fileVLA'_pix.sdf' fwhm=$smooth850 out=$home_dir/VLA/$fileVLA/$fileVLA'_smooth850.sdf'

echo primary beam convolution:

set smooth450_MB = `echo "$fwhm450_MB/$p450" | bc -l`
set smooth850_MB = `echo "$fwhm850_MB/$p850" | bc -l`

gausmooth in=$home_dir/VLA/$fileVLA/$fileVLA'_pix.sdf' fwhm=$smooth450_MB out=$home_dir/450/$file450/s450convolve_MB.sdf
gausmooth in=$home_dir/VLA/$fileVLA/$fileVLA'_pix.sdf' fwhm=$smooth850_MB out=$home_dir/850/$file850/s850convolve_MB.sdf

#Secondary Beam#
echo secondary beam convolution:

set smooth450_S = `echo "$fwhm450_S/$p450" | bc -l`
set smooth850_S = `echo "$fwhm850_S/$p850" | bc -l`

gausmooth in=$home_dir/VLA/$fileVLA/$fileVLA'_pix.sdf' fwhm=$smooth450_S out=$home_dir/450/$file450/s450convolve_S.sdf
gausmooth in=$home_dir/VLA/$fileVLA/$fileVLA'_pix.sdf' fwhm=$smooth850_S out=$home_dir/850/$file850/s850convolve_S.sdf

#Normilisation#
echo normalisation:

cmult in=$home_dir/450/$file450/s450convolve_MB.sdf scalar=$A_a out=$home_dir/450/$file450/s450norm_MB.sdf
cmult in=$home_dir/450/$file450/s450convolve_S.sdf scalar=$A_b out=$home_dir/450/$file450/s450norm_S.sdf
cmult in=$home_dir/850/$file850/s850convolve_MB.sdf scalar=$B_a out=$home_dir/850/$file850/s850norm_MB.sdf
cmult in=$home_dir/850/$file850/s850convolve_S.sdf scalar=$B_b out=$home_dir/850/$file850/s850norm_S.sdf

add in1=$home_dir/450/$file450/s450norm_MB.sdf in2=$home_dir/450/$file450/s450norm_S.sdf out=$home_dir/VLA/$fileVLA/$fileVLA'_smooth450.sdf'
add in1=$home_dir/850/$file850/s850norm_MB.sdf in2=$home_dir/850/$file850/s850norm_S.sdf out=$home_dir/VLA/$fileVLA/$fileVLA'_smooth850.sdf'


### 3) Alignment of VLA on to SCUBA-2 grids

echo 3: Alignment

ndftrace $home_dir/VLA/$fileVLA/$fileVLA'_smooth450.sdf' QUIET
set ndim450 = `parget ndim ndftrace` #get number of dimensions for maps to check if itermap cut or not
ndftrace $home_dir/VLA/$fileVLA/$fileVLA'_smooth850.sdf' QUIET
set ndim850 = `parget ndim ndftrace`

#COLLAPSE to remove 4th and 3rd dimensions
echo Collapsing images to 2D

collapse in=$home_dir/VLA/$fileVLA/$fileVLA'_smooth450.sdf' axis=4 out=$home_dir/VLA/$fileVLA/s450collapse.sdf QUIET
collapse in=$home_dir/VLA/$fileVLA/s450collapse.sdf axis=3 out=$home_dir/VLA/$fileVLA/s450collapse2.sdf QUIET

collapse in=$home_dir/VLA/$fileVLA/$fileVLA'_smooth850.sdf' axis=4 out=$home_dir/VLA/$fileVLA/s850collapse.sdf QUIET
collapse in=$home_dir/VLA/$fileVLA/s850collapse.sdf axis=3 out=$home_dir/VLA/$fileVLA/s850collapse2.sdf QUIET

#collapse ref. maps to 2D
collapse in=$REDUCED_DIR/IR1/$file450'.sdf' axis=3 out=$home_dir/VLA/$fileVLA/$file450'_2D.sdf' QUIET
collapse in=$REDUCED_DIR/IR1/$file850'.sdf' axis=3 out=$home_dir/VLA/$fileVLA/$file850'_2D.sdf' QUIET

#WCSALIGN to map VLA maps onto SCUAB-2 pixel coordinates
echo Image alignment 
#This SHOULD be TRUE. Currently not working, check with Jenny. 
wcsalign in=$home_dir/VLA/$fileVLA/s450collapse2.sdf out=$home_dir/VLA/$fileVLA/$fileVLA'_align450.sdf' ref=$home_dir/VLA/$fileVLA/$file450'_2D.sdf' method=nearest conserve=TRUE accept QUIET
wcsalign in=$home_dir/VLA/$fileVLA/s850collapse2.sdf out=$home_dir/VLA/$fileVLA/$fileVLA'_align850.sdf' ref=$home_dir/VLA/$fileVLA/$file850'_2D.sdf' method=nearest conserve=TRUE accept QUIET

### 4) Scaling of VLA to SMM along freefree power law

echo 4: Free-Free scaling
#Follows the law Sv proporto v^alpha - for each pixel in the map.

#currently using single scalar version.
set A_450 = 163.8
set A_850 = 82.5

cmult in=$home_dir/VLA/$fileVLA/$fileVLA'_align450.sdf' scalar=$A_450 out=$home_dir/450/$file450/mwc297_freefree_comp450convolve.sdf
cmult in=$home_dir/VLA/$fileVLA/$fileVLA'_align850.sdf' scalar=$A_850 out=$home_dir/850/$file850/mwc297_freefree_comp850convolve.sdf

### 5) convolve flux with the JCMT beam

#echo 5: convolve with the JCMT beam


echo Final Free component maps are stored 'in:'
echo $home_dir/450/$file450/mwc297_freefree_comp450convolve.sdf
echo $home_dir/850/$file850/mwc297_freefree_comp850convolve.sdf

### 6) subtract from IR1 data 

echo 5: Contamination subtraction

set val = 0

nomagic in=$home_dir/450/$file450/mwc297_freefree_comp450convolve.sdf out=$home_dir/450/$file450/mwc297_freefree_comp450convolve0.sdf repval=$val QUIET
nomagic in=$home_dir/850/$file850/mwc297_freefree_comp850convolve.sdf out=$home_dir/850/$file850/mwc297_freefree_comp850convolve0.sdf repval=$val QUIET

sub in1=$REDUCED_DIR/IR1/$file450'.sdf'  in2=$home_dir/450/$file450/mwc297_freefree_comp450convolve0.sdf out=SerpensMWC297_20140424_s450_IR1_6cmfreefree+jet103.sdf

sub in1=$REDUCED_DIR/IR1/$file850'.sdf'  in2=$home_dir/850/$file850/mwc297_freefree_comp850convolve0.sdf out=SerpensMWC297_20140424_s850_IR1_6cmfreefree+jet103.sdf

echo COMPLETED

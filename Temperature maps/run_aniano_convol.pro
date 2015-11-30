
im=mrdfits('~/Ophiuchus/new_maps/jypix/ophmain_850_noco_cal537_2014_extract_Jypix.fits',0,imhead)
kernel=mrdfits('~/san/Kernels/KERNEL_SCUBA2850-SPIRE250_6as.fits',0,khead)

.comp convolve_image

do_the_convolution,im,imhead,kernel,khead,convol,chead,rker,rkhead,1

mwrfits,convol,'~/Ophiuchus/new_maps/conv250/not_l1688/ophmain_850_noco_cal537_2014_extract_Jypix_conv250.fits',chead,/create


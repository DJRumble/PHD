;Wrapper for use with ratiomap450K850.py

; Read input from the terminal:
s450 = ''
read, s450, prompt='Please enter the file name of the 450micron map
;to be convolved with the Kernel: '

;Set input 450 map
im=mrdfits(s450,0,imhead)
;Set kernel
kernel=mrdfits('KERNEL_MODEL450-MODEL850_4as.fits',0,khead)

print,'Running kernel convolution through convolve_image.pro'
.comp convolve_image

do_the_convolution,im,imhead,kernel,khead,convol,chead,rker,rkhead,1

mwrfits,convol,'K450convolve.fits',chead,/create

exit

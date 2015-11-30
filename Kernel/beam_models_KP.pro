

pro beam_models

primary=psf_gaussian(npixel=501,fwhm=13.0,/double,/normalize)
secondary=psf_gaussian(npixel=501,fwhm=48.0,/double,/normalize)

beam=(0.75*primary)+(0.25*secondary)

mkhdr,head,beam,/image

sxaddpar,head,'CDELT1',1.0/(60.0^2.0)
sxaddpar,head,'CDELT2',1.0/(60.0^2.0)

print,head

writefits,'~/san/Kernels/model_850_beam.fits',beam,head

primary=psf_gaussian(npixel=501,fwhm=7.9,/double,/normalize)
secondary=psf_gaussian(npixel=501,fwhm=25.0,/double,/normalize)

beam=(0.6*primary)+(0.4*secondary)

mkhdr,head,beam,/image
sxaddpar,head,'NAXIS',2
sxaddpar,head,'CDELT1',1.0/(60.0^2.0)
sxaddpar,head,'CDELT2',1.0/(60.0^2.0)

writefits,'~/san/Kernels/model_450_beam.fits',beam,head

END

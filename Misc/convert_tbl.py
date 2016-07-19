#!/usr/local64/python2.6/bin/python
#Convert c2d ascii tables to FITS for faster reading

import atpy

#tbl1=atpy.Table('CATALOGS/catalog-SCO_1-v1.tbl',type='c2d')
#tbl1.write('Catalogues/catalog-SCO_1-v2.fits')
#tbl2=atpy.Table('CATALOGS/catalog-SCO_2-v1.tbl',type='c2d')
#tbl2.write('Catalogues/catalog-SCO_2-v2.fits')
#tbl3=atpy.Table('CATALOGS/catalog-SCO_3-v2.tbl',type='c2d')
#tbl3.write('Catalogues/catalog-SCO_3-v2.fits')
#tbl4=atpy.Table('CATALOGS/catalog-SCO_4-v1.tbl',type='c2d')
#tbl4.write('Catalogues/catalog-SCO_4-v2.fits')
#tbl5=atpy.Table('CATALOGS/catalog-SCO_5-v1.tbl',type='c2d')
#tbl5.write('Catalogues/catalog-SCO_5-v2.fits')
#tbl6=atpy.Table('CATALOGS/catalog-SCO_6-v1.tbl',type='c2d')
#tbl6.write('Catalogues/catalog-SCO_6-v2.fits')
#tbl7=atpy.Table('Catalogues/master4_edited2.tbl',type='c2d')
#tbl7.write('Catalogues/catalog-SCO_ALL.fits')

tbl=atpy.Table('C2DCORES/CB68/catalog-CB68-v2.tbl',type='c2d')
tbl.write('Catalogues/catalog-CB68-v2.fits')

tbl2=atpy.Table('C2DCORES/L234E/catalog-L234E-v2.tbl',type='c2d')
tbl2.write('Catalogues/catalog-L234E-v2.fits')



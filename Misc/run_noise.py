import noise


s450mwc297 = 'serpens/MWC297/IR2/SerpensMWC297_20141219_s450_IR2extmask_s2_cal_JypixJH.sdf'
s850mwc297 = 'serpens/MWC297/IR2/SerpensMWC297_20141219_s850_IR2extmask_s2_cal_JypixJH.sdf'
s450main = 'serpens/serpensmain/IR2/SerpensMain_20141223_s450_IR2extmask_s2_cal_JypixJH.sdf'
s850main = 'serpens/serpensmain/IR2/SerpensMain_20141223_s850_IR2extmask_s2_cal_JypixJH.sdf'
s850mainnoco = 'serpens/serpensmain/IR2/SerpensMain_20150326_850_IR2_noco_JypixHK.sdf'
s450E = 'serpens/serpensE/IR2/SerpensE_20141219_s450_IR2extmask_s2_cal_JypixJH.sdf'
s850E = 'serpens/serpensE/IR2/SerpensE_20141219_s850_IR2extmask_s2_cal_JypixJH.sdf'
s450N = 'serpens/serpensN/IR2/SerpensN_20141219_s450_IR2extmask_s2_cal_JypixJH.sdf'
s850N = 'serpens/serpensN/IR2/SerpensN_20141219_s850_IR2extmask_s2_cal_JypixJH.sdf'


noise.noise_by_data(s450main,'FASLE')
noise.noise_by_data(s850main,'FASLE')
noise.noise_by_data(s450mwc297,'FASLE')
noise.noise_by_data(s850mwc297,'FASLE')
noise.noise_by_data(s450E,'FASLE')
noise.noise_by_data(s850E,'FASLE')
noise.noise_by_data(s450N,'FASLE')
noise.noise_by_data(s850N,'FASLE')
noise.noise_by_data(s850mainnoco,'FASLE')

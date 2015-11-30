import math

sig_m4 = 7.9
sig_m8 = 13.0
sig_s4 = 25.0
sig_s8 = 48.0

a4 = 0.94
a8 = 0.98
b4 = 0.06
b8 = 0.02

sig_e = (((a4*a8*sig_m4*sig_m8)+(a8*b4*sig_s4*sig_m8)+(b8*a4*sig_m4+sig_s8)+(b4*b8*sig_s4*sig_s8))/((a4*a8*(sig_m8**2.0))+(b4*b8*(sig_s8**2.0))+(a4*b8*(sig_s8**2.0))+(a8*b4*(sig_m8**2.0))))**0.5

print sig_e

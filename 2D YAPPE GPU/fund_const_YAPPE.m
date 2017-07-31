%this function loads in a bunch of fundamental constants in SI units into
%s.fundSI and cgs versions in s.cgs.

function[s] = fund_const_YAPPE(s)

    s.SI.e = -1.6e-19; %electron charge in coulomb
    s.SI.m_e = 9.11e-31; %electron mass in kg
    s.SI.eps_0 = 8.854e-12; %permitivity of free space in Farads/meter
    s.SI.hbar = 1.05e-34; %planck's constant/2*pi in Joules*seconds
    s.SI.c = 3e8; %speed of light in m/s
    
    s.cgs.c = 3e10; %speed of light in cm/s
    
end
%this script calculates multiphoton ionization coefficients. I used
%equation 8 from the paper "Measurements of multiphoton ionization
%coefficients with ultrashort ultraviolet laser pulses" which gives the
%coefficient in terms optical wavelength and ionization potential of the
%atom/molecule + some fundamental constants. 

function[s] = MPI_YAPPE(s)

    delE_J = s.mat.Eg; %ionization potential in joules
    omg_cen = 2*pi*s.cgs.c/s.input.lambda_vac; %central angular frequency in rad/s

    s.mat.pow = ceil(delE_J/(s.SI.hbar*omg_cen));
    s.mat.sigk = omg_cen*(delE_J/(s.SI.hbar*omg_cen))^(3/2)*(1e4*s.SI.e^2/(4*s.SI.c*s.SI.eps_0*s.SI.m_e*omg_cen^2*delE_J))^s.mat.pow; %cross section in CGS\

end
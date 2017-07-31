%This precomputes useful functions and coefficients which will be used in
%the main computation loop, but are not a function of the field itself. Kz
%is the prefactor for the linear term in UPPE (diffraction and dispersion).
%Q is the prefactor for the nonlinear term in UPPE (self-focusing and
%plasma for this code). For more info on Kz and Q, see p. 15 of Kolesik's
%paper "Practioner's Guide to laser pulse propagation models and
%simulations". We also compute Kz_move, the moving frame version of Kz,
%the discrete Hankel transform H, and a bunch of coefficients used in
%calculating ionization and nonlinear susceptibility. 

function[s] = precompute_YAPPE(s)

    %calculate Kz and Q functions
    s.f.Kz = sqrt(bsxfun(@minus, s.g.k.^2, bsxfun(@plus, s.g.kx.^2,s.g.ky.^2)));
    s.f.Q = gpuArray(bsxfun(@times, s.g.omg.^2/s.cgs.c^2, 1./s.f.Kz));
    
    %calculate linear propagator (handles diffraction and dispersion)
    s.f.Kz_move = gpuArray(bsxfun(@minus, s.f.Kz, s.g.omg/s.g.vg));
    s.f.lin_prop = gather(exp(1i*s.input.outperiod*s.f.Kz_move));
    
    %coefficients used in calculating ionization
    m2tocm2 = 1e4;
    s.ion.sigma = s.SI.e^2*s.mat.tau_c*s.g.n0/(s.SI.m_e*s.SI.eps_0*s.SI.c)/(1 + (s.g.omg_cen*s.mat.tau_c)^2)*m2tocm2; % avalanche cross section in cm^2    
    s.ion.a1 = s.ion.sigma/(s.g.n0^2*s.mat.Eg);
    s.ion.a2 = s.mat.sigk*s.mat.N0;
    s.ion.a3 = s.mat.recomb;

    %coefficients used in calculating nonlinear susceptibility
    s.NL.betak = s.mat.N0*s.mat.sigk*s.mat.pow*s.SI.hbar*s.g.omg_cen;
    s.NL.b1 = 2*s.g.n0*s.mat.n2;
    s.NL.b2 = 1i*1e6*s.SI.e^2/(s.SI.eps_0*s.SI.m_e*s.g.omg_cen)/(1/s.mat.tau_c - 1i*s.g.omg_cen); %converting rho to m^-3 from cm^-3
    s.NL.b3 = 1i*s.g.n0^2*s.NL.betak/s.g.kvac(1); %not sure if the "k" in this expression should be medium or vaccuum central wavelength

end
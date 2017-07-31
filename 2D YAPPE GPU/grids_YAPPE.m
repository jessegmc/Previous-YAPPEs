%create grids and grid parameters and load into substruct s.g

function [s] = grids_YAPPE(s)
    
    %lab frame axis properties
    s.g.z_range = s.input.z_extent; %propagation range in cm
    s.g.zout = 0:s.input.outperiod:s.g.z_range; %output propagation locations in cm

    %local time axis grid and grid properties
    s.g.xi_range = s.input.xi_extent; %temporal range in seconds
    s.g.xi_pts = s.input.xi_pts; %number of grid pts along local time axis
    s.g.xi = linspace(0, s.g.xi_range, s.input.xi_pts); %local time axis in seconds
    s.g.dxi = s.g.xi(2) - s.g.xi(1); %local time axis grid spacing in seconds
    
    %spectral axis grid and grid properties
    s.g.omg_cen = 2*pi*s.cgs.c/s.input.lambda_vac; %central angular frequency in rad/s
    s.g.domg = 2*pi/(s.g.xi_range); %angular frequency spacing in rad/s
    s.g.omg = s.g.domg*((1:s.g.xi_pts) - round((s.g.xi_pts+1)/2)); %angular frequency axis in rad/s
    s.g.omg = s.g.omg + s.g.omg_cen; %center axis about central frequency
    s.g.omg = ifftshift(s.g.omg,2); %shift the axis
    
    %wavenumbers
    s.g.kvac = s.g.omg/s.cgs.c; %vaccuum wavenumbers
    s = dispersion_YAPPE(s); %get permittivity 
    s.g.n = sqrt(s.g.perm); %linear index 
    s.g.n0 = s.g.n(1); %central index
    s.g.k = s.g.omg.*s.g.n/s.cgs.c; %medium wavenumbers
    s.g.kcen = s.g.k(1); %central medium wavenumber
    s.g.vg = ( (s.g.k(2) - s.g.k(end))/(2*s.g.domg) )^(-1); %group velocity
    
    %transverse spatial grid and grid properties
    s.g.r_pts = s.input.r_pts; %radial window length in cm
    s.g.r_extent = s.input.r_extent;
    s = hankel_sample_YAPPE(s); %this generates the bessel zero-spaced radial grid

end
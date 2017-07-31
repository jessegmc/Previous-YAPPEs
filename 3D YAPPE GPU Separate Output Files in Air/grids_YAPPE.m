%3d

%create grids and grid parameters and load into substruct s.g

function [s] = grids_YAPPE(s)
    
    %lab frame axis properties
    s.g.z_range = s.input.z_extent; %propagation range in cm
    s.g.zout = 0:s.input.outperiod:s.g.z_range; %output propagation locations in cm

    %local time axis grid and grid properties
    s.g.xi_range = s.input.xi_extent; %temporal range in seconds
    s.g.xi_pts = s.input.xi_pts; %number of grid pts along local time axis
    
% %     s.g.xi = linspace(0, s.g.xi_range, s.input.xi_pts); %local time axis in seconds
% %     Changed xi to be along the 3rd dimension
    s.g.xi = linspace(0, s.g.xi_range, s.input.xi_pts); %local time axis in seconds
    s.g.xi = permute(s.g.xi,[1,3,2]);
    s.g.dxi = s.g.xi(2) - s.g.xi(1); %local time axis grid spacing in seconds
    
    %spectral axis grid and grid properties
    s.g.omg_cen = 2*pi*s.cgs.c/s.input.lambda_vac; %central angular frequency in rad/s
    s.g.domg = 2*pi/(s.g.xi_range); %angular frequency spacing in rad/s
    s.g.omg = s.g.domg*((1:s.g.xi_pts) - round((s.g.xi_pts+1)/2)); %angular frequency axis in rad/s
    s.g.omg = s.g.omg + s.g.omg_cen; %center axis about central frequency
% %     s.g.omg = ifftshift(s.g.omg,2); %shift the axis
    s.g.omg = permute(s.g.omg,[1,3,2]);
    s.g.omg = ifftshift(s.g.omg,3); %shift the axis
    
    %wavenumbers
    s.g.kvac = s.g.omg/s.cgs.c; %vaccuum wavenumbers
    s = dispersion_YAPPE(s); %get permittivity 
    s.g.n = sqrt(s.g.perm); %linear index 
    s.g.n0 = s.g.n(1); %central index
    s.g.k = s.g.omg.*s.g.n/s.cgs.c; %medium wavenumbers
    s.g.kcen = s.g.k(1); %central medium wavenumber
    s.g.vg = ( (s.g.k(2) - s.g.k(end))/(2*s.g.domg) )^(-1); %group velocity
    
% %     %transverse spatial grid and grid properties
% %     s.g.r_pts = s.input.r_pts; %radial window length in cm
% %     s.g.r_extent = s.input.r_extent;
% %     s = hankel_sample_YAPPE(s); %this generates the bessel zero-spaced radial grid

    %x-axis spatial grid and grid properties
    s.g.x_extent = s.input.x_extent; %x-axis window length from center of pulse to edge
    s.g.x_pts = s.input.x_pts; %x-axis window resolution
    s.g.x = linspace(-s.g.x_extent/2,s.g.x_extent/2,s.g.x_pts); %x-axis grid in cm
%     s.g.x = s.g.x'; %make x a column vector
    s.g.dx = s.g.x(2) - s.g.x(1); %x-axis grid spacing

% %    define s.g.kx
    s.g.dkx = 2*pi/s.g.x_extent;
    s.g.kx = s.g.dkx*((1:s.g.x_pts) - round((s.g.x_pts+1)/2));
    s.g.kx = ifftshift(s.g.kx,2);
    
    %y-axis spatial grid and grid properties
    s.g.y_extent = s.input.y_extent; %y-axis window length from center of pulse to edge
    s.g.y_pts = s.input.y_pts; %y-axis window resolution
    s.g.y = linspace(-s.g.y_extent/2,s.g.y_extent/2,s.g.y_pts); %y-axis grid in cm
    s.g.y = s.g.y';
    s.g.dy = s.g.y(2) - s.g.y(1); %y-axis grid spacing
    
% %    define s.g.ky
    s.g.dky = 2*pi/s.g.y_extent;
    s.g.ky = s.g.dky*((1:s.g.y_pts) - round((s.g.y_pts+1)/2));
    s.g.ky = s.g.ky';
    s.g.ky = ifftshift(s.g.ky,1);
    
end
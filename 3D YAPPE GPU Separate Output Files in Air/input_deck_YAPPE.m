%3d

%This function creates the input struct for YAPPE.

function [s] = input_deck_YAPPE()
    
    %grid parameters
    s.input.xi_extent = 600*1e-15; %axial window length in seconds
    s.input.xi_pts = 280; %number of pts along axial direction
% %     s.input.r_extent = .2; %radial window length in cm
% %     s.input.r_pts = 600; %number of pts along radial direction    
    s.input.x_extent = .06; %x-axis window length in cm from edge to edge
    s.input.x_pts = 100; %number of pts along x direction
    s.input.y_extent = .06; %y-axis window length in cm from edge to edge
    s.input.y_pts = 100; %number of pts along y direction
    s.input.z_extent = 5; %extent of propagation in cm

    %load E-field options
    s.input.lambda_vac = 800e-7; %vacuum central wavelength in cm
    s.input.infield.type = 'double gauss'; %choose between 'gauss' and 'custom' and 'double gauss'   
    s.input.infield.path = '~/matlab_scripts/hankel_FILA/E_in_example.mat'; %location of custom field    
% %     s.input.infield.waist = 40e-4; %gaussian waist in cm (used in gauss)
    s.input.infield.waistx = 40e-4; %gaussian x-waist in cm (used in gauss)
    s.input.infield.waisty = 40e-4; %gaussian y-waist in cm (used in gauss)
    s.input.infield.tfwhm = 40e-15; %intensity fwhm in s (used in gauss)
    s.input.infield.energ = 1e-6; %beam energy in J (used in gauss)
    s.input.infield.f = inf; %focusing length in cm (used in gauss)
    s.input.infield.beamOffset = 0.005; %double beam x offset from center (in cm)

    %choose propagation medium
    s.input.medium = 'water';
    
    %specify output path and output period
    s.input.outpath = 'R:\Jesse GPU\MATLAB\3D save, multiple beams YAPPE\Outputs\Run3\'; %output path
    s.input.outperiod = 0.05; %output period in cm
    
    %toggle dispersion, plasma and n2 propagation modules
    s.input.dispersion = 1;
    s.input.plasma = 1;
    s.input.n2 = 1;
    
    %solution tolerances for ODE calls
    s.input.RelTol = 1e-3;
    s.input.AbsTol = 1e-6;
    
    %absorbing boundaries
    s.input.freqbd = 1; %toggle the absorbing frequency boundary on and off
    s.input.freqbd_length = .005; %set absorption length in cm
    s.input.freqbd_width = .1; %fractional value of boundary width (boundary width/total bandwidth)
    

end
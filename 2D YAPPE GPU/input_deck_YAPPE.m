%This function creates the input struct for YAPPE.

function [s] = input_deck_YAPPE()
    
    %grid parameters
    s.input.xi_extent = 600*1e-15; %axial window length in seconds
    s.input.xi_pts = 100; %number of pts along axial direction
    s.input.r_extent = .1; %radial window length in cm
    s.input.r_pts = 200; %number of pts along radial direction    
    s.input.z_extent = 0.2; %extent of propagation in cm

    %load E-field options
    s.input.lambda_vac = 800e-7; %vacuum central wavelength in cm
    s.input.infield.type = 'gauss'; %choose between 'gauss' and 'custom'    
    s.input.infield.path = '~/matlab_scripts/hankel_FILA/E_in_example.mat'; %location of custom field    
    s.input.infield.waist = 40e-4; %gaussian waist in cm (used in gauss)
    s.input.infield.tfwhm = 40e-15; %intensity fwhm in s (used in gauss)
    s.input.infield.energ = .5e-6; %beam energy in J (used in gauss)
    s.input.infield.f = inf; %focusing length in cm (used in gauss)

    %choose propagation medium
    s.input.medium = 'water';
    
    %specify output path and output period
    s.input.outpath = 'C:\Users\Pinky\Documents\MATLAB\Jesse GPU\2D GPU YAPPE\Outputs\Run1'; %output path
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
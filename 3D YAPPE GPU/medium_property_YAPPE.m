%3d

%loads material parameters into substruct s.mat.

function [s] = medium_property_YAPPE(s)

    if s.input.medium == 'water'
       
        s.mat.tau_c = 1e-14; %electron collision time in seconds
        s.mat.N0 = 3.33e22; %number density of water in cm^-3
        s.mat.Eg = 1.12e-18; %ionization energy in J (7eV)      
        s.mat.recomb = 2e-15*1e6; %recombination rate in cm^3/s (was "s.a3")        
        s.mat.n2 = 1.9e-16; %nonlinear index - valid at 800nm in cm^2/W
        
    elseif s.input.medium == 'air'
       
        s.mat.tau_c = 3.5e-13; %electron collision time in seconds
        
        s.mat.N0 = 2.504e19; %number density of air in cm^-3
        s.mat.NO_N2 = 2e19; %number density of N2 in the air
        s.mat.NO_O2 = 5e18; %number density of O2 in the air
        
        s.mat.Eg_N2 = 2.496e-18; %ionization energy in J of N2 (15.58 eV)      
        s.mat.Eg_O2 = 1.934e-18; %ionization energy in J of O2 (12.07 eV)   
        
        s.mat.recomb = 5e-13*1e6; %recombination rate in cm^3/s (was "s.a3")
        
        s.mat.n2_N2 = 7.9e-20; %nonlinear index of N2 - valid at 800nm in cm^2/W
        s.mat.n2_O2 = 8.1e-20; %nonlinear index of O2 - valid at 800nm in cm^2/W
    end
    s = MPI_YAPPE(s); %calculate MPI rate and # of photons needs to ionize


end
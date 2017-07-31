%3d

%loads material parameters into substruct s.mat.

function [s] = medium_property_YAPPE(s)

    if s.input.medium == 'water'
       
        s.mat.tau_c = 1e-14; %electron collision time in seconds
        s.mat.N0 = 3.33e22; %number density of water in cm^-3
        s.mat.Eg = 1.12e-18; %ionization energy in J (7eV)      
        s.mat.recomb = 2e-15*1e6; %recombination rate in cm^3/s (was "s.a3")
        s = MPI_YAPPE(s); %calculate MPI rate and # of photons needs to ionize
        
%         s.mat.n2 = 2.7e-16; %nonlinear index - used by Kolesik at 527 nm in cm^2/W
        s.mat.n2 = 1.9e-16; %nonlinear index - valid at 800nm in cm^2/W
        
    end

end
%This function is for calculating the linear permittivity of a specified
%medium. Angular frequency s.omg must be passed in in units of rad/s.

function[s] = dispersion_YAPPE(s)

    c = 3e14; %light speed in microns/sec
    lam = 2*pi*c./s.g.omg; %wavelength (must be in microns)

%Here is the dispersion formula one Kolesik uses. It appears as though he
%uses a blend of this and some experimental GVD data. This empirical
%formula is used for angular frequencies greater than 4 PHz, and he uses
%the experimental data from 1 to 4 PHz. From: "Refraction correction for
%fluorescence spectra of squeous solutions".
    
    if s.input.medium == 'water'
        
        %water dispersion formula comes from "Refraction correction for
        %fluorescence spectra of aqueous solutions" - it is the same
        %formula that Kolesik uses for frequencies greater than 4 PHz.
        a1 = 1.7604457; %a1,b1,c1 are dimensionless
        b1 = 4.03368e-3; %a2, b2 and c2 are in micron^2
        c1 = -1.54182e-2;
        d1 = 6.44277e-3;
        d2 = -1.49119e-2;

        s.g.perm = a1 + b1*lam + c1*lam.^2 + d1./(lam.^2 + d2);
        
    elseif s.input.medium == 'SF11'    
    
        %SF11 dispersion comes from Sellmeier formula on
        %refractiveindex.info
        a1 = 1.73759695; %a1,b1,c1 are dimensionless
        a2 = 0.013188707; %a2, b2 and c2 are in micron^2
        b1 = 0.313747356;
        b2 = 0.0623068142;
        c1 = 1.89878101;
        c2 = 155.23629;

        s.g.perm = 1 + a1*lam.^2./(lam.^2 - a2) + b1*lam.^2./(lam.^2 - b2) + c1*lam.^2./(lam.^2 - c2);
        
    elseif s.input.medium == 'vaccuum'
        
        s.g.perm = ones(size(s.g.omg));
        
    end
    
    if s.input.dispersion == 0
    
        %flat dispersion with the permittivity of the central wavelength
        %applied to all wavelengths
        s.g.perm = s.g.perm(1)*ones(size(s.g.omg));
        
    end
    
    %apply absorbing boundary at edges of frequency domain
    if s.input.freqbd == 1
        
        kvac = fftshift(s.g.kvac); %vacuum wavenumbers arranged from lowest to highest
        n_absorb  = zeros(size(s.g.perm));
        width = round(s.input.freqbd_width*length(s.g.omg)); %length in pixels for boundary
        n_amp_low = 0.5/(kvac(1)*s.input.freqbd_length); %imaginary index amplitude for low frequencies
        n_amp_high = 0.5/(kvac(end)*s.input.freqbd_length); %imaginary index amplitude of high frequencies
        
        n_absorb(1:width) = n_amp_low*( 1+cos((0:width-1)*pi/width) ); %cosine shaped function to turn on absorption 
        n_absorb(end-width+1:end) = n_amp_high/n_amp_low*flip(n_absorb(1:width), 2);
        n_absorb = ifftshift(n_absorb); %ifftshift to align with omg/perm axes
        
        s.g.perm = ( sqrt(s.g.perm) + 1i*n_absorb ).^2; %modify permittivity to yield imaginary index
    
    end
    
    
    
end
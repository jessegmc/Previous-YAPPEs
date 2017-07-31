%This function resamples the input radial vector (s.input.r_in) into a
%radial vector (s.g.r) whose spacings are a scaled version of the zeros of
%the j_0 bessel function (necessary for Discrete Hankel Transformation). In
%addition we define the associated transverse wavenumbers s.g.krange.
function [s] = hankel_sample_YAPPE(s)

    bessj0 = @(x) besselj(0,x); %bessel j_0 function handle

    a = zeros(s.input.r_pts+1,1);
    for n = 1:(s.input.r_pts+1)
        a(n) = fzero(bessj0,[(n-1) n]*pi); %find the firts s.pts+1 zeros of j_0
    end

    s.g.aM = a(end); %this is the last 0 
    s.g.bzero = a(1:end-1); %all but last zero
    s.g.r = s.g.r_extent*s.g.bzero/s.g.aM; %radial axis sampled at j_0 zeros
    temp = [0; s.g.r];
    s.g.dr = diff(temp); %radial spacings
    
    %calculate the wavevectors of the nonlinearly sampled radial axis
    Rmax = s.g.r(1)*s.g.aM/s.g.bzero(1);
    s.g.krange = s.g.bzero/Rmax;    

end
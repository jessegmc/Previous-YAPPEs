%This function generates the discrete Hankel transform matrix.

function s = hankel_gen_YAPPE(s)

    temp = bsxfun(@times,s.g.bzero,s.g.bzero');
    s.f.H = bsxfun(@times,(2/s.g.aM)*besselj(0,temp/s.g.aM), 1./(besselj(1,s.g.bzero')).^2);
    
end
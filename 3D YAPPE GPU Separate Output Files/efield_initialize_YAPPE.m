%3d

%This function resamples the input electric field (s.input.E_in) from the
%input radial and local time axes (s.input.r_in and s.input.xi_in) onto the
%radial and local time axes used by the code in the main computation loop
%(s.g.r and s.g.xi).

function [s] = efield_initialize_YAPPE(s)
    
    %generate/load in linearly sampled input
    if s.input.infield.type == 'gauss'
        
        s = gauss_maker_YAPPE(s); %this makes a linearly sampled gaussian beam
        
        
    elseif s.input.infield.type == 'custom'
        
        v = load(s.input.infield.path); %this loads in some time domain electric field
        s.input.E_in = v.E_in;      
% %         s.input.r_in = v.r_in;
        s.input.x_in = v.x_in; 
        s.input.y_in = v.y_in;
        s.input.xi_in = v.xi_in;
        
        x_in_mat = bsxfun(@times,s.input.x_in,ones(1,max(size(s.input.y_in)),max(size(s.input.xi_in))));
        y_in_mat = bsxfun(@times,s.input.y_in,ones(max(size(s.input.x_in)),1,max(size(s.input.xi_in))));   
        xi_in_mat = bsxfun(@times,s.input.xi_in,ones(max(size(s.input.x_in)),max(size(s.input.y_in))),1);
        s.f.E = interp3(xi_in_mat,x_in_mat,y_in_mat,s.input.E_in,s.g.xi,s.g.x,s.g.y);
    end
    s.f.E = s.input.E_in;

    %also initialize the spectral electric field
    s.f.Ef = fftn(s.f.E);
end
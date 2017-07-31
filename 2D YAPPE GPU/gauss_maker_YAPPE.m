%This function creates a gaussian with a specified transverse size (waist),
%pulse length (tau), energy (energ) and quadratic focusing (f), and medium
%wavelength (kmed). If you don't want focusing, enter f= inf. The gaussian
%is defined on the grid r, xi, 

function[s] = gauss_maker_YAPPE(s)
    
    %create radial and temporal axes - (oversample)
    s.input.xi_in = linspace(0, s.input.xi_extent, 10*s.input.xi_pts); %local time axis in seconds
    s.input.r_in = linspace(0, s.input.r_extent, 10*s.input.r_pts)'; %local time axis in seconds
    dxi = s.input.xi_in(2)-s.input.xi_in(1);
    dr = s.input.r_in(2) - s.input.r_in(1);
    
    %gaussian time constant
    tau = s.input.infield.tfwhm /sqrt(2*log(2));  

    %create the desired gaussian shape
    Er(:,1) = exp(-(s.input.r_in/s.input.infield.waist).^2).*exp(-1i*s.g.k(1)*s.input.r_in.^2/(2*s.input.infield.f));
    Et(1,:) = exp(-((s.input.xi_in - s.input.xi_in(end)/2)/tau).^2);
    E = bsxfun(@times,Er,Et); 
    
    %normalize gaussian to create desired energy (assumes I = abs(E.^2) and integral of I = energy)
    I = abs(E.^2);
    nom_energ = 2*pi*dr*dxi*sum(sum(bsxfun(@times,I,s.input.r_in)));

    s.input.E_in = sqrt(s.input.infield.energ/nom_energ)*E;
    
end
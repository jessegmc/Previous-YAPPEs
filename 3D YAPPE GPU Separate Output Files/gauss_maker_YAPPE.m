%3d

%This function creates a gaussian with a specified transverse size (waist),
%pulse length (tau), energy (energ) and quadratic focusing (f), and medium
%wavelength (kmed). If you don't want focusing, enter f= inf. The gaussian
%is defined on the grid r, xi, 

function[s] = gauss_maker_YAPPE(s)
    

    %create tranverse and temporal axes - (oversample)
%     s.input.xi_in = permute(linspace(0, s.input.xi_extent, 1*s.input.xi_pts),[1,3,2]); %local time axis in seconds
%     s.input.x_in = linspace(0, s.input.x_extent, 1*s.input.x_pts);
%     s.input.y_in = linspace(0, s.input.y_extent, 1*s.input.y_pts)'; 
%     dxi = s.input.xi_in(2) - s.input.xi_in(1);
%     dx = s.input.x_in(2) - s.input.x_in(1);
%     dy = s.input.y_in(2) - s.input.y_in(1);
    
    %gaussian time constant
    tau = s.input.infield.tfwhm /sqrt(2*log(2));  

    %create the desired gaussian shape

    Ex(:,1,1) = exp(-(s.g.x/s.input.infield.waistx).^2).*exp(-1i*s.g.k(1)*s.g.x.^2/(2*s.input.infield.f));
    Ey(1,:,1) = exp(-(s.g.y/s.input.infield.waisty).^2).*exp(-1i*s.g.k(1)*s.g.y.^2/(2*s.input.infield.f));
    Et(1,1,:) = exp(-((s.g.xi - s.g.xi(end)/2)/tau).^2);

%     Ex(:,1,1) = exp(-(s.input.x_in/s.input.infield.waistx).^2).*exp(-1i*s.g.k(1)*s.g.x.^2/(2*s.input.infield.f));
%     Ey(1,:,1) = exp(-(s.input.y_in/s.input.infield.waisty).^2).*exp(-1i*s.g.k(1)*s.g.y.^2/(2*s.input.infield.f));
%     Et(1,1,:) = exp(-((s.input.xi_in - s.input.xi_in(end)/2)/tau).^2);
    
    E = bsxfun(@times, Et, bsxfun(@times,Ex,Ey));
    
    %normalize gaussian to create desired energy (assumes I = abs(E.^2) and integral of I = energy)

    I = abs(E.^2);
    nom_energ = s.g.dx*s.g.dy*s.g.dxi*sum(sum(sum(I)));
%     nom_energ = dx*dy*dxi*sum(sum(sum(I)));

    s.input.E_in = sqrt(s.input.infield.energ/nom_energ)*E;

end
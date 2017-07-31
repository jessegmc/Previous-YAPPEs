%3d

%This function creates a gaussian with a specified transverse size (waist),
%pulse length (tau), energy (energ) and quadratic focusing (f), and medium
%wavelength (kmed). If you don't want focusing, enter f= inf. The gaussian
%is defined on the grid r, xi,

function[s] = double_gauss_maker_YAPPE(s)


%gaussian time constant
tau = s.input.infield.tfwhm /sqrt(2*log(2));

%create the desired gaussian shape

Ex1(:,1,1) = exp(-((s.g.x-s.input.infield.beamOffset)/s.input.infield.waistx).^2).*exp(-1i*s.g.k(1)*s.g.x.^2/(2*s.input.infield.f));
Ey1(1,:,1) = exp(-(s.g.y/s.input.infield.waisty).^2).*exp(-1i*s.g.k(1)*s.g.y.^2/(2*s.input.infield.f));
Et1(1,1,:) = exp(-((s.g.xi - s.g.xi(end)/2)/tau).^2);

Ex2(:,1,1) = exp(-((s.g.x+s.input.infield.beamOffset)/s.input.infield.waistx).^2).*exp(-1i*s.g.k(1)*s.g.x.^2/(2*s.input.infield.f));
Ey2(1,:,1) = exp(-(s.g.y/s.input.infield.waisty).^2).*exp(-1i*s.g.k(1)*s.g.y.^2/(2*s.input.infield.f));
Et2(1,1,:) = exp(-((s.g.xi - s.g.xi(end)/2)/tau).^2);

%     Ex(:,1,1) = exp(-(s.input.x_in/s.input.infield.waistx).^2).*exp(-1i*s.g.k(1)*s.g.x.^2/(2*s.input.infield.f));
%     Ey(1,:,1) = exp(-(s.input.y_in/s.input.infield.waisty).^2).*exp(-1i*s.g.k(1)*s.g.y.^2/(2*s.input.infield.f));
%     Et(1,1,:) = exp(-((s.input.xi_in - s.input.xi_in(end)/2)/tau).^2);

E1 = bsxfun(@times, Et1, bsxfun(@times,Ex1,Ey1));

E2 = bsxfun(@times, Et2, bsxfun(@times,Ex2,Ey2));



%
E = E1+E2;
%     E = E1;
%normalize gaussian to create desired energy (assumes I = abs(E.^2) and integral of I = energy)

I = abs(E.^2);
nom_energ = s.g.dx*s.g.dy*s.g.dxi*sum(sum(sum(I)));
%     nom_energ = dx*dy*dxi*sum(sum(sum(I)));

s.input.E_in = sqrt(s.input.infield.energ/nom_energ)*E;

end
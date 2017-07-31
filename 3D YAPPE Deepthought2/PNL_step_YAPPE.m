%This function steps the solution forward w.r.t. the nonlinear
%polarizability. This is where almost all the computation time is. The load
%is fairly equally spread between 2 fft calls, 2 hankel transforms, the
%plasma integration (vector sums), and the 2 element-wise matrix
%exponentiation calls. Declaring the struct "s" as global does not have a
%performance cost. I observe roughly equal call times between the above
%operations for a 300x300 grid.

function[dEf_dz_vec] = PNL_step_YAPPE(z,Ef_vec)

    global s
    
    %reshape Ef back into a matrix
% %   Ef = reshape(Ef_vec,s.input.r_pts,s.input.xi_pts);
    Ef = reshape(Ef_vec,s.input.x_pts,s.input.y_pts,s.input.xi_pts);
    Ef = Ef.*exp(1i*s.f.Kz_move*z);
    
    %convert to spatiotemporal domain
% %     E = ifft(Ef,[],2);
% %     E = s.f.H*E;

    E = ifftn(Ef);
    
    %calculate intensity envelope
    s.f.I = abs(E).^2;
    Ipow = s.f.I.^(s.mat.pow);
    
    %calculate plasma density in #/cm^3
    s.f.rho = zeros(size(s.f.I));
 
     if s.input.plasma == 1 %this toggles the plasma module        
% %         for m = 2:size(s.f.rho,2)
% %             s.f.rho(:,m) = s.f.rho(:,m-1) + s.g.dxi*( s.ion.a1*s.f.I(:,m-1).*s.f.rho(:,m-1) + s.ion.a2*Ipow(:,m-1) + s.ion.a3*s.f.rho(:,m-1).^2 );           
% %         end
        for m = 2:size(s.f.rho,3)
            s.f.rho(:,:,m) = s.f.rho(:,:,m-1) + s.g.dxi*( s.ion.a1*s.f.I(:,:,m-1).*s.f.rho(:,:,m-1) + s.ion.a2*Ipow(:,:,m-1) + s.ion.a3*s.f.rho(:,:,m-1).^2 );
        end
     end
     clear Ipow;

    %calculate nonlinear susceptibility
    s.f.chiNL = s.input.n2*( s.NL.b1*s.f.I ) + s.input.plasma*( s.NL.b2*s.f.rho + s.NL.b3*s.f.I.^(s.mat.pow-1) );
    clear s.f.rho;
    
    %calculate nonlinear polarizability
    s.f.PNL = s.f.chiNL.*E; 
    clear s.f.chiNL;
    
    %convert nonlinear polarizability to spectral domain
% %     s.f.PNLf = fft(s.f.PNL,[],2);
% %     s.f.PNLf = s.f.H*s.f.PNLf;

    s.f.PNLf = fftn(s.f.PNL);
    clear s.f.PNL;
     
    %z derivative in matrix form
    dEf_dz =  0.5*1i*s.f.Q.*s.f.PNLf.*exp(-1i*s.f.Kz_move*z);
     
    %z derivative in vector form
    dEf_dz_vec = dEf_dz(:);
    s.count = s.count+1;
    
end

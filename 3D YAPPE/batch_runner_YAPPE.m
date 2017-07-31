%3d

% For use in a for loop. Define a global struct s constructed from the
% input deck, and this function will run YAPPE for the specified initial
% conditions. The struct s has to be global so that the matlab ODE function
% call can use information from s.

function[] = batch_runner_YAPPE()
    
    %%%%%% INITIALIZATIONS %%%%%%%%%

    global s
    s = fund_const_YAPPE(s); %fundamental constants
    s = medium_property_YAPPE(s); %medium parameters
    s = grids_YAPPE(s); %make grids and interpolate input E-field onto new grids
    s = precompute_YAPPE(s); %precompute useful arrays (Hankel matrix, Kz, Q, various coefficiencts...)
    s = efield_initialize_YAPPE(s); %generate input E-fields (spatiotemporal and spectral)     
    
    %initialize outputs (helps with speed and acts as a check that there is enough memory for outputs)
% %     s.f.E_out = zeros(s.input.r_pts,s.input.xi_pts,length(s.g.zout)); %initialize output E-field matrix
% %     s.f.E_out(:,:,1) = s.f.E; %first entry is input field    
% %     s.f.rho_out_r = zeros(s.input.r_pts,length(s.g.zout)); %initialize output plasma matrix
% %     s.f.rho_out_xi = zeros(s.input.xi_pts,length(s.g.zout)); %initialize output plasma matrix

    s.f.E_out = zeros(s.input.x_pts,s.input.y_pts,s.input.xi_pts,length(s.g.zout)); %initialize output E-field matrix
    s.f.E_out(:,:,:,1) = s.f.E; %first entry is input field    
    s.f.rho_out_x = zeros(s.input.x_pts,length(s.g.zout)); %initialize output plasma matrix
    s.f.rho_out_y = zeros(s.input.y_pts,length(s.g.zout)); %initialize output plasma matrix
    s.f.rho_out_xi = zeros(s.input.xi_pts,length(s.g.zout)); %initialize output plasma matrix
    
    s.count = 1; %used for monitoring # of ode45 function calls
    opts = odeset('RelTol',s.input.RelTol,'AbsTol',s.input.AbsTol); %set solution tolerances

    %%%%%% MAIN COMPUTATION LOOP %%%%%%%%%
    
    for m = 2:length(s.g.zout)

        initial_Ef = s.f.Ef(:); %convert e-field matrix into vector for ode45 solve
        [~,x] = ode45( @(t,y) PNL_step_YAPPE(t,y), [0 s.input.outperiod/2 s.input.outperiod], initial_Ef, opts); %solve ODE

% %         s.f.Ef = reshape(x(end,:)', [s.input.r_pts, s.input.xi_pts]); %convert output of ODE solve to matrix          
        s.f.Ef = reshape(x(end,:,:)', [s.input.x_pts, s.input.y_pts, s.input.xi_pts]); %convert output of ODE solve to matrix        
        clear x
        s.f.Ef = conj(s.f.Ef); %ODE45 has a "dagger" where they should have a "transpose", must conjugate
        s.f.Ef = s.f.Ef.*s.f.lin_prop; %apply linear propagator "realignment"
        
        %convert from spectral to spatiotemporal and write to output array
% %         s.f.E_out(:,:,m) = ifft(s.f.Ef,[],2);
% %         s.f.E_out(:,:,m) = s.f.H*s.f.E_out(:,:,m);


            s.f.E_out(:,:,:,m) = ifftn(s.f.Ef);

        
        %write out the the local time-integrated plasma
% %         s.f.rho_out_r(:,m) = squeeze(sum(s.f.rho*s.g.dxi,2)); %radially resolved time integrated plasma density
% %         s.f.rho_out_xi(:,m) = squeeze(s.f.rho(1,:)); %axial lineout of plasma density

        %Haven't changed this part yet - Nihal says the physics meaning is
        %nonsense, return to this later
% %         s.f.rho_out_x(:,m) = squeeze(sum(s.f.rho*s.g.dxi,2)); %x time integrated plasma density
% %         s.f.rho_out_y(:,m) = squeeze(sum(s.f.rho*s.g.dxi,2)); %y time integrated plasma density
% %         s.f.rho_out_xi(:,m) = squeeze(s.f.rho(1,:)); %axial lineout of plasma density
        
        disp(strcat( 'propagated to z =', 32, num2str(s.g.zout(m)), 32, 'cm'))

    end
    
    
    %%%%%% WRITING OUT %%%%%%%%%
    
    %save the outputs configuration
    outputnam = strcat(s.input.outpath,'full_output.mat');
    if exist(s.input.outpath,'file')==0
        mkdir(s.input.outpath);
    end   
    cd(s.input.outpath);
    save(outputnam, 's');

    disp('a million coupled ODEs cried out... and were solved')
    
end
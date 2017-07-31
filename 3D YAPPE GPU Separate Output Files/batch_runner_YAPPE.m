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
    
    %%%%%% SAVING GRIDS %%%%%%%%
    if (exist(s.input.outpath,'file')==0)
        mkdir(s.input.outpath);
    end
    input_outpath = strcat(s.input.outpath,'input.mat');
    save(input_outpath,'s');
    
    s = precompute_YAPPE(s); %precompute useful arrays (Hankel matrix, Kz, Q, various coefficiencts...)
    s = efield_initialize_YAPPE(s); %generate input E-fields (spatiotemporal and spectral)     
 
    %%%%%% SAVING FIRST OUTPUT %%%%%%
    E_out = s.f.E;
    Eout_general_path = strcat(s.input.outpath,'Outputs/');
    if (exist(Eout_general_path,'file')==0)
        mkdir(Eout_general_path);
    end
    Eout_specific_path = strcat(Eout_general_path, 'Step',num2str(1));
    save(Eout_specific_path,'E_out');
     
    s.count = 1; %used for monitoring # of ode45 function calls
    opts = odeset('RelTol',s.input.RelTol,'AbsTol',s.input.AbsTol); %set solution tolerances
    
    disp('Entering main computation loop')
    %%%%%% MAIN COMPUTATION LOOP %%%%%%%%%
    
    for m = 2:length(s.g.zout)
        initial_Ef = s.f.Ef(:); %convert e-field matrix into vector for ode45 solve
        [~,x] = ode45( @(t,y) PNL_step_YAPPE(t,y), [0 s.input.outperiod/2 s.input.outperiod], initial_Ef, opts); %solve ODE
        s.f.Ef = reshape(x(end,:,:)', [s.input.x_pts, s.input.y_pts, s.input.xi_pts]); %convert output of ODE solve to matrix        
        clear x
        s.f.Ef = conj(s.f.Ef); %ODE45 has a "dagger" where they should have a "transpose", must conjugate
        s.f.Ef = s.f.Ef.*s.f.lin_prop; %apply linear propagator "realignment"
        
        disp('saving output');
        %%%%%% SAVING EACH OUTPUT %%%%%%%
        E_out = ifftn(s.f.Ef);
        Eout_specific_path = strcat(Eout_general_path, 'Step',num2str(m));
        save(Eout_specific_path, 'E_out');
        
        disp(strcat( 'propagated to z =', 32, num2str(s.g.zout(m)), 32, 'cm'))
    end

    disp('a million coupled ODEs cried out... and were solved')
    
end
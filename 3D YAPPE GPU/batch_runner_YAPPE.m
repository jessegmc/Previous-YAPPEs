%3d

% For use in a for loop. Define a global struct s constructed from the
% input deck, and this function will run YAPPE for the specified initial
% conditions. The struct s has to be global so that the matlab ODE function
% call can use information from s.

function[] = batch_runner_YAPPE()

%%%%%% INITIALIZATIONS %%%%%%%%%

global s
s.gpu = gpuDevice(1);
reset(s.gpu);
s = fund_const_YAPPE(s); %fundamental constants
s = medium_property_YAPPE(s); %medium parameters
s = grids_YAPPE(s); %make grids and interpolate input E-field onto new grids
s = precompute_YAPPE(s); %precompute useful arrays (Hankel matrix, Kz, Q, various coefficiencts...)
s = efield_initialize_YAPPE(s); %generate input E-fields (spatiotemporal and spectral)

disp('Precompute finished... entering ode45')

% s.f.E_out = zeros(s.input.x_pts,s.input.xi_pts,length(s.g.zout)); %initialize output E-field matrix
s.f.E_out = zeros(s.input.x_pts,s.input.y_pts,s.input.xi_pts,length(s.g.zout)); %initialize output E-field matrix
% s.f.E_out(:,:,1) = s.f.E(:,round(size(s.f.E,2)/2),:); %first entry is input field
s.f.E_out(:,:,:,1) = s.f.E; %first entry is input field

s.count = 1; %used for monitoring # of ode45 function calls
opts = odeset('RelTol',s.input.RelTol,'AbsTol',s.input.AbsTol); %set solution tolerances

%%%%%% MAIN COMPUTATION LOOP %%%%%%%%%

for m = 2:length(s.g.zout)
    
    [~,s.f.Ef] = ode45GPU( @(t,y) PNL_step_YAPPE(t,y), [0 s.input.outperiod/2 s.input.outperiod], s.f.Ef(:), opts); %solve ODE
    s.f.Ef = reshape(s.f.Ef', [s.input.x_pts, s.input.y_pts, s.input.xi_pts]); %convert output of ODE solve to matrix
    s.f.Ef = conj(s.f.Ef); %ODE45 has a "dagger" where they should have a "transpose", must conjugate
    s.f.Ef = s.f.Ef.*s.f.lin_prop; %apply linear propagator "realignment"
    s.f.E = gather(ifftn(s.f.Ef));
    s.f.E_out(:,:,:,m) = s.f.E;

    disp(strcat( 'propagated to z =', 32, num2str(s.g.zout(m)), 32, 'cm'))
end


%%%%%% WRITING OUT %%%%%%%%%

%save the outputs configuration
outputnam = strcat(s.input.outpath,'full_output.mat');
inputnam = strcat(s.input.outpath,'input_deck.mat');
if exist(s.input.outpath,'file')==0
    mkdir(s.input.outpath);
end
cd(s.input.outpath);
s.f = rmfield(s.f,'Kz');
s.f = rmfield(s.f,'Q');
s.f = rmfield(s.f,'Kz_move');
s.f = rmfield(s.f,'lin_prop');
s.f = rmfield(s.f,'E');
s.f = rmfield(s.f,'Ef');
s.f = rmfield(s.f,'I');
save(outputnam, 's','-v7.3');
input = s.input;
% save(inputnam, 'input','-v7.3');

disp('a million coupled ODEs cried out... and were solved')

end
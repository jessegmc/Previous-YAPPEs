function [tf,ynew] = ode45GPU(ode,tspan,y,options,varargin)

% Check inputs
if nargin < 4
    options = [];
    if nargin < 3
        y = [];
        if nargin < 2
            tspan = [];
            if nargin < 1
                error(message('MATLAB:ode45:NotEnoughInputs'));
            end
        end
    end
end

%Stats
nsteps = 0;
nfailed = 0;
nfevals = 0;

%Output
FcnHandlesUsed  = isa(ode,'function_handle');
output_ty  = 1;  % [t,y,...] = odeXX(...)
% There might be no output requested...


% Handle solver arguments - odearguments
if FcnHandlesUsed  % function handles used
    if isempty(tspan) || isempty(y)
        error(message('MATLAB:odearguments:TspanOrY0NotSupplied', solver));
    end
    if length(tspan) < 2
        error(message('MATLAB:odearguments:SizeTspan', solver));
    end
    htspan = abs(tspan(2) - tspan(1));
    tspan = tspan(:);
    ntspan = length(tspan);
    t0 = tspan(1);
    next = 2;       % next entry in tspan
    tfinal = tspan(end);
    args = {};
    
else  % ode-file used   (ignored when solver == ODE15I)
    % Get default tspan and y0 from the function if none are specified.
    if isempty(tspan) || isempty(y)
        if exist(ode)==2 && ( nargout(ode)<3 && nargout(ode)~=-1 )
            error(message('MATLAB:odearguments:NoDefaultParams', funstring( ode ), solver, funstring( ode )));
        end
        [def_tspan,def_y0,def_options] = feval(ode,[],[],'init',extras{:});
        if isempty(tspan)
            tspan = def_tspan;
        end
        if isempty(y)
            y = def_y0;
        end
        options = odeset(def_options,options);
    end
    tspan = tspan(:);
    ntspan = length(tspan);
    if ntspan == 1    % Integrate from 0 to tspan
        t0 = 0;
        next = 1;       % Next entry in tspan.
    else
        t0 = tspan(1);
        next = 2;       % next entry in tspan
    end
    htspan = abs(tspan(next) - t0);
    tfinal = tspan(end);
    
    % The input arguments of f determine the args to use to evaluate f.
    if (exist(ode)==2)
        if (nargin(ode) == 2)
            args = {};                   % f(t,y)
        else
            args = [{''} extras];        % f(t,y,'',p1,p2...)
        end
    else  % MEX-files, etc.
        try
            args = [{''} extras];        % try f(t,y,'',p1,p2...)
            feval(ode,tspan(1),y(:),args{:});
        catch ME
            args = {};                   % use f(t,y) only
        end
    end
end

y = y(:);
neq = length(y);

% Test that tspan is internally consistent.
if t0 == tfinal
    error(message('MATLAB:odearguments:TspanEndpointsNotDistinct'));
end
tdir = sign(tfinal - t0);
if any( tdir*diff(tspan) <= 0 )
    error(message('MATLAB:odearguments:TspanNotMonotonic'));
end

f = gpuArray(zeros(neq,7,'double'));
% f = zeros(neq,7,'double');

f(:,1) = feval(ode,t0,y,args{:});   % ODE15I sets args{1} to yp0.;
% f(:,1) = gather(feval(ode,t0,y,args{:}));   % ODE15I sets args{1} to yp0.;
[m,n] = size(f(:,1));
if n > 1
    error(message('MATLAB:odearguments:FoMustReturnCol', funstring( ode )));
elseif m ~= neq
    error(message('MATLAB:odearguments:SizeIC', funstring( ode ), m, neq, funstring( ode )));
end

dataType = 'double';

% Get the error control options, and set defaults.
rtol = odeget(options,'RelTol',1e-3,'fast');
if (length(rtol) ~= 1) || (rtol <= 0)
    error(message('MATLAB:odearguments:RelTolNotPosScalar'));
end
if rtol < 100 * eps(dataType)
    rtol = 100 * eps(dataType);
    warning(message('MATLAB:odearguments:RelTolIncrease', sprintf( '%g', rtol )))
end
atol = odeget(options,'AbsTol',1e-6,'fast');
if any(atol <= 0)
    error(message('MATLAB:odearguments:AbsTolNotPos'));
end
normcontrol = strcmp(odeget(options,'NormControl','off','fast'),'on');
if normcontrol
    if length(atol) ~= 1
        error(message('MATLAB:odearguments:NonScalarAbsTol'));
    end
    normy = norm(y);
else
    if (length(atol) ~= 1) && (length(atol) ~= neq)
        error(message('MATLAB:odearguments:SizeAbsTol', funstring( ode ), neq));
    end
    atol = atol(:);
    normy = [];
end
threshold = atol / rtol;

% By default, hmax is 1/10 of the interval.
hmax = min(abs(tfinal-t0), abs(odeget(options,'MaxStep',0.1*(tfinal-t0),'fast')));
if hmax <= 0
    error(message('MATLAB:odearguments:MaxStepLEzero'));
end
htry = abs(odeget(options,'InitialStep',[],'fast'));
if ~isempty(htry) && (htry <= 0)
    error(message('MATLAB:odearguments:InitialStepLEzero'));
end

odeFcn = ode;


% Handle the output
if nargout > 0
    outputFcn = odeget(options,'OutputFcn',[],'fast');
else
    outputFcn = odeget(options,'OutputFcn',@odeplot,'fast');
end
outputArgs = {};
if isempty(outputFcn)
    haveOutputFcn = false;
else
    haveOutputFcn = true;
    outputs = odeget(options,'OutputSel',1:neq,'fast');
    if isa(outputFcn,'function_handle')
        % With MATLAB 6 syntax pass additional input arguments to outputFcn.
        outputArgs = varargin;
    end
end


%check refined points
refine = max(1,odeget(options,'Refine',4,'fast'));
if ntspan > 2
    outputAt = 'RequestedPoints';         % output only at tspan points
elseif refine <= 1
    outputAt = 'SolverSteps';             % computed points, no refinement
else
    outputAt = 'RefinedSteps';            % computed points, with refinement
    S = (1:refine-1) / refine;
end
printstats = strcmp(odeget(options,'Stats','off','fast'),'on');

% No event functions allowed

% No mass matrices allowed

% Non-negative solution components
idxNonNegative = odeget(options,'NonNegative',[],'fast');
nonNegative = ~isempty(idxNonNegative);
if nonNegative  % modify the derivative function
    [odeFcn,thresholdNonNegative] = odenonnegative(odeFcn,y,threshold,idxNonNegative);
    f(:,1) = feval(odeFcn,t0,y,odeArgs{:});
%     f(:,1) = gather(feval(odeFcn,t0,y,odeArgs{:}));
    nfevals = nfevals + 1;
end


%Evaluate first point
t = t0;
% y = y0;

%Allocate memory for output
nout = 0;
tout = []; 
% yout = [];
if ntspan>2
    tout = zeros(1,ntspan,dataType);
else
    chunk = min(max(100,50*refine), refine+floor((2^13)/neq));
    tout = zeros(1,chunk,dataType);
    yout = gpuArray(zeros(neq,chunk,dataType));
end

nout = 1;
tout(nout) = t;
% yout(:,nout) = y;

%Initialize method parameters
pow = 1/5;
A = [1/5, 3/10, 4/5, 8/9, 1, 1];
B = [
    1/5         3/40    44/45   19372/6561      9017/3168       35/384
    0           9/40    -56/15  -25360/2187     -355/33         0
    0           0       32/9    64448/6561      46732/5247      500/1113
    0           0       0       -212/729        49/176          125/192
    0           0       0       0               -5103/18656     -2187/6784
    0           0       0       0               0               11/84
    0           0       0       0               0               0
    ];
E = [71/57600; 0; -71/16695; 71/1920; -17253/339200; 22/525; -1/40];

hmin = 16*eps(t);
if isempty(htry)
    % Compute an initial step size h using y'(t).
    absh = min(hmax, htspan);
    if normcontrol
        rh = (norm(f(:,1)) / max(normy,threshold)) / (0.8 * rtol^pow);
    else
        rh = norm(f(:,1) ./ max(abs(y),threshold),inf) / (0.8 * rtol^pow);
    end
    if absh * rh > 1
        absh = 1 / rh;
    end
    absh = max(absh, hmin);
else
    absh = min(hmax, max(hmin, htry));
end

% THE MAIN LOOP

done = false;
while ~done
    %Calcualte new hmin based on t
    hmin = 16*eps(t);
    absh = min(hmax, max(hmin, absh));    % couldn't limit absh until new hmin
    h = tdir * absh;
    
    if 1.1*absh >= abs(tfinal - t)
        h = tfinal - t;
        absh = abs(h);
        done = true;
    end
    
    % LOOP FOR ADVANCING ONE STEP.
    nofailed = true;    %No failed attempts yet
    while true
        hA = h * A;
        hB = h * B;
        for i = 2:6
            f(:,i) = feval(odeFcn,t+hA(i-1),y+f*hB(:,i-1));
%           f(:,i) = gather(feval(odeFcn,t+hA(i-1),y+f*hB(:,i-1)));
        end
        
        tnew = t + hA(6);
        if done
            tnew = tfinal;
        end
        h = tnew-t; %Create new h
        
        ynew = y + f*hB(:,6);
        f(:,7) = feval(odeFcn,tnew,ynew);
%         f(:,7) = gather(feval(odeFcn,tnew,ynew));
        nfevals = nfevals+6;
        
        % Estimate the error.
        NNrejectStep = false;
        if normcontrol
            normynew = norm(ynew);
            errwt = max(max(normy,normynew),threshold);
            err = absh * (norm(f * E) / errwt);
            if nonNegative && (err <= rtol) && any(ynew(idxNonNegative)<0)
                errNN = norm( max(0,-ynew(idxNonNegative)) ) / errwt ;
                if errNN > rtol
                    err = errNN;
                    NNrejectStep = true;
                end
            end
        else
            err = absh * norm((f * E) ./ max(max(abs(y),abs(ynew)),threshold),inf);
            if nonNegative && (err <= rtol) && any(ynew(idxNonNegative)<0)
                errNN = norm( max(0,-ynew(idxNonNegative)) ./ thresholdNonNegative, inf);
                if errNN > rtol
                    err = errNN;
                    NNrejectStep = true;
                end
            end
        end
        
        % Accept the solution only if the weighted error is no more than the
        % tolerance rtol.  Estimate an h that will yield an error of rtol on
        % the next step or the next try at taking this step, as the case may be,
        % and use 0.8 of this value to avoid failures.
        if err > rtol     %Failed step
            if absh <= hmin
                warning(message('MATLAB:ode45:IntegrationTolNotMet', sprintf( '%e', t ), sprintf( '%e', hmin )));
                return;
            end
            
            if nofailed
                nofailed = false;
                if NNrejectStep
                    absh = max(hmin, 0.5*absh);
                else
                    absh = max(hmin, absh * max(0.1, 0.8*(rtol/err)^pow));
                end
            else
                absh = max(hmin, 0.5 * absh);
            end
            h = tdir *absh;
            done = false;
        else  %Successful step
            NNreset_f7 = false;
            if nonNegative && any(ynew(idxNonNegative)<0)
                ynew(idxNonNegative) = max(ynew(idxNonNegative),0);
                if normcontrol
                    normynew = norm(ynew);
                end
                NNreset_f7 = true;
            end
            
            break;
        end
    end
    nsteps = nsteps + 1;
    
    if output_ty || haveOutputFcn
        switch outputAt
            case 'SolverSteps'        % computed points, no refinement
                nout_new = 1;
                tout_new = tnew;
                yout_new = ynew;
            case 'RefinedSteps'       % computed points, with refinement
                tref = t + (tnew-t)*S;
                nout_new = refine;
                tout_new = [tref, tnew];
                %yout_new = [ntrp45(tref,t,y,[],[],h,f,[]), ynew];
                %ntrp45
                BI = [
                    1       -183/64      37/12       -145/128
                    0          0           0            0
                    0       1500/371    -1000/159    1000/371
                    0       -125/32       125/12     -375/64
                    0       9477/3392   -729/106    25515/6784
                    0        -11/7        11/3        -55/28
                    0         3/2         -4            5/2
                    ];
                S = (tref - t)/h;
                yinterp = y(:,ones(size(tref))) + f*(h*BI)*cumprod([S;S;S;S]);
                %end ntrp45
                yout_new = [yinterp, ynew];
            case 'RequestedPoints'    % output only at the final tspan point
%                 nout_new =  0;
%                 tout_new = [];
%                 yout_new = [];
%                 while next <= ntspan
%                     if tdir * (tnew - tspan(next)) < 0
%                         break;
%                     end
%                     nout_new = nout_new + 1;
%                     tout_new = [tout_new, tspan(next)];
%                     if tspan(next) == tnew
%                         yout_new = [yout_new, ynew];
%                     else
%                         %ntrp45
%                         BI = [
%                             1       -183/64      37/12       -145/128
%                             0          0           0            0
%                             0       1500/371    -1000/159    1000/371
%                             0       -125/32       125/12     -375/64
%                             0       9477/3392   -729/106    25515/6784
%                             0        -11/7        11/3        -55/28
%                             0         3/2         -4            5/2
%                             ];
%                         S = (tspan(next) - t)/h;
%                         yinterp = y(:,ones(size(tspan(next)))) + f*(h*BI)*cumprod([S;S;S;S]);
%                         %end ntrp45
%                         yout_new = [yout_new, yinterp];
%                     end
%                     next = next + 1;
%                 end
        end
        
%         if nout_new > 0
%             oldnout = nout;
%             nout = nout + nout_new;
%             if nout > length(tout)
%                 tout = [tout, zeros(1,chunk,dataType)];  % requires chunk >= refine
%                 yout = [yout, gpuArray(zeros(neq,chunk,dataType))];
%             end
%             idx = oldnout+1:nout;
%             tout(idx) = tout_new;
%             yout(:,idx) = yout_new;
%         end
%     end
    
    if done
        break
    end
    
    %If there were no failures, compute new h
    if nofailed
        % Note that absh may shrink by 0.8, and that err may be 0.
        temp = 1.25*(err/rtol)^pow;
        if temp > 0.2
            absh = absh / temp;
        else
            absh = 5.0*absh;
        end
    end
    
    %Advance integration one step
    t = tnew;
    y = ynew;
    if normcontrol
        normy = normynew;
    end
    f(:,1) = f(:,7);
end

tf = tout(1:nout).';
% yf = yout(:,1:nout).';

end
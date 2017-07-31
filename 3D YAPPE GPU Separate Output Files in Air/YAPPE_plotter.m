clear all
close all

dir = '~/YAPPE_outputs/water_cell_convergence_testing';
cd(dir)

runs = [1:10];

for m = 1:length(runs)
    
    %load in matspace
    runstr{runs(m)} = ['r', num2str(runs(m))];
    load([dir,'/run',32, num2str(runs(m)),'/full_output.mat']);
    v.(runstr{runs(m)}) = s;
    clear s;
    
    %reconstruct the spectral field
    v.(runstr{runs(m)}).f.Ef_out = ifft(v.(runstr{runs(m)}).f.E_out, [] ,2);    
    for n = 1:size(v.(runstr{runs(m)}).f.E_out,3)
        v.(runstr{runs(m)}).f.Ef_out(:,:,n) = v.(runstr{runs(m)}).f.H*v.(runstr{runs(m)}).f.Ef_out(:,:,n);
    end            

end

disp('mat spaces are loaded and spectral fields have been constructed')

%%
for m = 1:length(runs)
    legend_cell{m} = strcat('r',num2str(runs(m)));
end

for n = 1:size(v.r1.f.E_out,3)
    figure(1)
    clf
    for m = 1:length(runs)
    
        s = v.(runstr{runs(m)});
        I = abs(s.f.E_out).^2;
        
        figure(1)
        hold all
        plot(s.g.xi,(I(1,:,n)))
        title(v.r1.g.zout(n))
%         pause
        
    end
    legend(legend_cell)
    
    figure(2)
    clf
    for m = 1:length(runs)
    
        s = v.(runstr{runs(m)});
        xi_ind = round(s.g.xi_pts/2);
        I = abs(s.f.E_out).^2;
        
        figure(2)
        hold all
        plot(s.g.r(s.g.r<.05),(I(s.g.r<.05,xi_ind,n)))
        title(s.g.zout(n))
%         pause
        
    end
    legend(legend_cell)
    
    
    
    pause
end
        
    

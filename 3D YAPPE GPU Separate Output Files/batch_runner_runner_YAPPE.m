%Script for calling batch_runner_YAPPE (to actually run the code). In the
%run loop, a master struct is loaded in (input_master.mat) which contains
%the requisite inputs to run batch_runner_YAPPE. Before calling
%batch_runner_YAPPE, various inputs can be modified so that different
%computations can be performed in a loop.

clear all
close all

global s


%generate the master input struct
s = input_deck_YAPPE();
maindir = s.input.outpath;

% % xipts = [1 2 .5 2 .5 1 1 1 1 1];
% % Txi = [1 1 1 2 .5 1 1 1 1 1];
% % 
% % rpts = [1 1 1 1 1 2 .5 2 .5 1];
% % Tr = [1 1 1 1 1 1 1 2 .5 1];0.1
% 
% energ = [1 1.5 2 2.5 3 3.5 4 4.5 5]*1e-6;
% 
% % % run YAPPE in a loop
% runs = 1:length(energ);
% tic;
% for m = runs
% % 
% % %     load('input_master.mat')
%     s.input.infield.energ = energ(m);
% %     
% % %     s.input.xi_pts = xipts(m)*s.input.xi_pts;
% % %     s.input.xi_extent = Txi(m)*s.input.xi_extent;
% % %     
% % %     s.input.r_pts = rpts(m)*s.input.r_pts;
% % %     s.input.r_extent = Tr(m)*s.input.r_extent;
% %     
% % %     if m == 10
% % %         s.input.freqbd = 0;
% % %     end
% %     
%     s.input.outpath = strcat(maindir,'run',32, num2str(m),'/');    
%     batch_runner_YAPPE();
%     m
% %     
% end

% load('input_master.mat')
batch_runner_YAPPE();

h = toc
s.count
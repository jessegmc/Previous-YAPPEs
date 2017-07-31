%Script for calling batch_runner_YAPPE (to actually run the code). In the
%run loop, a master struct is loaded in (input_master.mat) which contains
%the requisite inputs to run batch_runner_YAPPE. Before calling
%batch_runner_YAPPE, various inputs can be modified so that different
%computations can be performed in a loop.
function s = batch_runner_runner_YAPPE

clear
close all

global s

s = input_deck_YAPPE();

batch_runner_YAPPE();


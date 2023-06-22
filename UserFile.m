%%% UserFile
% This is the main interface file for MOMIBB. Please make sure that you 
% have installed and started Gurobi and INTLAB before calling AdEnA.
% Call initSession.m prior to first use

%% Some clean-up first
clear;
close all;
clc;

%% Please enter your parameters below
% Your Problem
problem = 'T4';
param = [2;2];

% Provide initial bounds yourself or leave empty to auto-compute (INTLAB is
% required for auto-compute)
L = [];
U = [];

% Set quality epsilon and offset
EPSILON = 0.1;
OFFSET = EPSILON*1e-3;

% Specify cutting strategy [0 == none, 1 == continuous, 2 == continuous+mixed integer]
CUT_MODE = 1;

% Specify handling of quadratic instances [1 == normal, 2 == directly using Gurobi]
SOL_MODE = 1;

% Should the result be plotted (m = 2 and m = 3 only) [1 == yes, 0 == no]
plot_result = 1;

%% Call solver
[L,U,N,box_struct,it,exitflag,time] = callSolver(problem,param,L,U,CUT_MODE,SOL_MODE,EPSILON,OFFSET,plot_result);
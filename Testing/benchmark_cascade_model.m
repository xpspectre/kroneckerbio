%% Benchmark cascade model building and simulation
clear; close all; clc
rng('default');

N = 100;

%% Build
tic
m = cascade_model(N);
disp(['Time to build: ' num2str(toc)])

%% Simulate
% Construct experiment
con = experimentInitialValue(m, [], [], [], 'CascadeExpt');
tF = 1; % final time

tic
sim = SimulateSystem(m, con, tF);
disp(['Time to simulate: ' num2str(toc)])

% Plot result
figure
times = linspace(0, tF, 100);
plot(times, sim.x(times))
xlabel('Time')
ylabel('Amount')
title('Cascade Model Simulation')

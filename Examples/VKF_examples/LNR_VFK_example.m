%% LNR example
% This example employs the model lnr_hgf
% The first section simulates beliefs and responses of a putative agent.
% The second section shows the use of PAM
%% Initialize path
clear
addpath(genpath(fullfile('..','..','PAM_master')))

load u.mat

%% SECTION 1: Simulate dataset

% TODO: Set parameters ------------------------------------------------
lambda = .5; % learning rate
v_0 = .5;
omega = .1;
a = -0.53; % mu of the lognormal distribution
b_val = -0.47; % validity effect on a
b = -0.19; % belief effect
sigma = .25; % LNR variance
Ter = 0; % non-decision time
%----------------------------------------------------------------------

% SIMULATE BELIEFS
muhat = vkf_bin(u,lambda,v_0,omega);
muhat = 1./(1+exp(-muhat));

% SIMULATE RESPONSES
% Define the mus of the two accumulators
mu_c1 = a + b_val.*double(u==1) + b.*(.5-muhat);
mu_c0 = a + b_val.*double(u==0) + b.*(.5-(1-muhat));
rt = nan(length(u),1); resp = nan(length(u),1);
for n = 1:length(u) % looping over the trial list
    % Defining from the first accumulator probability distribution
    P1 =    lognpdf(0.01:0.01:3,mu_c1(n),sigma);
    % Sampling from the first accumulator probability distribution
    P1 = randsample(0.01:0.01:3, 1, true,P1);
    % Defining from the second accumulator probability distribution
    P2 =    lognpdf(0.01:0.01:3,mu_c0(n),sigma);
    % Sampling from the second accumulator probability distribution
    P2 = randsample(0.01:0.01:3, 1, true,P2);
    P = [P1 P2];
    % Taking the minimum probability from the two
    % probability distribution
    [rt(n,1), resp(n,1)] = min(P);
end

rt = rt+Ter;
resp(resp==2)=0;% Setting the response of the second accumulator to 0
y = [rt resp];

%% SECTION 2: fit the model
% Customize Configuration
c = PAM_lnr_vkf_config;
% Uncomment lines below To estimate Ter
% c.Ter_fix = 0;
% c.Ter_x0 = 0;

c = PAM_lnr_vkf_fitModel(u,y,c);
% c = PAM_ddm_vkf_fitModel(u,y,'PAM_lnr_vkf_config');
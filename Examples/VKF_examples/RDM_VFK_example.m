%% Initialize path
clear
addpath(genpath(fullfile('..','..','PAM_master')))

load u.mat
%% SECTION 1: Simulate dataset

% TODO: Set parameters ------------------------------------------------
lambda = .5; % VKF: lambda
v0 = .5; % VKF: v0
omega = .1; % VKF: omega
a_a = 2; % RDM: Intercept of decision threshold "a"
b_a = 1.2; % RDM: Muhat slope effect for "a"
a_v = 2.5; % RDM: Intercept of drift rate "v"
b_val = 2.5; % RDM: Effect of validity (resp = input) on the drift
b_v = 0; % RDM: Muhat slope effect for "v"
Ter = 0; % Non decision time
%----------------------------------------------------------------------

% SIMULATE BELIEFS
muhat = vkf_bin(u,lambda,v0,omega);
muhat = 1./(1+exp(-muhat));

% SIMULATE RESPONSES
% Calculate trial-wise threshold for both the accumulators
a_c1 = a_a + b_a.*(.5-muhat);
a_c0 = a_a + b_a.*(.5-(1-muhat));

% Calculate drift rate for both accumulators
drift_c1 = a_v + b_val.*(u==1) +  b_v .* (muhat - .5);
drift_c0 = a_v + b_val.*(u==0) +  b_v .* ((1-muhat) - .5);

rt = nan(length(u),1); resp = nan(length(u),1);
for n = 1:length(u) % looping over the trial list
    probs_1 =    RDM_pdf(0.01:0.01:3,drift_c1(n),a_c1(n));
    P1 = randsample(0.01:0.01:3, 1, true,probs_1);
    probs_2 =    RDM_pdf(0.01:0.01:3,drift_c0(n),a_c0(n));
    P2 = randsample(0.01:0.01:3, 1, true,probs_2);
    P = [P1 P2];
    [rt(n,1), resp(n,1)] = min(P);
end
rt = rt+Ter;
resp(resp==2)=0;
y = [rt resp];

%% SECTION 2: fit the model
% Customize Configuration
c = PAM_rdm_vkf_config;
c.b_v_fix = 1; % To avoid fitting RDM b_v
c = PAM_rdm_vkf_fitModel(u,y,c);
% or run default configuration:
% c = PAM_rdm_vkf_fitModel(u,y,'PAM_rdm_vkf_config');
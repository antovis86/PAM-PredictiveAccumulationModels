%% Tutorial - Fitting of the model (no HGF) 

% This tutorial demonstrates how to fit the Wiener First Passage of time model (WFPT) to simulated data using MATLAB. 
% The WFPT is a popular model in cognitive neuroscience and psychology for modeling decision-making processes. 
% It describes how noisy evidence is accumulated over time until a decision boundary is reached, resulting in a response. 

% In this tutorial, we will:

% 1) Simulate data using specified DDM parameters.
% 2) Fit the model to the simulated data using maximum likelihood estimation via fminunc.
% 3) Compare the fitted parameters to the original parameters used for data generation. 

clear all
load u.mat % Load trial list
u = double(u==0);

% starting point 
a = 1.2; 

% drift rate
v = 1.9; 


% non-decision time
T = .150; 

% bias value
w = 0; 

% Initalize an empty array for response times  
rt = nan(length(u),nsubj);

% Initialize an empty array for participant responses
resp = nan(length(u),nsubj);

for n = 1:length(u) % looping over the trial list
    P1 = utl_wfpt(.001:.001:3, -v, a, 1-w);
    P2 = utl_wfpt(.001:.001:3, v, a, w);
    P = [P2(end:-1:1) P1];
    P = randsample([-3:.001:-.001 .001:.001:3], 1, true,P);
    rt(n) = abs(P);
    resp(n) = double(P>0);
end

rt = rt + T; 


initial_params = [0.5 0.5 0.5 .15];
original_parameters = [a v bw T];

% Optimization options
options = optimoptions('fminunc', 'Display', 'iter', 'MaxIterations', 1000, 'OptimalityTolerance', 1e-6);

% Run the optimization
[best_params, fval] = fminunc(@(params) ddm_negLogLikelihood(params, rt, resp), initial_params, options);

% Display the fitted parameters
disp('Recovered parameters:');
disp(best_params);

% Display the parameters used to generate the data
disp('Original parameters')
disp(original_parameters);

%% Tutorial - Fitting the model in combination with HGF and TAPAS library
% Influence the starting point "a" with priors from the  perceptual model
% This part of the script uses the Hierarchical Gaussian Filter model (HGF)
% and has to be run in combination with the TAPAS library, specifically 
% with the HGF Toolbox, open source code available as part of the TAPAS
% software collection: Frässle, S., et al. (2021). TAPAS: An Open-Source Software Package 
% for Translational Neuromodeling and Computational Psychiatry. Frontiers in Psychiatry, 12:680811. 
% https://www.translationalneuromodeling.org/tapas


%% Simulate the perceptual model
clear
tapas_init('HGF')

load u.mat % Load trial list
u = double(u==0);
om2 = -4;
% Simulate HGF
esim = tapas_simModel(u,...
    'tapas_ehgf_binary',...
    [NaN 0 1 NaN 1 1 NaN 0 0 1 1 NaN om2 -Inf]);

% extract the priors from the HGF
muhat = esim.traj.muhat(:,1);
%% Simulation

% starting point 
a = 1.2; 

% drift rate
v = 1.9; 

% scaling parameter for muhat influence on "a"
ba = .7; 

% non-decision time
T = .150; 

% scaling parameter for muhat influence on "w"
bw = 0; 

% scaling parameter for muhat influence on "v"
bv = 0; 

% just one subject
nsubj = 1; 

% cutting out the last level of hgf
prc_model = tapas_ehgf_binary_config;
prc_model.ommu(end) = -Inf;
prc_model.omsa(end) = 0;

% Initialize the priors of the observed model 
obs_model = ddm_hgf_config;
obs_model.bvmu = 0;
obs_model.bvsa = 0;
obs_model.bwmu = 0;
obs_model.bwsa = 0;
obs_model = tapas_align_priors(obs_model);

% Compute the parameters 
w = .5 + bw.*(muhat - .5);
a = a + ba.*(abs(.5-muhat)).*a;
v = u.*(v + bv.*(muhat - .5).*2.*v) - (1-u).*(v + bv.*((1-muhat)- .5).*2.*v);

% Initialize the arrays of response times and response of the subject
rt = nan(length(u),nsubj);
resp = nan(length(u),nsubj);

for n = 1:length(u) % looping over the trial list
    P1 = utl_wfpt(.001:.001:3, -v(n), a(n), 1-w(n));
    P2 = utl_wfpt(.001:.001:3, v(n), a(n), w(n));
    P = [P2(end:-1:1) P1];
    P = randsample([-3:.001:-.001 .001:.001:3], 1, true,P);
    rt(n) = abs(P);
    resp(n) = double(P>0);
end

% add non-decision time to the response time
rt = rt+T; 

% concatenate response times and responses
y = [rt resp];

% fit the model on the recovered data 
m = tapas_fitModel(y,...
    u,...
    prc_model, ...
    obs_model,...
    'tapas_quasinewton_optim_config');

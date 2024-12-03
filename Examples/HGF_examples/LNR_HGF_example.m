%% LNR example
% This example employs the model lnr_hgf
% The first section simulates beliefs and responses of a putative agent.
% The second section shows the use of PAM
%% Initialize path
clear
tapas_init('HGF')
addpath(genpath(fullfile('..','..','PAM_master')))

load u.mat

%% SECTION 1: Simulate dataset

% TODO: Set parameters ------------------------------------------------
om2 = -4; % learning rate
a = -0.53; % mu of the lognormal distribution
b_val = -0.47; % validity effect on a
b = -0.19; % belief effect
sigma = .25; % LNR variance
Ter = 0; % non-decision time
%----------------------------------------------------------------------

% SIMULATE BELIEFS
prc_model = tapas_ehgf_binary_config;
bo = tapas_fitModel([], u, prc_model, 'tapas_bayes_optimal_binary_config');
priormus = bo.p_prc.p; priormus(end-1)=om2;
esim = tapas_simModel(u,...
    'tapas_ehgf_binary',...
    priormus);
% Extract trial-wise beliefs about u
muhat = esim.traj.muhat(:,1);

% SIMULATE RESPONSES

% Define the mus of the two accumulators
mu_c1 = a + b_val.*double(u==1) + b.*(muhat -.5);
mu_c0 = a + b_val.*double(u==0) + b.*((1-muhat) - .5);
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

%% Fit ddm_hgf
% CONFIGURE PERCEPTUAL MODEL
prc_model = tapas_ehgf_binary_config;
bo = tapas_fitModel([], u, prc_model, 'tapas_bayes_optimal_binary_config');
prc_model.ommu=bo.p_prc.om; % set priors of perceptual parameters to their Bayes-optimal values
prc_model = tapas_align_priors(prc_model);

% CONFIGURE DECISION MODEL
obs_model = lnr_hgf_config;
% Uncomment the lines below to fit Ter:
% obs_model.Termu = 0;
% obs_model.Tersa = 4;


% Fit the PAM model
m = tapas_fitModel(y,...
    u,...
    prc_model, ...
    obs_model,...
    'tapas_quasinewton_optim_config');

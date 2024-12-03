%% RDM example
% This example employs the model rdm_hgf
% The first section simulates beliefs and responses of a putative agent.
% The second section shows the use of PAM
%% Initialize path
clear
tapas_init('HGF')
addpath(genpath(fullfile('..','..','PAM_master')))

load u.mat
%% SECTION 1: Simulate dataset

% TODO: Set parameters ------------------------------------------------
om2 = -4; % HGF learning rate
a_a = 2; % RDM: Intercept of decision threshold "a"
b_a = -1.2; % RDM: Muhat slope effect for "a"
a_v = 2.5; % RDM: Intercept of drift rate "v"
b_val = 2.5; % RDM: Effect of validity (resp = input) on the drift
b_v = 0; % RDM: Muhat slope effect for "v"
Ter = 0; % Non decision time
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
% Calculate trial-wise threshold for both the accumulators
a_c1 = a_a + b_a.*(muhat - .5);
a_c0 = a_a + b_a.*((1-muhat) - .5);

% Calculate drift rate for both accumulators
drift_c1 = a_v + b_val.*(u==1) +  b_v .* (muhat - .5);
drift_c0 = a_v + b_val.*(u==0) +  b_v .* ((1-muhat) - .5);

rt=nan(length(u),1);resp=nan(length(u),1);
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


%% Fit ddm_hgf_a_surprise_config
% CONFIGURE PERCEPTUAL MODEL
prc_model = tapas_ehgf_binary_config;
bo = tapas_fitModel([], u, prc_model, 'tapas_bayes_optimal_binary_config');
prc_model.ommu=bo.p_prc.om; % set priors of perceptual parameters to their Bayes-optimal values
prc_model = tapas_align_priors(prc_model);

% CONFIGURE DECISION MODEL
obs_model = RDM_hgf_config;

% Configure a reduced model (Comment the lines below to fit a full model):
obs_model.b_vsa = 0; 
% obs_model.b_asa = 0;
obs_model = tapas_align_priors(obs_model);

% Uncomment the lines below to fit Ter:
% obs_model.Termu = 0;
% obs_model.Tersa = 4;

% Fit the PAM model
m = tapas_fitModel(y,...
    u,...
    prc_model, ...
    obs_model,...
    'tapas_quasinewton_optim_config');



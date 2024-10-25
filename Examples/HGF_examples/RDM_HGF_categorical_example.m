%% Example of RDM for multiple-choice task
%% Initialize path
clear
tapas_init('HGF')
addpath(genpath(fullfile('..','..','PAM_master')))

load example_categorical_input
%% Simulate dataset

% Set parameters
om2 = -4; % HGF learning rate
a_a = 2; % RDM: Intercept of decision threshold "a"
b_a = 1.2; % RDM: Muhat slope effect for "a"
a_v = 2.5; % RDM: Intercept of drift rate "v"
b_val = 2.5; % RDM: Effect of validity (resp = input) on the drift
b_v = 0; % RDM: Muhat slope effect for "v"
Ter = 0; % Non decision time
n_choices = 3;

% SIMULATE BELIEFS
prc_model = tapas_hgf_categorical_config;
prc_model.n_outcomes = n_choices; % adjust to the number of inputs/choices
prc_model.mu2_0mu = repmat(tapas_logit(1/prc_model.n_outcomes,1),1,prc_model.n_outcomes);
prc_model.mu2_0sa = zeros(1,prc_model.n_outcomes);
prc_model.logsa2_0mu = repmat(log(1),1,prc_model.n_outcomes);
prc_model.logsa2_0sa = zeros(1,prc_model.n_outcomes);
prc_model = tapas_align_priors(prc_model);
bo = tapas_fitModel([], u_cat, prc_model, 'tapas_bayes_optimal_categorical_config');
priormus = bo.p_prc.p; priormus(end-1)=om2;
esim = tapas_simModel(u_cat,...
    'tapas_hgf_categorical',...
    priormus);
% Extract trial-wise beliefs about u_cat
muhat = squeeze(esim.traj.muhat(:,1,:)); muhat=muhat./sum(muhat,2);



% SIMULATE RESPONSES
% Define trial-wise thresholds for the accumulators
a_c1 = a_a + b_a.*(1/n_choices-muhat(:,1));
a_c2 = a_a + b_a.*(1/n_choices-muhat(:,2));
a_c3 = a_a + b_a.*(1/n_choices-muhat(:,3));
% Define trial-wise drift rates for the accumulators
drift_c1 = a_v + b_val.* double(u_cat == 1) +  b_v.*(muhat(:,1) - 1/n_choices);
drift_c2 =a_v + b_val.* double(u_cat == 2) +  b_v.*(muhat(:,2) - 1/n_choices);
drift_c3 = a_v + b_val.* double(u_cat == 3) +  b_v.*(muhat(:,3) - 1/n_choices);

for n = 1:length(u_cat) % looping over the trial list
    probs_1 =    RDM_pdf(0.01:0.01:3,drift_c1(n),a_c1(n));
    P1 = randsample(0.01:0.01:3, 1, true,probs_1);
    probs_2 =    RDM_pdf(0.01:0.01:3,drift_c2(n),a_c2(n));
    P2 = randsample(0.01:0.01:3, 1, true,probs_2);
    probs_3 =    RDM_pdf(0.01:0.01:3,drift_c3(n),a_c3(n));
    P3 = randsample(0.01:0.01:3, 1, true,probs_3);
    P = [P1 P2 P3];
    % Taking the minimum probability from the two
    % probability distribution
    [rt(n,1), resp(n,1)] = min(P);
end
rt = rt+Ter;
y = [rt(:,1) resp(:,1)];

%% Fit RDM_hgf_categorical

% Config model for the number of choices
prc_model = tapas_hgf_categorical_config;
prc_model.n_outcomes = n_choices; % adjust to the number of inputs/choices
prc_model.mu2_0mu = repmat(tapas_logit(1/prc_model.n_outcomes,1),1,prc_model.n_outcomes);
prc_model.mu2_0sa = zeros(1,prc_model.n_outcomes);
prc_model.logsa2_0mu = repmat(log(1),1,prc_model.n_outcomes);
prc_model.logsa2_0sa = zeros(1,prc_model.n_outcomes);
% determine the Bayes optimal perceptual parameters to use as priors
bo = tapas_fitModel([], u_cat, prc_model, 'tapas_bayes_optimal_categorical_config');
prc_model.ommu=bo.p_prc.om;
prc_model = tapas_align_priors(prc_model);

% Configure decision model
obs_model = RDM_hgf_categorical_config;
    
% to avoid to estimate the influence of beliefs on the drift rate set bvmu
% and bvsa to 0 (comment the lines to fit a full model in which influence
% of beliefs are estimated for both threshold and drift)
obs_model.b_vmu = 0; 
obs_model.b_vsa = 0;

% to avoid to estimate the influence of beliefs on the threshold set bamu
% and basa to 0: 
% obs_model.b_amu = 0; 
% obs_model.b_asa = 0;
obs_model = tapas_align_priors(obs_model);




% fit the model
m = tapas_fitModel(y,...
    u_cat,...
    prc_model, ...
    obs_model,...
    'tapas_quasinewton_optim_config');

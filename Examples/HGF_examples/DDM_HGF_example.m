%% DDM example
% This example employs the model ddm_hgf
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
a_a = 1.2; % boundary separation
a_v = 2; % drift rate
b_a = 0; % influence of beliefs on a
b_w = .7; % influence of beliefs on w
b_v = 0; % no influence of beliefs on v
Ter = .150; % non-decision time
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
precision = tapas_sgm(1./(muhat.*(1-muhat))-4,1)-.5;

% SIMULATE RESPONSES
a = a_a + b_a.*(precision);
w = .5 + b_w.*(muhat - .5);
v = u.*(a_v + b_v.*(muhat - .5)) - (1-u).*(a_v + b_v.*((1-muhat)- .5));
for n = 1:length(u) % looping over the trial list
    P1 = utl_wfpt(.001:.001:3, -v(n), a(n), 1-w(n));
    P2 = utl_wfpt(.001:.001:3, v(n), a(n), w(n));
    P = [P2(end:-1:1) P1];
    P = randsample([-3:.001:-.001 .001:.001:3], 1, true,P);
    rt(n,1) = abs(P);
    resp(n,1) = double(P>0);
end
rt = rt+Ter;
y = [rt(:,1) resp(:,1)];


%% Fit ddm_hgf
% CONFIGURE PERCEPTUAL MODEL
prc_model = tapas_ehgf_binary_config;
bo = tapas_fitModel([], u, prc_model, 'tapas_bayes_optimal_binary_config');
prc_model.ommu=bo.p_prc.om; % set priors of perceptual parameters to their Bayes-optimal values
prc_model = tapas_align_priors(prc_model);

% CONFIGURE DECISION MODEL
obs_model = ddm_hgf_config;
% Configure a reduced model (Comment the lines below to fit a full model):
obs_model.b_asa = 0;
obs_model.b_vsa = 0;
obs_model = tapas_align_priors(obs_model);

% Fit the PAM model
m = tapas_fitModel(y,...
    u,...
    prc_model, ...
    obs_model,...
    'tapas_quasinewton_optim_config');



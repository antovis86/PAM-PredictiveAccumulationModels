%% Example of DDM in which the boundary separation is modulated by stimulus surprise
% This example employs the model ddm_hgf_surprise_a, which assumes that the
% boundary separation "a" is modulated by the surprise elicited by the
% observed stimulus. Modulations of "w" and "v" are as in the ddm_hgf
%% Initialize path
clear
tapas_init('HGF')
addpath(genpath(fullfile('..','..','PAM_master')))

load u.mat
%% Simulate dataset

% Set parameters
om2 = -4;
a_a = 1.2;
a_v = 2;
b_a = .7;
b_w = 0; % no influence of beliefs
b_v = 0; % no influence of beliefs
Ter = .150;

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
surprise = u.*(1./(1-log2(muhat))) + (1-u).*(1./(1-log2(1-muhat)));
surprise = -(surprise-.5);
a = a_a + b_a.*(surprise);
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


%% Fit ddm_hgf_a_surprise_config
% Config perceptual model
prc_model = tapas_ehgf_binary_config;
bo = tapas_fitModel([], u, prc_model, 'tapas_bayes_optimal_binary_config');
prc_model.ommu=bo.p_prc.om;
prc_model = tapas_align_priors(prc_model);

% Configure decision mdoel
obs_model = ddm_hgf_a_surprise_config;
% Comment the lines below to fit a full model:
obs_model.b_vmu = 0;
obs_model.b_vsa = 0;
obs_model.b_wmu = 0;
obs_model.b_wsa = 0;
obs_model = tapas_align_priors(obs_model);

% Fit the model
m1 = tapas_fitModel(y,...
    u,...
    prc_model, ...
    obs_model,...
    'tapas_quasinewton_optim_config');



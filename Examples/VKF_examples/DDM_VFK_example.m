%% Initialize path
clear
addpath(genpath(fullfile('..','..','PAM_master')))

load u.mat
%% SECTION 1: Simulate dataset

% TODO: Set parameters ------------------------------------------------
lambda = .5; % learning rate
v_0 = .5;
omega = .1;
a_a = 1.2; % boundary separation
a_v = 2; % drift rate
b_a = 0; % influence of beliefs on a
b_w = .7; % influence of beliefs on w
b_v = 0; % no influence of beliefs on v
Ter = .150; % non-decision time
%----------------------------------------------------------------------

% SIMULATE BELIEFS
muhat = vkf_bin(u,lambda,v_0,omega);
muhat = 1./(1+exp(-muhat));

% SIMULATE RESPONSES
a = a_a + b_a.*(abs(.5-muhat));
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

%% SECTION 2: fit the model
% Customize Configuration
c = PAM_ddm_vkf_config;
c.b_a_fix = 1; % To avoid fitting DDM ba
c.b_v_fix = 1; % To avoid fitting DDM bv
c = PAM_ddm_vkf_fitModel(u,y,c);
% c = PAM_ddm_vkf_fitModel(u,y,'PAM_ddm_vkf_config');
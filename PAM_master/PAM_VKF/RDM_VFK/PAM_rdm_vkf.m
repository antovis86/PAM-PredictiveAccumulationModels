function logp = PAM_rdm_vkf(x,u,y,c)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function implement the PAM framework using the following models:
% PERCEPTUAL MODEL: Volatile Kalman Filter (https://doi.org/10.1371/journal.pcbi.1007963)
% DECISION MODEL: Racing Diffusion Model (https://doi.org/10.3758/s13423-020-01719-6)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% --------------------------------------------------------------------------------------------------
% Copyright (C) 2024 Antonino Visalli
%
% This file is part of the PAM toolbox, which is released under the terms of the GNU General Public
% Licence (GPL), version 3. You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version). For further details, see the file
% COPYING or <http://www.gnu.org/licenses/>.

ptrans = c.initpt;
ptrans(~c.fix)=x;
ptrans = PAM_rdm_vkf_transp(ptrans,y);

% VKF parameters
lambda = ptrans(1);
v0 = ptrans(2);
omega = ptrans(3);
muhat = vkf_bin(u,lambda,v0,omega);
muhat = 1./(1+exp(-muhat)); % Belief trajectory

% RDM parameters
a_a = ptrans(4);
b_a = ptrans(5);
a_v = ptrans(6);
b_val = ptrans(7);
b_v = ptrans(8);


% irregualar trials
irr = isnan(y(:,1));

% Weed irregular trials out from inferred states, responses, and inputs
muhat(irr) = [];

rt = y(:,1);
rt(irr) = [];

resp = y(:,2);
resp(irr)  = [];

% Fitting the non-decision time with the minimum value of estimated
% non-decision time
Ter = ptrans(9);
rt = max(eps,rt-Ter);

% extract the trial list and remove the irregular trials
u(irr) = [];


% Calculate trial-wise threshold for both the accumulators
a_c1 = a_a + b_a.*(muhat - .5);
a_c0 = a_a + b_a.*((1-muhat) - .5);

% Calculate drift rate for both accumulators
drift_c1 = a_v + b_val.*(u==1) +  b_v .* (muhat - .5);
drift_c0 = a_v + b_val.*(u==0) +  b_v .* ((1-muhat) - .5);

% Calculate predicted log-likelihood
logp_reg = NaN(length(u),1);
for ntrial = 1:length(u)
    if rt(ntrial)>0
        % Extract the drifts for the specific trial
        drift_pdf = (resp(ntrial) == 1)*drift_c1(ntrial) + (resp(ntrial) == 0) * drift_c0(ntrial);
        drift_cdf = (resp(ntrial) == 1)*drift_c0(ntrial) + (resp(ntrial) == 0) * drift_c1(ntrial);

        % Extract the threshold for the specific trial
        a_pdf = (resp(ntrial) == 1)*a_c1(ntrial) + (resp(ntrial) == 0) * a_c0(ntrial);
        a_cdf = (resp(ntrial) == 1)*a_c0(ntrial) + (resp(ntrial) == 0) * a_c1(ntrial);        
        % Compute the defective distribution
        P = utl_inverse_gaussian_defective(rt(ntrial),drift_pdf,drift_cdf,a_pdf,a_cdf);

        if P>0
            logp_reg(ntrial) = log(P+eps);
        else
            logp_reg(ntrial) = NaN;
        end
    end
end


logp = -sum(logp_reg);
return;

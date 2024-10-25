function logp = PAM_lnr_vkf(x,u,y,c)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function implement the PAM framework using the following models:
% PERCEPTUAL MODEL: Volatile Kalman Filter (https://doi.org/10.1371/journal.pcbi.1007963)
% DECISION MODEL: Log-normal race model (https://doi.org/10.3389/fpsyg.2012.00292)
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
ptrans = PAM_lnr_vkf_transp(ptrans,y);

% VKF parameters
lambda = ptrans(1);
v0 = ptrans(2);
omega = ptrans(3);
muhat = vkf_bin(u,lambda,v0,omega);
muhat = 1./(1+exp(-muhat)); % Belief trajectory

% LNR parameters
a = ptrans(4);
b_val = ptrans(5);
b = ptrans(6);
sigma = ptrans(7);



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
Ter = ptrans(8);
rt = max(eps,rt-Ter);

% extract the trial list and remove the irregular trials
u(irr) = [];


% Calculate trial-wise drift rates for the two accumulators
mu_c1 = a + b_val.* double(u == 1) + b .*(.5 -muhat);
mu_c0 = a + b_val.* double(u == 0) + b .*(.5 - (1 -muhat));



% Calculate predicted log-likelihood
logp_reg = NaN(length(u),1);
for ntrial = 1:length(u)
    if rt(ntrial)>0
       mu_pdf = (resp(ntrial) == 1)*mu_c1(ntrial) + (resp(ntrial) == 0) * mu_c0(ntrial);
        mu_cdf = (resp(ntrial) == 1)*mu_c0(ntrial) + (resp(ntrial) == 0) * mu_c1(ntrial);
        P = utl_lnr_pdf(rt(ntrial),mu_pdf,mu_cdf,sigma);
        
        if P>0
            logp_reg(ntrial) = log(P+eps);
        else
            logp_reg(ntrial) = NaN;
        end
    end
end


logp = -sum(logp_reg);
return;

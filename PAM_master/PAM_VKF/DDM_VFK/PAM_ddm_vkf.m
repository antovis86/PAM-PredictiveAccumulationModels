function logp = PAM_ddm_vkf(x,u,y,c)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function implement the PAM framework using the following models:
% PERCEPTUAL MODEL: Volatile Kalman Filter (https://doi.org/10.1371/journal.pcbi.1007963)
% DECISION MODEL: Wiener diffusion model (https://doi.org/10.1016/j.jmp.2009.02.003)
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
ptrans = PAM_ddm_vkf_transp(ptrans,y);

% VKF parameters
lambda = ptrans(1);
v0 = ptrans(2);
omega = ptrans(3);
muhat = vkf_bin(u,lambda,v0,omega);
muhat = 1./(1+exp(-muhat)); % Belief trajectory

% DDM parameters
a_a = ptrans(4);
a_v = ptrans(5);
b_w = ptrans(6);
b_a = ptrans(7);
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


% Calculate trial-wise starting point
w = .5 + b_w.*(muhat - .5);


% Calculate trial-wise absorbing barrier
a = a_a + b_a.*(abs(.5-muhat));

% Calculate trial-wise drift
v = u.*(a_v + b_v.*(muhat - .5)) ...
    - (1-u).*(a_v + b_v.*((1-muhat) - .5));


% Calculate predicted log-likelihood
logp_reg = NaN(length(u),1);
for ntrial = 1:length(u)
    if rt(ntrial)>0
        P = utl_wfpt(rt(ntrial), -v(ntrial), a(ntrial), (1-w(ntrial)))*resp(ntrial) + ...
            utl_wfpt(rt(ntrial), v(ntrial), a(ntrial), w(ntrial))*(1-resp(ntrial));

        if P>0
            logp_reg(ntrial) = log(P+eps);
        else
            logp_reg(ntrial) = NaN;
        end
    end
end


logp = -sum(logp_reg);
return;

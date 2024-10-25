function [logp, yhat, res] = ddm_hgf_a_surprise(r, infStates, ptrans)
% Calculates the log-probability of responses
%
%
%
% --------------------------------------------------------------------------------------------------
% Copyright (C) 2012-2013 Christoph Mathys, TNU, UZH & ETHZ
%
% This file is part of the HGF toolbox, which is released under the terms of the GNU General Public
% Licence (GPL), version 3. You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version). For further details, see the file
% COPYING or <http://www.gnu.org/licenses/>.

% Transform parameters to their native space
a_a = exp(ptrans(1));
a_v = exp(ptrans(2));
b_w = 2/(1+exp(-ptrans(3)))-1;
b_a = ptrans(4);
b_v = ptrans(5);


% Initialize returned log-probabilities, predictions,
% and residuals as NaNs so that NaN is returned for all
% irregualar trials
n = size(infStates,1);
logp = NaN(n,1);
yhat = NaN(n,1); % not used
res  = NaN(n,1); % not used

% Weed irregular trials out from inferred states, responses, and inputs
mu1hat = infStates(:,1,1);
mu1hat(r.irr) = [];

rt = r.y(:,1);
rt(r.irr) = [];

resp = r.y(:,2);
resp(r.irr)  = [];

% Fitting the non-decision time with the minimum value of estimated
% non-decision time
Ter = min(rt)/(1+exp(-ptrans(6)));
rt = max(eps,rt-Ter);

% extract the trial list and remove the irregular trials
u = r.u(:,1);
u(r.irr) = [];

% Calculate trial-wise starting point
w = .5 + b_w.*(mu1hat - .5);

% Calculate trial-wise absorbing barrier
surprise = u.*(1./(1-log2(mu1hat))) + (1-u).*(1./(1-log2(1-mu1hat)));
surprise = -(surprise-.5); 
a = a_a + b_a.*(surprise);

% Calculate trial-wise drift
v = u.*(a_v + b_v.*(mu1hat - .5)) ...
    - (1-u).*(a_v + b_v.*((1-mu1hat) - .5));

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

reg = ~ismember(1:n,r.irr);
logp(reg) = logp_reg;
return;

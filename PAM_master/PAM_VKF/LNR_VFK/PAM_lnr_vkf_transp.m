function [pvec] = PAM_lnr_vkf_transp(par_vec,y)
% Copyright (C) 2024 Antonino Visalli
%
% This file is part of the PAM toolbox, which is released under the terms of the GNU General Public
% Licence (GPL), version 3. You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version). For further details, see the file
% COPYING or <http://www.gnu.org/licenses/>.
% --------------------------------------------------------------------------------------------------
% Rdm_vkf_transp: Function to perform transformation of parameters
% [Inputs]
% - par_vec: vector of the parameters 
% - y: responses

% Initialize an empty array to store the transformed values
pvec    = NaN(1,length(par_vec));

%% VKF Parameters
pvec(1) = 1/(1+exp(-par_vec(1))); % lamba
pvec(2) = exp(par_vec(2)); % v0
pvec(3) = exp(par_vec(3)); % omega

%% LNR Parameters 
pvec(4) = par_vec(4); % a
pvec(5) = par_vec(5); % b_val
pvec(6) = par_vec(6); % b
pvec(7) = exp(par_vec(7)); % sigma
pvec(8) = min(y(:,1))/(1+exp(-par_vec(8))); % Ter
return;
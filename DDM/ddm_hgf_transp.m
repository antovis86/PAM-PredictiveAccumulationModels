function [pvec, pstruct] = ddm_hgf_transp(r, ptrans) 
% --------------------------------------------------------------------------------------------------


% ddm_hgf_transp: Function to perform transformation of parameters

% [Inputs]
% - r: array of responses
% - ptrans: structure that contains the parameters


% The structure and methodologies of this file are inspired
% from the HGF Toolbox, open source code available as part of the TAPAS
% software collection: Frässle, S., et al. (2021). TAPAS: An Open-Source Software Package 
% for Translational Neuromodeling and Computational Psychiatry. Frontiers in Psychiatry, 12:680811. 
% https://www.translationalneuromodeling.org/tapas


% extracting times
y = r.y(:,1);
y(r.irr) = [];

% Initialize an empty array to store the transformed values
pvec    = NaN(1,length(ptrans));
pstruct = struct;

% Apply exponential transformation to "a"
pvec(1)      = exp(ptrans(1));
pstruct.a   = pvec(1);

% Apply exponential transformation to "v"
pvec(2)      = exp(ptrans(2));
pstruct.v   = pvec(2);


% Transform "bw" values in the range (-1,1)
pvec(3)      = 2/(1+exp(-ptrans(3)))-1;
pstruct.bw   = pvec(3);

% Transform "ba" values in the range (-1,1)
pvec(4)      = 2/(1+exp(-ptrans(4)))-1;
pstruct.ba   = pvec(4);

% Transform "bv" values in the range (-1,1)
pvec(5)      = 2/(1+exp(-ptrans(5)))-1;
pstruct.bv   = pvec(5);

% Transform "T" using a sigmoid centered at the minimum value of y 
pvec(6)      = min(y)/(1+exp(-ptrans(6)));
pstruct.Tmu = pvec(6);

return;
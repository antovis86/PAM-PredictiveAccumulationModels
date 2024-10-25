function c = PAM_lnr_vkf_fitModel(u,y,c)
% This is the main function for fitting the parameters of the "PAM_lnr_vkf"
% model
% --------------------------------------------------------------------------------------------------
% Copyright (C) 2024 Antonino Visalli
%
% This file is part of the PAM toolbox, which is released under the terms of the GNU General Public
% Licence (GPL), version 3. You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version). For further details, see the file
% COPYING or <http://www.gnu.org/licenses/>.
if ischar(c)
    c = eval(c);
end

% Replace default initpt and fix vectos with newly ones
initpt = []; % Initialize new initial points
fix = logical([]); % Initialize new fix flag values
names = fieldnames(c); % Get fieldnames.
for i = 1:length(names) % Loop over fields
    if regexp(names{i}, '_x0$')
        initpt = [initpt, c.(names{i})];
    elseif regexp(names{i}, '_fix$')
        fix = [fix, c.(names{i})];
    end
end
c.initpt = initpt;
c.fix = fix;


% FIT MODEL
x0 = c.initpt(~c.fix);
fun = @(x) PAM_lnr_vkf(x,u,y,c);
options = optimoptions('fminunc', 'Algorithm', 'quasi-newton');
x = fminunc(fun, x0, options);

c.final = c.initpt;
c.final(~c.fix) = x;
c.final_transp = PAM_lnr_vkf_transp(c.final,y);

c.final_prc.lambda = c.final_transp(1);
c.final_prc.v0 = c.final_transp(2);
c.final_prc.omega = c.final_transp(3);

c.final_obs.a = c.final_transp(4);
c.final_obs.b_val = c.final_transp(5);
c.final_obs.b = c.final_transp(6);
c.final_obs.sigma = c.final_transp(7);
c.final_obs.Ter = c.final_transp(8);


disp('Results:');
disp(' ')
disp('Parameter estimates for the perceptual model:');
disp(c.final_prc)
disp(' ')
disp('Parameter estimates for the observation model:');
disp(c.final_obs)


end


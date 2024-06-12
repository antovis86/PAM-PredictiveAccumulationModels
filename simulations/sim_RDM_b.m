%% Fitting the perceptual model
clear all
tapas_init('HGF')
% Load trial list
load u.mat
u = double(u==0);
om2 = -4;

% Simulate HGF
esim = tapas_simModel(u,...
    'tapas_ehgf_binary',...
    [NaN 0 1 NaN 1 1 NaN 0 0 1 1 NaN om2 -Inf]);
muhat = esim.traj.muhat(:,1);


%% Parameter Recovery %%

% Define search space of Inverse Gaussian values
b0B_all = [1.5 3];
b1_all = [.3 .7];
Vv_all = [4];
bv_all = 0;
Vi_all = [1 2 3];

% Define the number of subject to simulate
nsubj = 100;

% cutting out the last level of hgf
prc_model = tapas_ehgf_binary_config;
prc_model.ommu(end) = -Inf;
prc_model.omsa(end) = 0;


% Cancel out estimation of bv
obs_model = RDM_hgf_config;
obs_model.bvmu = 0;
obs_model.bvsa = 0;

% Generate empty arrays to collect results
tot = length(b0B_all)*length(b1_all)*length(Vv_all)*length(Vi_all)*nsubj;
p_obs = nan(length(b0B_all),length(b1_all),length(Vv_all),length(Vi_all),nsubj,5);
p_prc = nan(length(b0B_all),length(b1_all),length(Vv_all),length(Vi_all),nsubj,2);

% Initialize a table to collect the different models
models = table();

for b0Bx = 1:length(b0B_all) % looping over "a" values
    for b1x = 1:length(b1_all) % looping over "b1mu" values
        for Vvx = 1:length(Vv_all)
            for Vix = 1:length(Vi_all)
                for bvx = 1:length(bv_all)

                    % extracting "b0" value
                    b0B = b0B_all(b0Bx);

                    % extracting "b1" value
                    b1 = b1_all(b1x);

                    % extracting "b2" values
                    Vv = Vv_all(Vvx);

                    % extracting "Vi" values
                    Vi = Vi_all(Vix);

                    % extracting "bv" values
                    bv = bv_all(bvx);

                    % Calculate threshold
                    B_dx = b0B + b1.*(.5-muhat).*2.*b0B;
                    B_sx = b0B + b1.*(.5-(1-muhat)).*2*b0B;


                    % Calculate drift rate
                    drift_r1 = Vv .* double(u == 1) + Vi .* double(u == 0);
                    drift_r1 = drift_r1 + (bv .* (muhat - .5).*2.*drift_r1);

                    drift_r2 = Vv .* double(u == 0) + Vi .* double(u == 1);
                    drift_r2 = drift_r2 + (bv.*((1-muhat) - .5).*2.*drift_r2);

                    % Initialize empty arrays
                    rt = nan(length(u),1);
                    resp = nan(length(u),1);


                    % Simulate participants
                    for idx = 1:nsubj
                        for n = 1:length(u) % looping over the trial list

                            probs_1 =    RDM_pdf(0.001:0.001:3,drift_r1(n),B_dx(n));
                            P1 = randsample(0.001:0.001:3, 1, true,probs_1);
                            probs_2 =    RDM_pdf(0.001:0.001:3,drift_r2(n),B_sx(n));
                            P2 = randsample(0.001:0.001:3, 1, true,probs_2);
                            P = [P1 P2];
                            [rt(n,idx), resp(n,idx)] = min(P);

                        end
                    end

                    resp(resp==2)=0;


                    % Parameter recovery

                    for idx = 1:nsubj % fit the model on every simulated subject
                        % modify the priors and parameter's research space
                        y = [rt(:,idx) resp(:,idx)];

                        % cutting out the last level of hgf
                        m = tapas_fitModel(y,...
                            u,...
                            prc_model, ...
                            'RDM_hgf_config',...
                            'tapas_quasinewton_optim_config');
                        p_obs(b0Bx,b1x,Vvx,Vix,bvx,idx,:)=m.p_obs.p;
                        p_prc(b0Bx,b1x,Vvx,Vix,bvx,idx,:)=m.p_prc.om(2:3);

                        % Populate the table

                        % 1) create a temporary table
                        temp_table = table();

                        % 2) Add the parameters that have to be recovered
                        temp_table.b0B = repmat(b0B_all(b0Bx), numel(idx), 1); % Replicate 'b0B' values for each subject
                        temp_table.b1 = repmat(b1_all(b1x), numel(idx), 1); % Replicate 'b1' values for each subject
                        temp_table.Vv = repmat(Vv_all(Vvx), numel(idx), 1); % Replicate 'Vv' values for each subject
                        temp_table.Vi = repmat(Vi_all(Vix), numel(idx), 1); % Replicate 'Vi' values for each subject
                        temp_table.bv = repmat(bv_all(bvx), numel(idx), 1); % Replicate 'bv' values for each subject

                        % Add the subject index
                        temp_table.subject = ones(numel(idx), 1) * idx; % Assign the subject index

                        % Add the model
                        temp_table.model = repmat({m}, numel(idx), 1); % Replicate 'm' object for each subject

                        % Concatenate the temporary table with the models table
                        models = [models; temp_table];
                    end
                end
            end
        end
    end
end



% Saving the outputs
save(fullfile('results','SIM_threshold.mat'),'p_obs','p_prc','b0B_all','b1_all',"Vv_all","Vi_all","bv_all")
save(fullfile('results','table_threshold.mat'),'models')
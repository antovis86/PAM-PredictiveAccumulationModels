%% Initiate paths
addpath(genpath('../../PAM/RDM'))
addpath('..')
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
b0B_all = [2 3];
b1B_all = 0;
Vv_all = [5 4];
bv_all =  [.3 .7];
Vi_all = [.5 .71];
T = .150;

% Define the number of subject to simulate
nsubj = 100;

% cutting out the last level of hgf
prc_model = tapas_ehgf_binary_config;
prc_model.ommu(end) = -Inf;
prc_model.omsa(end) = 0;
prc_model.logkamu(end) = log(0);
bopars = tapas_fitModel([],...
    u,...
    prc_model, ...
    'tapas_bayes_optimal_binary_config',...
    'tapas_quasinewton_optim_config');
prc_model.ommu = bopars.p_prc.om;

% Generate empty arrays to collect results
tot = length(b0B_all)*length(b1B_all)*length(Vv_all)*length(Vi_all)*nsubj;
p_obs = nan(length(b0B_all),length(b1B_all),length(Vv_all),length(Vi_all),nsubj,6);
p_prc = nan(length(b0B_all),length(b1B_all),length(Vv_all),length(Vi_all),nsubj,2);

% Cancel estimation of b1
obs_model = RDM_hgf_config;
obs_model.b1Bmu = 0;
obs_model.b1Bsa = 0;

obs_model = tapas_align_priors(obs_model);

% Initialize a table to collect the different models
models = table();

for sx = 1:length(b0B_all) % looping over "a" values
    for sx2=1:length(Vi_all)
        for sx3=1:length(bv_all)

            % extracting "b0" value
            b0B = b0B_all(sx);

            % extracting "b1" value
            b1B = b1B_all;

            % extracting "b2" values
            Vv = Vv_all(sx);

            % extracting "Vi" percentage
            Vi = Vi_all(sx2);

            % Calculating Vi
            Vi = (Vi*Vv);

            % extracting "bv" values
            bv = bv_all(sx3);

            % Calculate threshold
            B_dx = b0B + b1B.*(.5-muhat).*2.*b0B;
            B_sx = b0B + b1B.*(.5-(1-muhat)).*2*b0B;


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

                    probs_1 =    RDM_pdf(0.01:0.01:3,drift_r1(n),B_dx(n));
                    P1 = randsample(0.01:0.01:3, 1, true,probs_1);
                    probs_2 =    RDM_pdf(0.01:0.01:3,drift_r2(n),B_sx(n));
                    P2 = randsample(0.01:0.01:3, 1, true,probs_2);
                    P = [P1 P2];
                    [rt(n,idx), resp(n,idx)] = min(P);

                end
            end

            resp(resp==2)=0;
            rt = rt +T;


            % Parameter recovery

            parfor idx = 1:nsubj % fit the model on every simulated subject
                % modify the priors and parameter's research space
                y = [rt(:,idx) resp(:,idx)];

                % cutting out the last level of hgf
                m1 = tapas_fitModel(y,...
                    u,...
                    prc_model, ...
                    obs_model,...
                    'tapas_quasinewton_optim_config');
                m2 = tapas_fitModel(y,...
                    u,...
                    prc_model, ...
                    RDM_hgf_config,...
                    'tapas_quasinewton_optim_config');
                p_obs(sx,sx,sx,sx2,idx,:)=m1.p_obs.p;
                p_prc(sx,sx,sx,sx2,idx,:)=m1.p_prc.om(2:3);

                % Populate the table

                % 1) create a temporary table
                temp_table = table();

                % 2) Add the parameters that have to be recovered
                temp_table.b0B = repmat(b0B, numel(idx), 1); % Replicate 'b0B' values for each subject
                temp_table.b1 = repmat(b1B, numel(idx), 1); % Replicate 'b1' values for each subject
                temp_table.Vv = repmat(Vv, numel(idx), 1); % Replicate 'Vv' values for each subject
                temp_table.Vi = repmat(Vi, numel(idx), 1); % Replicate 'Vi' values for each subject
                temp_table.bv = repmat(bv, numel(idx), 1); % Replicate 'bv' values for each subject
                temp_table.T = repmat(T, numel(idx), 1); % Replicate 'T' values for each subject


                % Add the subject index
                temp_table.subject = ones(numel(idx), 1) * idx; % Assign the subject index

                % Add the model
                temp_table.model_1 = repmat({m1}, numel(idx), 1); % Replicate 'm' object for each subject
                temp_table.model_2 = repmat({m2}, numel(idx), 1); % Replicate 'm' object for each subject


                % Concatenate the temporary table with the models table
                models = [models; temp_table];
            end
        end
    end
end


% Saving the outputs
save D:\PAM\PAM_02_08_2024\results\rdm_m1_m2\rdm_HGF_V_sim_results_seconds_non_decision_time.mat models

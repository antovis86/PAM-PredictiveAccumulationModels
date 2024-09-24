%% Initiate paths
% set path to HGF
% addpath('/home/ccnl/neuroscience/codici/tapas-master')

addpath(genpath('../../PAM/LNR'))
%% Simulate the perceptual model
clear
tapas_init('HGF')

load ../u.mat % Load trial list
u = double(u==0);
om2 = -4;
% Simulate HGF
esim = tapas_simModel(u,...
    'tapas_ehgf_binary',...
    [NaN 0 1 NaN 1 1 NaN 0 0 1 1 NaN om2 -Inf]);
muhat = esim.traj.muhat(:,1);

%% Simulation
% Parameter space
bv_all = [-1 -1 -.3 -.3];
bi_all = [-.53 -.75 .19 -0.1];
b1_all = [0.3 0.7];
sigma = .25;
T = .150;

% Define the number of subjects
nsubj = 100;

%
prc_model = tapas_ehgf_binary_config;
prc_model.ommu(end) = -Inf;
prc_model.omsa(end) = 0;
prc_model.logkamu(end) = log(0);
bopars = tapas_fitModel([],...
    u,...
    prc_model, ...
    'tapas_bayes_optimal_binary_config',...
    'tapas_quasinewton_optim_config');
% cutting out the last level of ehgf
prc_model.ommu = bopars.p_prc.om;



% Initialize empty arrays to store results
p_obs = nan(length(bv_all),length(bi_all),length(b1_all),nsubj,4);
p_prc = nan(length(bv_all),length(bi_all),length(b1_all),nsubj,1);

% Initialize an empty table to store models and results
models = table();

% looping over scenarios
for sx = 1:length(bv_all)
    for b1x = 1:length(b1_all)


        % extract the right parameter
        bv = bv_all(sx);
        bi = bi_all(sx);

        b1 = b1_all(b1x);

        % Define the mus of the two accumulators
        mu_r1 = bv.*double(u==1) + bi.*double(u==0)+(b1.*(.5-muhat));
        mu_r2 = bv.*double(u==0) + bi.*double(u==1)+(b1.*(.5-(1-muhat)));

        % Initialize an empty array to store the rts
        rt = nan(length(u),nsubj); resp = nan(length(u),nsubj);

        for idx = 1:nsubj  % looping over the number of subjects
            for n = 1:length(u) % looping over the trial list

                % Defining from the first accumulator probability distribution
                P1 =    lognpdf(0.01:0.01:3,mu_r1(n),sigma);
                % Sampling from the first accumulator probability distribution
                P1 = randsample(0.01:0.01:3, 1, true,P1);
                % Defining from the second accumulator probability distribution
                P2 =    lognpdf(0.01:0.01:3,mu_r2(n),sigma);
                % Sampling from the second accumulator probability distribution
                P2 = randsample(0.01:0.01:3, 1, true,P2);
                P = [P1 P2];
                % Taking the minimum probability from the two
                % probability distribution
                [rt(n,idx), resp(n,idx)] = min(P);
            end
        end
        % Setting the response of the second accumulator to 0
        resp(resp==2)=0;
        rt = rt +T;
        % looping over the number of simulated subjects
        
        parfor idx = 1:nsubj
            y = [rt(:,idx) resp(:,idx)];
            obs_model = lnr_hgf_config;
            %obs_model.bvmu=log((nanmean(y(:,1))));
            %obs_model.bimu=log((nanmean(y(:,1))));
            %obs_model.sigmamu = log(log(std(y(:,1))));
            obs_model = tapas_align_priors(obs_model);
            % fitting the model
            m = tapas_fitModel(y,...
                u,...
                prc_model, ...
                obs_model,...
                'tapas_quasinewton_optim_config');

            p_obs(sx,sx,b1x,idx,:)=m.p_obs.p(1:end-1);
            p_prc(sx,sx,b1x,idx,:)=m.p_prc.om(2);
            % populate the table
            % 1) create a temporary table
            temp_table = table();
            % 2) Add the parameters that have to be recovered
            temp_table.bv = repmat(bv, numel(idx), 1); % Replicate 'bv' values for each subject
            temp_table.bi = repmat(bi, numel(idx), 1); % Replicate 'bi' values for each subject
            temp_table.b1 = repmat(b1, numel(idx), 1); % Replicate 'b1' values for each subject
            temp_table.sigma = repmat(sigma, numel(idx), 1); % Replicate 'sigma' values for each subject
            temp_table.T = repmat(T, numel(idx), 1); % Replicate 'T' values for each subject
            temp_table.om2 = repmat(om2, numel(idx), 1); % Replicate 'f' values for each subject
            % Add the subject index
            temp_table.subject = ones(numel(idx), 1) * idx; % Assign the subject index
            % Add the model
            temp_table.model = repmat({m}, numel(idx), 1); % Replicate 'm' object for each subject
            % Concatenate the temporary table with the models table
            models = [models; temp_table];

        end
    end
end

% save result of the simulation
save D:\PAM\PAM_02_08_2024\results\lnr_HGF_sim_results_seconds_non_decision_time.mat models

%% Simulate the perceptual model
clear
tapas_init('HGF')

load u.mat % Load trial list
u = double(u==0);
om2 = -4;
% Simulate HGF
esim = tapas_simModel(u,...
    'tapas_ehgf_binary',...
    [NaN 0 1 NaN 1 1 NaN 0 0 1 1 NaN om2 -Inf]);
muhat = esim.traj.muhat(:,1);

%% Simulation
% Parameter space
bv_all = [5.5 6.5];
bi_all = [.05 .1];
b1_all = [.3 .7];
sigma = .25;

% Define the number of subjects
nsubj = 100;


% cutting out the last level of ehgf
prc_model = tapas_ehgf_binary_config;
prc_model.ommu(end) = -Inf;
prc_model.omsa(end) = 0;

% Initialize empty arrays to store results 
p_obs = nan(length(bv_all),length(bi_all),length(b1_all),nsubj,4);
p_prc = nan(length(bv_all),length(bi_all),length(b1_all),nsubj,1);

% Initialize an empty table to store models and results
models = table();

% looping over bv parameter
for bvx = 1:length(bv_all)
    % looping over bv parameter
    for bix = 1:length(bi_all)
        % looping over b1 parameter
        for b1x = 1:length(b1_all)
            
            % extract the right parameter 
            bv = bv_all(bvx);
            bi = bi_all(bix);
            bi = bv+bi*bv;
            b1 = b1_all(b1x);
            
            % Define the mus of the two accumulators 
            mu_r1 = bv.*double(u==1) + bi.*double(u==0)+(b1.*(.5-muhat));
            mu_r2 = bv.*double(u==0) + bi.*double(u==1)+(b1.*(.5-(1-muhat)));
            
            % Initialize an empty array to store the rts
            rt = nan(length(u),nsubj); resp = nan(length(u),nsubj);
            
            for idx = 1:nsubj  % looping over the number of subjects 
                for n = 1:length(u) % looping over the trial list
                    
                    % Defining from the first accumulator probability distribution
                    P1 =    lognpdf(1:3000,mu_r1(n),sigma);
                    % Sampling from the first accumulator probability distribution 
                    P1 = randsample(1:3000, 1, true,P1);
                    % Defining from the second accumulator probability distribution
                    P2 =    lognpdf(1:3000,mu_r2(n),sigma);
                    % Sampling from the second accumulator probability distribution 
                    P2 = randsample(1:3000, 1, true,P2);
                    P = [P1 P2];
                    % Taking the minimum probability from the two
                    % probability distribution
                    [rt(n,idx), resp(n,idx)] = min(P);
                end
            end
            % Setting the response of the second accumulator to 0
            resp(resp==2)=0;
            % looping over the number of simulated subjects
            parfor idx = 1:nsubj
                y = [rt(:,idx) resp(:,idx)];
                
                % fitting the model 
                m = tapas_fitModel(y,...
                    u,...
                    prc_model, ...
                    'lnr_hgf_config',...
                    'tapas_quasinewton_optim_config');

                p_obs(bvx,bix,b1x,idx,:)=m.p_obs.p(1:end-1);
                p_prc(bvx,bix,b1x,idx,:)=m.p_prc.om(2);
                % populate the table
                % 1) create a temporary table
                temp_table = table();
                % 2) Add the parameters that have to be recovered
                temp_table.bv = repmat(bv, numel(idx), 1); % Replicate 'bv' values for each subject
                temp_table.bi = repmat(bi, numel(idx), 1); % Replicate 'bi' values for each subject
                temp_table.b1 = repmat(b1, numel(idx), 1); % Replicate 'b1' values for each subject
                temp_table.sigma = repmat(sigma, numel(idx), 1); % Replicate 'sigma' values for each subject
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
end

save lnrHGF_sim_results.mat models

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
a_all = [1.2 1.2 2 2];
v_all = [1.9 1 1.1 .6];
bv_all = [.3 .7];
T = .150;
bw = 0;
ba = 0;
nsubj = 100;


% cutting out the last level of hgf
prc_model = tapas_ehgf_binary_config;
prc_model.ommu(end) = -Inf;
prc_model.omsa(end) = 0;

obs_model = ddm_hgf_config;
obs_model.bamu = 0;
obs_model.basa = 0;
obs_model.bwmu = 0;
obs_model.bwsa = 0;
obs_model = tapas_align_priors(obs_model);

p_obs = nan(length(a_all),length(bv_all),nsubj,4);
p_prc = nan(length(a_all),length(bv_all),nsubj,4);
models = table();

for ax = 1:length(a_all)
    for bx = 1:length(bv_all)
        a = a_all(ax);
        v = v_all(ax);
        bv = bv_all(bx);

        w = .5 + bw.*(muhat - .5);
        a = a + ba.*(abs(.5-muhat)).*a;
        v = u.*(v + bv.*(muhat - .5).*2.*v) - (1-u).*(v + bv.*((1-muhat)- .5).*2.*v);


        rt = nan(length(u),nsubj); resp = nan(length(u),nsubj);

        for idx = 1:nsubj
            for n = 1:length(u) % looping over the trial list
                P1 = utl_wfpt(.001:.001:3, -v(n), a(n), 1-w(n));
                P2 = utl_wfpt(.001:.001:3, v(n), a(n), w(n));
                P = [P2(end:-1:1) P1];
                P = randsample([-3:.001:-.001 .001:.001:3], 1, true,P);
                rt(n,idx) = abs(P);
                resp(n,idx) = double(P>0);
            end
        end
        rt = rt+T;
        for idx = 1:nsubj
            y = [rt(:,idx) resp(:,idx)];
            m = tapas_fitModel(y,...
                u,...
                prc_model, ...
                obs_model,...
                'tapas_quasinewton_optim_config');

            p_obs(ax,bx,idx,:)=m.p_obs.p([1:2 5 end]);
            p_prc(ax,bx,idx,:)=m.p_prc.om(2);
            % populate the table
            % 1) create a temporary table
            temp_table = table();
            % 2) Add the parameters that have to be recovered
            temp_table.a = repmat(a_all(ax), numel(idx), 1); % Replicate 'bv' values for each subject
            temp_table.v = repmat(v_all(ax), numel(idx), 1); % Replicate 'bi' values for each subject
            temp_table.bv = repmat(bv_all(bx), numel(idx), 1); % Replicate 'b1' values for each subject
            temp_table.T = repmat(T, numel(idx), 1); % Replicate 'sigma' values for each subject
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


% save ddmHGF_a_sim_results.mat models

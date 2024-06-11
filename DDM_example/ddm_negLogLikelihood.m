function negLogLikelihood = ddm_negLogLikelihood(params, rt, resp)
    a = params(1);
    v = params(2);
    w = params(3);
    T  = params(4); 
    
    n = length(rt);
    logLikelihood = zeros(n, 1);
    
    for i = 1:n
        t = rt(i) - T;
        if t <= 0
            logLikelihood(i) = -Inf; % Non-decision time cannot be greater than RT
            continue;
        end
        
        prob = utl_wfpt(t, -v, a, 1-w)+ ...
            utl_wfpt(t, v, a, w);


        
        % Avoid log(0) by setting a small lower bound for probabilities
        prob = max(prob, 1e-10);
        logLikelihood(i) = log(prob);
    end
    
    negLogLikelihood = -sum(logLikelihood);
end

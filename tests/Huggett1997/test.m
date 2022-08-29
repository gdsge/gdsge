%% Compile
gdsge_codegen('huggett1997');

%% Solve a WarmUp problem
IterRslt = iter_huggett1997;

%% A fixed-point loop to solve the initial steady state
tolEq = 1e-5; metric = inf; iter = 0;
UPDATE_SPEED = 0.01;
K = IterRslt.var_others.kSs;
alpha = IterRslt.var_others.alpha; delta = IterRslt.var_others.delta;
% Non-stochastic simulation, prepare distribution grid
kFinePts = 1000; shockPts = IterRslt.shock_num;
kFine = linspace(min(IterRslt.var_state.k),max(IterRslt.var_state.k),kFinePts)';
kFineRight = [kFine(2:end);inf];
[kFineGrid,shockGrid] = ndgrid(kFine,1:shockPts);
% Parameters to simulate only one step
simuOptions.num_periods = 1;
simuOptions.num_samples = numel(kFineGrid);
simuOptions.init.k = kFineGrid(:);
simuOptions.init.shock = shockGrid(:);
while metric > tolEq
    % Solve at prices implied by current K
    options = struct;
    options.TASK = 'ss';
    options.r = alpha*K^(alpha-1) - delta;
    options.w = (1-alpha)*K^alpha;
    options.WarmUp = IterRslt;
    IterRslt = iter_huggett1997(options);
 
    % Non-stochastic simulation. Simulate one-step to get the state transition
    % functions over kFine
    SimuRslt = simulate_huggett1997(IterRslt,simuOptions);
    % Construct the Markov transition implied by policy functions
    kp = SimuRslt.k(:,2);
    [~,kpCell] = histc(kp, [kFine;inf]);
    leftWeights = (kFineRight(kpCell)-kp) ./ (kFineRight(kpCell)-kFine(kpCell));
    leftWeights(kpCell>=kFinePts) = 1;
    rowVec = [1:shockPts*kFinePts]';
    transToKp = sparse(rowVec,kpCell,leftWeights,shockPts*kFinePts,kFinePts) ...
        + sparse(rowVec(kpCell<kFinePts),kpCell(kpCell<kFinePts)+1,1-leftWeights(kpCell<kFinePts),...
        shockPts*kFinePts,kFinePts);
    % Accomodate the exogenous transition
    transFull = repmat(transToKp,[1,2]) * 0.5;
    % Simulate
    [stationaryDist,~] = eigs(transFull',1,1);
    stationaryDist = reshape(stationaryDist / sum(stationaryDist(:)),[kFinePts,shockPts]);
    % Statistics
    K_new = sum(reshape(stationaryDist.*reshape(kFine,[kFinePts,1]), [], 1));

    % Update
    metric = abs(log(K) - log(K_new));
    iter = iter + 1;
    fprintf('Steady-state iterations: %d, %g\n',iter, metric);
    fprintf('===============================\n');
    K = K_new*UPDATE_SPEED + K*(1-UPDATE_SPEED);
end

%% Solve the transition path
T = 1000;
K_t = K*ones(1,T);
K_t_new = K*ones(1,T);
tolEq = 1e-3; metric = inf; iter = 0;
UPDATE_SPEED = 0.01;
% Initial distribution in Huggett (1997)
dist0 = stationaryDist;
dist0(1,:) = 0.2/2;
kBar = K/0.8*2;
kBarIndex = find(kFine>kBar,1);
dist0(2:kBarIndex,:) = 0.8 / numel(dist0(2:kBarIndex,:));
dist0(kBarIndex+1:end,:) = 0;
while metric > tolEq
    % Backward loop
    options = struct;
    options.TASK = 'transition';
    options.PrintFreq = inf;
    options.MaxIter = T;
    options.T = T;
    options.TolEq = 0;  % Do not check TolEq
    options.r_t = alpha*K_t.^(alpha-1) - delta;
    options.w_t = (1-alpha)*K_t.^alpha;
    options.WarmUp = IterRslt;  % Start from steady state
    options.WarmUp.Iter = 0;    % Start with iter 0;
    IterRslt_t = iter_huggett1997(options);

    % Forward simulation
    dist = dist0;
    for t=1:1:T
        K_t_new(t) = sum(reshape(dist.*reshape(kFine,[kFinePts,1]), [], 1));
        % Simulate using period-t policies
        IterRslt.output_interp = IterRslt_t.var_others.output_interp_t{t};
        SimuRslt_t = simulate_huggett1997(IterRslt,simuOptions);
        % Construct the Markov transition implied by policy functions
        kp = SimuRslt_t.k(:,2);
        [~,kpCell] = histc(kp, [kFine;inf]);
        leftWeights = (kFineRight(kpCell)-kp) ./ (kFineRight(kpCell)-kFine(kpCell));
        leftWeights(kpCell>=kFinePts) = 1;
        rowVec = [1:shockPts*kFinePts]';
        transToKp = sparse(rowVec,kpCell,leftWeights,shockPts*kFinePts,kFinePts) ...
            + sparse(rowVec(kpCell<kFinePts),kpCell(kpCell<kFinePts)+1,1-leftWeights(kpCell<kFinePts),...
            shockPts*kFinePts,kFinePts);
        % Accomodate the exogenous transition
        transFull = [transToKp,transToKp] * 0.5;
        dist = reshape(dist(:)'*transFull,[kFinePts,shockPts]);
    end
    
    % Update K_t
    metric = max(abs(log(K_t) - log(K_t_new)));
    iter = iter + 1;
    fprintf('Transition path iterations: %d, %g\n',iter, metric);
    fprintf('==================================\n');
    if metric<2e-2
        UPDATE_SPEED = 0.03;
    end
    K_t = K_t_new*UPDATE_SPEED + K_t*(1-UPDATE_SPEED);
end

%% Plot
figure; hold on;
plot(K_t(1:500),'LineWidth',1.5);
plot([0,500],[K,K],'--','LineWidth',1.5);
title('Transition Path');
xlabel('TIME');
legend({'Equilibrium Path','Steady State Path'});
ylabel('Aggregate Capital Stock');
print('transition_path.png','-dpng');

figure; hold on;
plot(kFine,stationaryDist,'LineWidth',1.5);
title('Stationary Distribution');
xlim([0,20]);
xlabel('Capital');
ylabel('Fractions');
print('stationary_dist.png','-dpng');

figure; hold on;
plot(IterRslt.var_state.k,IterRslt.var_policy.k_next,'LineWidth',1.0);
plot(IterRslt.var_state.k,IterRslt.var_state.k,'k-');
xlim([0,10]);
ylim([0,10]);
legend('a''(k,0.8)','a''(k,1.2)','45 Degreee Line','Location','SouthEast');
xlabel('Capital');
ylabel('Capital Next Period');
title('Decision Rules for Saving');
print('policy_function_kp.png','-dpng');


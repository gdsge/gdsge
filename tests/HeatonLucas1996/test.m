gdsge_codegen('HL1996');
options = struct;
options.SaveFreq = inf;
IterRslt = iter_HL1996(options);
rng(0823);
SimuRslt = simulate_HL1996(IterRslt);

% Policy functions
%{
figure;
plot(IterRslt.var_state.w1, IterRslt.var_policy.c1,'LineWidth',1.5);
%}

figure; subplot(3,1,1);
plot(IterRslt.var_state.w1, IterRslt.var_policy.ms1(1,:),'LineWidth',1.5);
title('Multiplier for No-short Constraint, Agent 1');
subplot(3,1,2);
plot(IterRslt.var_state.w1, IterRslt.var_policy.mb1(1,:),'LineWidth',1.5);
title('Multiplier for Borrowing Constraint, Agent 1');
subplot(3,1,3);
plot(IterRslt.var_state.w1, IterRslt.var_aux.equity_premium(1,:)*100,'LineWidth',1.5);
title('Equity Premium');
xlabel('Wealth Share of Agent 1');

figure;
histogram(SimuRslt.w1(:,1000:end),50,'Normalization','probability');
title('Histogram of Wealth Share in the Ergodic Distribution');
xlabel('Wealth Share of Agent 1');
ylabel('Fractions');
% export_fig('figures/histogram_w1.png','-r300','-transparent');

IterRslt = load('IterRslt_HL1996_209.mat');
IterRslt = IterRslt.IterRslt;

figure;
plot(IterRslt.var_state.w1, IterRslt.var_aux.equity_premium*100,'LineWidth',1.5);
title('Equity Premium');
xlabel('Wealth Share of Agent 1');
ylabel('%');
% export_fig('figures/policy_equity_premium.png','-transparent');

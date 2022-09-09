gdsge_codegen('mendoza2010');

IterRslt = iter_mendoza2010;
% save('IterRslt.mat','IterRslt');

% Resolve with finer grids
options = struct;
cTilde_pts_fine = 80;
options.cTilde = linspace(IterRslt.var_state.cTilde(1),IterRslt.var_state.cTilde(end),...
    cTilde_pts_fine);
k_pts_fine = 80;
options.k = linspace(IterRslt.var_state.k(1),IterRslt.var_state.k(end),k_pts_fine);
options.WarmUp = IterRslt;
IterRslt = iter_mendoza2010(options);
% save('IterRslt.mat','IterRslt');

% load('IterRslt.mat');
SimuRslt = simulate_mendoza2010(IterRslt);

% Some moments
statsList = {'mean','std'};
varList = {'gdp','c','inv','nx_gdp','k','b_gdp','q','lev','v','wkcptl'};
ratioList = {'nx_gdp','b_gdp','lev'};
meanStats = struct;
stdStats = struct;
corrStats = struct;
autoCorrStats = struct;
startPeriod = 1001;
endPeriod = 50000;
for i_var = 1:length(varList)
    varName = varList{i_var};
    meanStats.(varName) = mean(reshape(SimuRslt.(varName)(:,startPeriod:endPeriod),[],1));
    if max(strcmp(varName,ratioList))==1
        stdStats.(varName) = std(reshape(SimuRslt.(varName)(:,startPeriod:endPeriod),[],1));
    else
        stdStats.(varName) = std(reshape(SimuRslt.(varName)(:,startPeriod:endPeriod),[],1)) ...
            / meanStats.(varName);
    end
    corrStats.(varName) = corr(reshape(SimuRslt.(varName)(:,startPeriod:endPeriod),[],1), ...
        reshape(SimuRslt.gdp(:,startPeriod:endPeriod),[],1));
    autoCorrStats.(varName) = corr(reshape(SimuRslt.(varName)(1,startPeriod+1:endPeriod),[],1), ...
        reshape(SimuRslt.(varName)(1,startPeriod:endPeriod-1),[],1));
end

% Fraction in crisis
muGreaterThanZero = SimuRslt.mu > 1e-7;
fractionCrisis = sum(muGreaterThanZero(:)) / numel(muGreaterThanZero(:));
fprintf('The fraction in crisis: %g\n',fractionCrisis);

% Write to tables
% writetable(struct2table(meanStats),'tables/mean_stats.xlsx');
% writetable(struct2table(stdStats),'tables/std_stats.xlsx');
% writetable(struct2table(corrStats),'tables/corr_stats.xlsx');
% writetable(struct2table(autoCorrStats),'tables/autoCorr_stats.xlsx');

figure;
hist3([SimuRslt.b(:),reshape(SimuRslt.k(:,1:end-1),[],1)],[50,50],'CDataMode','auto','FaceColor','interp');
view([15,30]);
xlabel('$b$','interpreter','latex','FontSize',16);
ylabel('$k$','interpreter','latex','FontSize',16);
title('Histogram of Capital and Bond','interpreter','latex','FontSize',14);
% print('figures/histogram_b_k.png','-dpng');

% Policy function
% State transition function
figure;
[cTilde_mesh,k_mesh] = ndgrid(IterRslt.var_state.cTilde,IterRslt.var_state.k);
surf(cTilde_mesh,k_mesh,squeeze(IterRslt.var_aux.b(1,:,:)));
xlabel('$\tilde{c}$','interpreter','latex','FontSize',16);
ylabel('$k$','interpreter','latex','FontSize',16);
zlabel('$b$','interpreter','latex','FontSize',16);
title('State Transformation of Current Bond');
% export_fig('figures/policy_b_original_state.png','-r300');


% Transformmed state space
figure;
[cTilde_mesh,k_mesh] = ndgrid(IterRslt.var_state.cTilde,IterRslt.var_state.k);
surf(squeeze(IterRslt.var_aux.b(1,:,:)),k_mesh,cTilde_mesh);
xlabel('$b$','interpreter','latex','FontSize',16);
ylabel('$k$','interpreter','latex','FontSize',16);
% zlabel('$\tilde{c}$','interpreter','latex','FontSize',16);
title('Policy Function for Consumption-labor Bundle $\tilde{c}$','interpreter','latex','FontSize',14);
xlim([-200,500]);
zlim([100,250]);
% view(-0.06,90);
% print('figures/policy_cTilde_feasible.png','-dpng','-r300');

%%%%%%%% Nonlinearality

figure;
[cTilde_mesh,k_mesh] = ndgrid(IterRslt.var_state.cTilde,IterRslt.var_state.k);
surf(squeeze(IterRslt.var_aux.b(1,:,:)),k_mesh,squeeze(IterRslt.var_aux.bNext(1,:,:)));
xlabel('b');
ylabel('k');
title('bNext');
xlim([-200,500]);
view(-0.06,30);
% print('figures/policy_bNext.png');

figure;
[cTilde_mesh,k_mesh] = ndgrid(IterRslt.var_state.cTilde,IterRslt.var_state.k);
surf(squeeze(IterRslt.var_aux.b(1,:,:)),k_mesh,squeeze(IterRslt.var_aux.q(1,:,:)));
xlabel('$b$','interpreter','latex','FontSize',16);
ylabel('$k$','interpreter','latex','FontSize',16);
xlim([-200,500]);
title('Capital Price, $q$','interpreter','latex','FontSize',14);
% print('figures/policy_q.png','-dpng','-r300');

figure;
[cTilde_mesh,k_mesh] = ndgrid(IterRslt.var_state.cTilde,IterRslt.var_state.k);
surf(squeeze(IterRslt.var_aux.b(1,:,:)),k_mesh,squeeze(IterRslt.var_policy.mu(1,:,:)));
xlabel('$b$','interpreter','latex','FontSize',16);
ylabel('$k$','interpreter','latex','FontSize',16);
xlim([-200,500]);
title('Multiplier for the Collateral Constraint, $\mu$','interpreter','latex','FontSize',14);
view(-16,26);
% print('figures/policy_mu.png','-dpng','-r300');

figure;
[cTilde_mesh,k_mesh] = ndgrid(IterRslt.var_state.cTilde,IterRslt.var_state.k);
surf(squeeze(IterRslt.var_aux.b(1,:,:)),k_mesh,squeeze(IterRslt.var_aux.Y(1,:,:)));
xlabel('b');
ylabel('k');
title('Output');
view(-0.06,42.8);
% print('figures/policy_Y.eps');

%%%%%%% Histogram
firstPeriod = 5001;

figure; hold on;
histogram(SimuRslt.k(:,firstPeriod:end),'Normalization','pdf');
[density,grid] = ksdensity(reshape(SimuRslt.k(:,firstPeriod:end),1,[]));
plot(grid,density,'r-','LineWidth',2);
ylabel('Probability');
title('Histogram of k');
% print('figures/histogram_K.eps');

figure; hold on;
histogram(SimuRslt.b(:,firstPeriod:end),'Normalization','pdf');
[density,grid] = ksdensity(reshape(SimuRslt.b(:,firstPeriod:end),1,[]));
plot(grid,density,'r-','LineWidth',2);
ylabel('Probability');
title('Histogram of b');
% print('figures/histogram_b.eps');


figure; hold on;
histogram(SimuRslt.mu(:,firstPeriod:end),'Normalization','pdf');
% [density,grid] = ksdensity(reshape(simuRslt.mu(:,firstPeriod:end),1,[]));
% plot(grid,density,'r-','LineWidth',2);
ylabel('Probability');
title('Histogram of multiplier');
% print('figures/histogram_mu.eps');


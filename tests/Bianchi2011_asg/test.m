gdsge_codegen('bianchi2011');

options = struct;
shock_process = load('shock_process.mat');
options.shock_trans = shock_process.shock_trans;
options.yT = shock_process.yT;
options.yN = shock_process.yN;
options.MaxIter = 50;
IterRslt = iter_bianchi2011(options);

options.MaxIter = inf;
options.WarmUp = IterRslt;
options.bMin = -1.1;
options.bMax = 0.0;
options.b = [options.bMin,options.bMax];
options.SkipModelInit = 1;
IterRslt = iter_bianchi2011(options);

save('IterRslt.mat','IterRslt');

SimuRslt = simulate_bianchi2011(IterRslt);

GNDSGE_ASG_INTERP = asg.construct_from_struct(IterRslt.asg_output_struct);
grids = GNDSGE_ASG_INTERP.get_grids_info;
bNext_idx = 1; % index in var_output
pN_idx = 2;
for j=1:16
    grid = grids{j};
    lenGrid = length(grid);
    bNext_fval{j} = GNDSGE_ASG_INTERP.eval(j*ones(1,lenGrid),bNext_idx*ones(1,lenGrid),grid);
    pN_fval{j} = GNDSGE_ASG_INTERP.eval(j*ones(1,lenGrid),pN_idx*ones(1,lenGrid),grid);
end

figure; 
subplot(2,1,1); hold on;
xy = sortrows([grids{1}',bNext_fval{1}']);
plot(xy(:,1),xy(:,2),'ro-');
xy = sortrows([grids{4}',bNext_fval{4}']);
plot(xy(:,1),xy(:,2),'kx-');
title('Policy Functions for Next Period Bond Holding, $b''$','interpreter','latex','FontSize',12);
xlabel('Current Bond Holding, $b$','FontSize',12,'interpreter','latex');

subplot(2,1,2); hold on;
xy = sortrows([grids{1}',pN_fval{1}']);
plot(xy(:,1),xy(:,2),'ro-');
xy = sortrows([grids{4}',pN_fval{4}']);
plot(xy(:,1),xy(:,2),'kx-');
title('Policy Functions for Non-tradable Good Price, $p^N$','interpreter','latex','FontSize',12);
xlabel('Current Bond Holding, $b$','FontSize',12,'interpreter','latex');
legend({'$y_t^T$ Lowest, $y_t^N$ Lowest','$y_t^T$ Highest, $y_t^N$ Lowest'},'Location','SouthEast','interpreter','latex','FontSize',12);
% print('figures/policy_combined.png','-dpng');

%%%%%%%%%%%%%%%%% Histogram of states %%%%%%%%%%%%%%%%
figure; hold on;
hh = histogram(SimuRslt.b(:,500:end),50,'Normalization','probability');
[density,grid] = ksdensity(reshape(SimuRslt.b(:,500:end),1,[]));
save('density.mat','density','grid');
density = density*hh.BinWidth;
plot(grid,density,'r-','LineWidth',2);
title('Histogram and Kernel Density of Bond Holding','interpreter','latex','FontSize',12);
xlabel('Bond Holding, $b$','FontSize',12,'interpreter','latex');
ylabel('Fractions','interpreter','latex');
xlim([-1.1,-0.5]);
% print('figures/histogram_b.png','-dpng');



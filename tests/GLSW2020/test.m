gdsge_codegen('GLSW_interp');
IterRslt = iter_GLSW_interp;

options = struct;
options.SimuPrintFreq = 1000;
rng(0823);
SimuRslt = simulate_GLSW_interp(IterRslt,options);

% Get Ergodi Set
startIdx = 5000;
endIdx = 10000;
shock = SimuRslt.shock(:,startIdx:endIdx);
shockOn = shock;
shockOn(:) = 2;
a1 = SimuRslt.a1(:,startIdx:endIdx);

lenSamples = numel(shock);
options.init.shock = [shock(:);shockOn(:)];
options.init.a1 = [a1(:);a1(:)];

% Simulate forward
options.num_samples = lenSamples*2;
options.num_periods = 100;
SimuRsltIrf = simulate_GLSW_interp(IterRslt,options);

% Construct and plot the IRF
irf.r = mean(SimuRsltIrf.r(lenSamples+1:2*lenSamples,:) - SimuRsltIrf.r(1:lenSamples,:));
irf.a1 = mean(SimuRsltIrf.a1(lenSamples+1:2*lenSamples,:) - SimuRsltIrf.a1(1:lenSamples,:));

figure; subplot(2,1,1);
plot(irf.r(1:10)*100,'LineWidth',2);
xlabel('Quarters','interpreter','latex','FontSize',12);
ylabel('%');
title('Impulse Response of Interest Rate to a Pandemic Shock','interpreter','latex','FontSize',11);

subplot(2,1,2);
plot(irf.a1(1:100),'LineWidth',2);
xlabel('Quarters','interpreter','latex','FontSize',12);
ylabel('Level');
title('Impulse Response of Contact-Intensive Workers'' Wealth to a Pandemic Shock','interpreter','latex','FontSize',11);
print('irf_r_a1.png','-dpng','-r300');


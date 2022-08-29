gdsge_codegen('CaoKS2016');

IterOptions.PrintFreq = 10;
IterOptions.SaveFreq = 100;

IterRslt = iter_CaoKS2016(IterOptions);
SimuRslt = simulate_CaoKS2016(IterRslt);

GNDSGE_ASG_INTERP = asg.construct_from_struct(IterRslt.asg_interp_struct);

[grids,surplus] = GNDSGE_ASG_INTERP.get_grids_info;

for j=1:4
    grid = grids{j};
    fval{j} = GNDSGE_ASG_INTERP.eval_vec(j*ones(1,size(grid,2)),grid);
end

figure;hold on;
for j=1:1
    scatter3(grids{j}(1,:),grids{j}(2,:),fval{j}(1,:));
end
legend('Shock\_index=1');
xlabel('K');
ylabel('X');
title('c1Future');
view(-122.8,42);
% export_fig('figures/c1Future.eps');

figure;hold on;
for j=1:1
    scatter3(grids{j}(1,:),grids{j}(2,:),fval{j}(2,:));
end
legend('Shock\_index=1');
xlabel('K');
ylabel('X');
title('c2Future');
view(-122.8,42);
% export_fig('figures/c2Future.eps');

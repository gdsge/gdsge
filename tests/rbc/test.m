% Global var_aux demo: MPK = GDSGE_DER(Y, K) is the spline derivative of
% output w.r.t. capital, constructed after each iteration; MPK_analytical
% is the closed form alpha*z*K^(alpha-1), so the two must coincide.
gdsge_codegen('rbc');
IterRslt = iter_rbc;
err = max(abs(IterRslt.var_aux.MPK(:) - IterRslt.var_aux.MPK_analytical(:)));
fprintf('MPK (GDSGE_DER) vs closed form: max abs err = %.3e\n', err);
if usejava('desktop')   % demo figure: interactive (GUI) sessions only
    figure; hold on;
    plot(IterRslt.var_state.K, IterRslt.var_aux.MPK(1,:), '-');
    plot(IterRslt.var_state.K, IterRslt.var_aux.MPK_analytical(1,:), '--');
    legend({'MPK (GDSGE\_DER)', 'MPK analytical'});
    xlabel('K'); title('Global var\_aux: spline derivative vs closed form');
end

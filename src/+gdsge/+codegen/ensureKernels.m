function ensureKernels(ir)
% ENSUREKERNELS  Compile (cache-gated) the interpolation kernels this model's C++
%   backend needs: asg_mex for ASG models; the fused spline constructor + the generic
%   evaluator for cartesian spline/linear models. Extracted from generateCxx so the
%   codegen driver can announce a "Compiling kernels" phase and reuse the same logic.
if strcmp(ir.options.interpMethod, 'asg')
    gdsge.codegen.ensureAsgMex();
else
    gdsge.codegen.ensureSplineConstructMex();
    gdsge.codegen.ensureInterpEvalMex();
end
end

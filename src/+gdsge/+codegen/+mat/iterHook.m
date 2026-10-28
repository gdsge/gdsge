function code = iterHook(ir, which)
% ITERHOOK  The gmod's pre_iter / post_iter MATLAB text (WHICH = 'pre' |
%   'post'), '' when not declared, spliced into the iter loop at the old
%   template's PRE_ITER_CODE / POST_ITER_CODE positions. The old placeholder
%   OUTPUT_CONSTRUCT_CODE inside a hook expands to the output-interp
%   construction (spline only).
f = [which 'Iter'];
code = '';
if ~isfield(ir, 'hooks') || ~isfield(ir.hooks, f) || isempty(ir.hooks.(f)); return; end
code = ir.hooks.(f);
tok = '(?<![A-Za-z0-9_])OUTPUT_CONSTRUCT_CODE(?![A-Za-z0-9_])';
if ~isempty(regexp(code, tok, 'once'))
    if strcmp(ir.options.interpMethod, 'asg')
        error('gdsge:codegen:outputConstructAsgUnsupported', ...
            'OUTPUT_CONSTRUCT_CODE in a %s_iter hook is supported for spline models only', which);
    end
    rep = gdsge.codegen.mat.emitOutputConstruct(ir);
    code = regexprep(code, tok, strrep(strrep(rep, '\', '\\'), '$', '\$'));
end
end

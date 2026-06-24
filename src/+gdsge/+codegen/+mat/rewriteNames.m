function out = rewriteNames(txt, names, replacements)
% REWRITENAMES  Word-boundary identifier replacement in opaque MATLAB text.
%   Replaces whole-identifier occurrences of names{i} with replacements{i}.
%   Boundaries are MATLAB identifier characters, so 'w1' does not match inside
%   'w1n' or 'Kw1'. Replacements containing identifiers are safe because the
%   inserted text always embeds the original name with a '_' prefix boundary.
out = txt;
for i = 1:numel(names)
    pat = ['(?<![A-Za-z0-9_])' regexptranslate('escape', names{i}) '(?![A-Za-z0-9_])'];
    rep = strrep(replacements{i}, '\', '\\');
    rep = strrep(rep, '$', '\$');
    out = regexprep(out, pat, rep);
end
end

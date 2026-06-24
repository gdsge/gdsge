function txt = fillTemplate(txt, pairs)
% FILLTEMPLATE  Sequential word-boundary placeholder substitution.
%   pairs: N-by-2 cell {placeholderName, filledText; ...}, applied in order
%   (assemble innermost templates first so later keys never hit filled text
%   unintentionally). Word boundary treats '_' as a word character, so
%   NUM_EQUATIONS does not match inside MY_NUM_EQUATIONS. Errors if a key
%   matches nothing: gdsge:codegen:placeholderNotFound.
if numel(unique(pairs(:, 1))) ~= size(pairs, 1)
    error('gdsge:codegen:duplicateKey', 'duplicate placeholder keys in pairs');
end
nonWord = '[^a-zA-Z0-9_]';
for i = 1:size(pairs, 1)
    key = pairs{i, 1};
    pat = ['((?<=^)|(?<=' nonWord '))' key '((?=$)|(?=' nonWord '))'];
    if isempty(regexp(txt, pat, 'once'))
        error('gdsge:codegen:placeholderNotFound', ...
            'placeholder %s not found in template', key);
    end
    rep = strrep(strrep(pairs{i, 2}, '\', '\\'), '$', '\$');   % regexprep escapes
    txt = regexprep(txt, pat, rep);
end
end

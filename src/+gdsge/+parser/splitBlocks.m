function out = splitBlocks(cleanText)
% SPLITBLOCKS  Partition preprocessed gmod into declaration lines and named
%   blocks. Returns struct with:
%     .declText           the non-block lines (joined with \n)
%     .blocks.<name>      raw text of each top-level block (opener/closer removed)
%   Block openers (depth-tracked, so nested equations; ... end; is retained
%   inside model): model simulate model_init pre_model pre_iter post_iter
%   pre_jac_code post_jac_code equations. Closer: end; .
%   NOTE (Phase-7c hardening): a standalone MATLAB control-flow `end;` inside a
%   hook/model body would also match the closer; the core models do not hit this.
openers = {'model','simulate','model_init','pre_model','pre_iter', ...
           'post_iter','pre_jac_code','post_jac_code','equations'};
lines = regexp(cleanText, '\n', 'split');
declLines = {};
blocks = struct();
modelRegions = {};
stack = {};
curTop = '';
curCond = '';
curBody = {};
for i = 1:numel(lines)
    raw = lines{i};
    t = strtrim(raw);
    opener = matchOpener(t, openers);
    isEnd = strcmp(t, 'end;');
    if isempty(stack)
        if ~isempty(opener)
            stack{end+1} = opener; curTop = opener; curBody = {}; %#ok<AGROW>
            if strcmp(opener, 'model')
                cm = regexp(regexprep(t, ';\s*$', ''), '^model\((.*)\)$', 'tokens', 'once');
                if isempty(cm); curCond = ''; else; curCond = cm{1}; end
            end
        elseif isEnd
            error('gdsge:parser:unterminatedBlock', 'Unexpected "end;" at line %d', i);
        else
            declLines{end+1} = raw; %#ok<AGROW>
        end
    else
        if ~isempty(opener)
            stack{end+1} = opener; curBody{end+1} = raw; %#ok<AGROW>
        elseif isEnd
            stack(end) = [];
            if isempty(stack)
                if strcmp(curTop, 'model')
                    modelRegions{end+1} = struct('condition', curCond, ...
                        'body', strjoin(curBody, sprintf('\n'))); %#ok<AGROW>
                else
                    blocks.(curTop) = strjoin(curBody, sprintf('\n'));
                end
                curTop = ''; curCond = ''; curBody = {};
            else
                curBody{end+1} = raw; %#ok<AGROW>
            end
        else
            curBody{end+1} = raw; %#ok<AGROW>
        end
    end
end
if ~isempty(stack)
    error('gdsge:parser:unterminatedBlock', 'Block "%s" not closed by "end;"', stack{end});
end
out = struct('declText', strjoin(declLines, sprintf('\n')), ...
    'blocks', blocks, 'modelRegions', {modelRegions});
end

function name = matchOpener(t, openers)
name = '';
tt = regexprep(t, ';\s*$', '');
if ismember(tt, openers)
    name = tt;
elseif ~isempty(regexp(tt, '^model\(.*\)$', 'once'))
    name = 'model';           % conditional region: model(<cond>)
elseif ~isempty(regexp(tt, '^model_init\(.*\)$', 'once'))
    error('gdsge:parser:conditionalModelInitUnsupported', ...
        'conditional model_init regions are not supported');
end
end

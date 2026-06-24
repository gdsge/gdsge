function md = gendoc(targetPath)
% GENDOC  Render the IR schema descriptor as Markdown.
%   md = gdsge.ir.gendoc()        returns the markdown string.
%   gdsge.ir.gendoc(targetPath)   also writes it to TARGETPATH (utf-8, LF).
s = gdsge.ir.schema();
L = {};
L{end+1} = '# IR Schema (generated — do not edit by hand)';
L{end+1} = '';
L{end+1} = sprintf('Generated from `src/+gdsge/+ir/schema.m` by `gdsge.ir.gendoc`. irVersion: `%s`.', s.irVersion);
L{end+1} = '';
L{end+1} = '## Document sections';
L{end+1} = '';
L = renderStruct(L, s.root, 0);
L{end+1} = '';
L{end+1} = '## AST nodes';
L{end+1} = '';
L = renderRegistry(L, s.nodes);
L{end+1} = '';
L{end+1} = '## Model statements';
L{end+1} = '';
L = renderRegistry(L, s.stmts);
L{end+1} = '';
L{end+1} = '## Equation kinds';
L{end+1} = '';
L = renderRegistry(L, s.eqkinds);
md = strjoin(L, sprintf('\n'));
if nargin >= 1 && ~isempty(targetPath)
    fid = fopen(targetPath, 'w', 'n', 'UTF-8');
    if fid < 0; error('gdsge:ir:gendoc:open', 'cannot write %s', targetPath); end
    fwrite(fid, md);
    fclose(fid);
end
end

function L = renderRegistry(L, reg)
fn = fieldnames(reg);
for i = 1:numel(fn)
    L{end+1} = sprintf('### `%s`', fn{i}); %#ok<AGROW>
    L{end+1} = ''; %#ok<AGROW>
    L = renderStruct(L, reg.(fn{i}), 0);
    L{end+1} = ''; %#ok<AGROW>
end
end

function L = renderStruct(L, spec, depth)
fn = fieldnames(spec.fields);
for i = 1:numel(fn)
    nm = fn{i}; fs = spec.fields.(nm);
    L{end+1} = line(depth, nm, fs); %#ok<AGROW>
    L = renderChildren(L, fs, depth + 1);
end
end

function L = renderChildren(L, fs, depth)
switch fs.kind
    case 'struct'
        L = renderStruct(L, fs, depth);
    case 'list'
        if strcmp(fs.item.kind, 'struct')
            L{end+1} = pad(depth, '- *(each item)*'); %#ok<AGROW>
            L = renderStruct(L, fs.item, depth + 1);
        end
    case 'map'
        if strcmp(fs.value.kind, 'struct')
            L{end+1} = pad(depth, '- *(each value)*'); %#ok<AGROW>
            L = renderStruct(L, fs.value, depth + 1);
        end
end
end

function s = line(depth, nm, fs)
extra = '';
switch fs.kind
    case 'enum'; extra = sprintf(' {%s}', strjoin(fs.values, ', '));
    case 'ref';  extra = sprintf(' -> pool:%s', fs.pool);
    case 'list'
        if strcmp(fs.item.kind, 'ref')
            extra = sprintf(' of ref -> pool:%s', fs.item.pool);
        else
            extra = sprintf(' of %s', fs.item.kind);
        end
    case 'map'
        if strcmp(fs.value.kind, 'ref')
            extra = sprintf(' of ref -> pool:%s', fs.value.pool);
        else
            extra = sprintf(' of %s', fs.value.kind);
        end
    case {'tagged','taggedList'}; extra = sprintf(' (%s)', fs.registry);
end
req = 'required'; if isfield(fs,'required') && ~fs.required; req = 'optional'; end
s = pad(depth, sprintf('- `%s` : %s%s (%s)', nm, fs.kind, extra, req));
end

function s = pad(depth, txt)
s = [repmat('  ', 1, depth), txt];
end

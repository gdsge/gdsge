function txt = readTemplate(name)
% READTEMPLATE  Read templates/cxx/<name> (repo-root-relative, resolved from
%   this file's location: src/+gdsge/+codegen/+cxx -> repo root is 4 up).
here = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(fileparts(fileparts(here))));
p = fullfile(repoRoot, 'templates', 'cxx', name);
if ~exist(p, 'file')
    error('gdsge:codegen:templateMissing', 'template not found: %s', p);
end
txt = fileread(p);
end

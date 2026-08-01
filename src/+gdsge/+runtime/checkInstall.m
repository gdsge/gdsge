function checkInstall(srcDir)
% CHECKINSTALL  Verify every toolbox source file actually declares a function.
%   A truncated or empty .m file makes MATLAB classify it as a script; calling
%   it then fails with MATLAB:scriptNotAFunction naming an internal file, which
%   tells the user nothing about what is wrong or how to fix it. This turns
%   that into one actionable error naming every damaged file.
%
%   checkInstall()        scan the installed tree, once per MATLAB session.
%   checkInstall(srcDir)  scan srcDir unconditionally (tests).
%
%   Contents.m is exempt: those files are comment-only by design. The scan
%   costs a few milliseconds over ~70 small files and runs from ensurePath,
%   ahead of any generated iter_/simulate_ work.
persistent checked
if nargin < 1
    if ~isempty(checked); return; end
    checked = true;    % set before scanning: a failure must not re-scan forever
    srcDir = fileparts(fileparts(fileparts(mfilename('fullpath'))));   % -> src/
end

files = [ dir(fullfile(srcDir, '*.m'));
          dir(fullfile(srcDir, 'kernels', '*.m'));
          dir(fullfile(srcDir, '+gdsge', '**', '*.m')) ];

bad = {};
for i = 1:numel(files)
    if strcmp(files(i).name, 'Contents.m'); continue; end
    p = fullfile(files(i).folder, files(i).name);
    txt = fileread(p);
    % \> is MATLAB's end-of-word anchor (\b is a backspace here, not a
    % boundary). ^\s* cannot skip a leading %, so a comment that merely
    % contains the word "function" does not count as a declaration.
    if isempty(regexp(txt, '^\s*(function|classdef)\>', 'once', 'lineanchors'))
        bad{end+1} = p; %#ok<AGROW>
    end
end

if ~isempty(bad)
    error('gdsge:install:corrupt', ...
        ['GDSGE install looks incomplete or corrupted - these files declare no ' ...
         'function:\n  %s\n' ...
         'Restore them (e.g. `git checkout -- src/`) or re-clone the toolbox.'], ...
        strjoin(bad, sprintf('\n  ')));
end
end

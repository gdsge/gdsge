function [ir, gen] = codegen(modelName, options)
% CODEGEN  Public codegen driver: <model>.gmod in cwd -> generated files + MEX.
%   Writes iter_<model>.m, simulate_<model>.m, mex_<model>.cpp,
%   compile_<model>.m and <model>.gdsge.json into the current directory, then
%   compiles the MEX — skipped when the C++ matches mex_<model>.cache. The
%   cache is written only after a successful compile, so failed builds retry.
%   Returns the IR struct. With a second output, also returns gen, a struct of
%   the generated code strings (iterCode/simulateCode/cppCode/compileCode), the
%   pre-overwrite cache content (cppCache), and codeSegment (the real generated
%   strings — not the old template fragments). Called by the flat shim
%   gdsge_codegen, which re-exposes these as its legacy positional outputs.
if nargin < 2; options = []; end
validateOptions(options);

gmodFile = fullfile(pwd, [modelName '.gmod']);
if ~exist(gmodFile, 'file')
    error('gdsge:codegen:gmodNotFound', ...
        '%s.gmod not found in %s', modelName, pwd);
end

% Register src/kernels on the path up front so a single `addpath('src')` suffices:
% ASG codegen calls asg.get_mex_constants() (-> asg_mex) inside generateCxx/emitCompile,
% well before the compile step. ensure*Mex build the kernels but do not addpath them.
gdsge.runtime.ensurePath();

TOTAL = 5;
banner(1, TOTAL, sprintf('Parsing gmod (%s)', modelName));
ir = gdsge.parser.parseFrontEnd(fileread(gmodFile), modelName, fileparts(gmodFile));
gdsge.codegen.writeText(fullfile(pwd, [modelName '.gdsge.json']), gdsge.ir.encode(ir));

banner(2, TOTAL, 'Generating MATLAB (iter + simulate)');
matFiles = gdsge.codegen.generateMatlab(ir, pwd);

banner(3, TOTAL, 'Compiling kernels (cache-gated)');
gdsge.codegen.ensureKernels(ir);

dec = gdsge.codegen.resolveBackend(ir);
banner(4, TOTAL, sprintf('Generating C++ (mex_%s.cpp)', modelName));
fprintf('      %s\n', gdsge.codegen.backendMessage(dec));
files = gdsge.codegen.generateCxxWithFallback(ir, pwd, dec);

cppText   = fileread(files.cppFile);
if isfield(files, 'simuMexFile') && ~isempty(files.simuMexFile)
    cppText = [cppText, fileread(files.simuMexFile)];   % recompile if simulate MEX changed
end
cacheFile = fullfile(pwd, ['mex_' modelName '.cache']);
cppCache  = '';
if exist(cacheFile, 'file'); cppCache = fileread(cacheFile); end
if gdsge.codegen.needsCompile(cppText, cacheFile)
    banner(5, TOTAL, sprintf('Compiling solver MEX (mex_%s)', modelName));
    feval(['compile_' modelName]);   % ensurePath already ran at the top of codegen
    gdsge.codegen.writeText(cacheFile, cppText);
else
    banner(5, TOTAL, sprintf('Solver MEX up to date (mex_%s)', modelName));
    fprintf('      C++ source unchanged — skipped.\n');
end

if nargout > 1
    iterCode     = fileread(matFiles.iterFile);
    simulateCode = fileread(matFiles.simulateFile);
    compileCode  = fileread(files.compileFile);
    gen = struct('iterCode', iterCode, 'simulateCode', simulateCode, ...
        'cppCode', cppText, 'compileCode', compileCode, 'cppCache', cppCache);
    gen.codeSegment = struct('iterCode', iterCode, 'simulateCode', simulateCode, ...
        'cppCode', cppText, 'compileCode', compileCode);
end
end

function validateOptions(options)
if isempty(options); return; end
if ~isstruct(options)
    error('gdsge:codegen:optionsNotAStruct', 'options must be a struct.');
end
fn = fieldnames(options);
% GenCodeSegment is a recognized legacy option with no counterpart in the
% refactored architecture; give it a specific diagnostic before the generic
% unknown-field sweep.
if any(strcmp(fn, 'GenCodeSegment'))
    error('gdsge:codegen:genCodeSegmentUnsupported', ...
        ['options.GenCodeSegment is not supported: the refactored thin-file ' ...
         'architecture (explicit named locals over gdsge.runtime helpers, no ' ...
         'v2struct) has no equivalent segment decomposition. Inspect the ' ...
         'generated iter_<model>.m and mex_<model>.cpp directly.']);
end
if ~isempty(fn)
    error('gdsge:codegen:unknownOption', ...
        'Unknown option field(s): %s. gdsge_codegen accepts no options yet.', ...
        strjoin(fn', ', '));
end
end

function banner(idx, total, label)
% One numbered phase header line for the codegen driver's progress output.
fprintf('[%d/%d] %s\n', idx, total, label);
end


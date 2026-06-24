function files = generateCxx(ir, outDir, opts)
% GENERATECXX  IR -> mex_<model>.cpp / compile_<model>.m written to outDir.
%   Adept-autodiff backend; spline or ASG interpolation (pchip: Phase 7c).
%   ASG models ensure the asg_mex kernel is built (emitCompile reads its
%   constants).  opts.includeDir overrides the default <repoRoot>/include.
if nargin < 2; outDir = pwd; end
if nargin < 3; opts = struct(); end
if ~exist(outDir, 'dir'); mkdir(outDir); end
gdsge.codegen.assertSupportedIR(ir);

% --- refuse what this phase does not support (clear errors, no wrong code)
if strcmp(ir.options.interpMethod, 'pchip')
    error('gdsge:codegen:unsupported', 'interpMethod pchip: Phase 7c');
end
gdsge.codegen.ensureKernels(ir);   % asg_mex, or spline constructor + evaluator
if ~isempty(ir.hooks.cxx)
    % preIter/postIter are MATLAB-side hooks (no C++ impact); preModel and
    % pre/postJac flow into the templates verbatim. Only cxx blocks are out.
    error('gdsge:codegen:unsupported', 'cxx hook blocks: Phase 7c');
end
% Phase 9a: MATLAB-side var_tensor (bounds/initial) never enters the C++ body,
% so no C++ guard is needed; assertSupportedIR rejects the ASG+tensor combo and
% analyzeModel rejects a tensor used in the model body.

% --- select the C++ model backend (driver may force via opts.backend; else the
%     IR jacobianBackend field decides, defaulting to autodiff) ---
if isfield(opts, 'backend')
    useSympy = strcmp(opts.backend, 'sympy');
else
    useSympy = isfield(ir.options, 'jacobianBackend') ...
        && strcmp(ir.options.jacobianBackend, 'sympy');
end
if useSympy
    if ~ismember(ir.options.interpMethod, {'spline', 'asg'})
        error('gdsge:codegen:sympyInterpUnsupported', ...
            'sympy backend supports spline or asg interpolation (got %s)', ir.options.interpMethod);
    end
    if ~isempty(ir.hooks.cxx)
        error('gdsge:codegen:cxxUnderSympy', ...
            'sympy backend cannot differentiate cxx blocks');
    end
    if ~gdsge.codegen.sympy.ensurePyenv()
        error('gdsge:codegen:sympyPythonUnavailable', ...
            ['sympy backend requires the uv Python env. Run: uv sync --project pyext']);
    end
    taskFn = @gdsge.codegen.cxx.sympymodel.emitTask;
else
    taskFn = @gdsge.codegen.cxx.emitTask;
end

if isfield(opts, 'includeDir')
    includeDir = opts.includeDir;
else
    here = fileparts(mfilename('fullpath'));       % src/+gdsge/+codegen
    repoRoot = fileparts(fileparts(fileparts(here)));
    includeDir = fullfile(repoRoot, 'include');
end

import gdsge.codegen.cxx.fillTemplate
import gdsge.codegen.cxx.readTemplate
taskInitCode = '';
if isfield(ir, 'modelInit')
    taskInitCode = taskFn(gdsge.codegen.initView(ir), 'task_init');
end
cpp = fillTemplate(readTemplate('mex.tpl.cpp'), { ...
    'GDSGE_OTHER_INCLUDE',    cxxIncludeText(ir); ...
    'TASK_INIT_CODE',         taskInitCode; ...
    'TASK_INF_HORIZON_CODE',  taskFn(ir)});

files = struct();
files.cppFile    = fullfile(outDir, ['mex_' ir.modelName '.cpp']);
gdsge.codegen.writeText(files.cppFile, cpp);

% Whole-loop SIMU_INTERP simulate MEX (spline + MEX-expressible simulate block).
files.simuMexFile = '';
if gdsge.codegen.cxx.isSimuMexExpressible(ir)
    files.simuMexFile = fullfile(outDir, ['simulate_' ir.modelName '_mex.cpp']);
    gdsge.codegen.writeText(files.simuMexFile, gdsge.codegen.cxx.emitSimulateMex(ir));
end

files.compileFile = fullfile(outDir, ['compile_' ir.modelName '.m']);
gdsge.codegen.writeText(files.compileFile, ...
    gdsge.codegen.cxx.emitCompile(ir, includeDir, ~isempty(files.simuMexFile)));
end

function s = cxxIncludeText(ir)
% Join hooks.cxxIncludes into a newline-separated block (empty when absent).
s = '';
if isfield(ir.hooks, 'cxxIncludes') && ~isempty(ir.hooks.cxxIncludes)
    s = strjoin(ir.hooks.cxxIncludes, newline);
end
end


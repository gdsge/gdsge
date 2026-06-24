function ir = parseFrontEnd(gmodText, modelName, baseDir)
% PARSEFRONTEND  gmod text + model name -> validated, semantically checked IR.
%   Optional baseDir resolves include()/cinclude() file names (default pwd).
if nargin < 3 || isempty(baseDir); baseDir = pwd; end
macro = gdsge.parser.expandMacros(gmodText, baseDir);
clean = gdsge.parser.preprocess(macro.text);
sb    = gdsge.parser.splitBlocks(clean);

decl    = gdsge.parser.parseDeclarations(sb.declText);
options = gdsge.parser.resolveOptions(decl.ws);

if ~isfield(sb.blocks, 'simulate')
    error('gdsge:parser:missingSimulate', 'No simulate block found');
end
sim = gdsge.parser.parseSimulate(sb.blocks.simulate);
decl.variables = gdsge.parser.resolveOutputs(decl.variables, sim.varSimu, decl.states.names);

if ~isfield(sb, 'modelRegions') || isempty(sb.modelRegions)
    error('gdsge:parser:missingModel', 'No model block found');
end
regions = cell(1, numel(sb.modelRegions));
for k = 1:numel(sb.modelRegions)
    body = gdsge.parser.parseModel(sb.modelRegions{k}.body);
    regions{k} = struct('condition', sb.modelRegions{k}.condition, ...
        'statements', {body.statements}, 'equations', {body.equations});
end
model = struct('regions', {regions});

modelInit = [];
hasInitVars = ~isempty(decl.initVariables.policyInit);
if isfield(sb.blocks, 'model_init')
    if ~hasInitVars
        error('gdsge:parser:modelInitNoVars', ...
            'model_init block present but no var_policy_init declared');
    end
    initBody = gdsge.parser.parseModel(sb.blocks.model_init);
    modelInit = struct( ...
        'variables',  decl.initVariables, ...
        'bounds',     {decl.boundsInit}, ...
        'statements', {initBody.statements}, ...
        'equations',  {initBody.equations});
elseif hasInitVars
    error('gdsge:parser:initVarsNoBlock', ...
        'var_policy_init declared but no model_init block found');
end

hooks = hooksFromBlocks(sb.blocks);
if ~isempty(macro.cxxIncludes)
    hooks.cxxIncludes = macro.cxxIncludes;
end
ir = gdsge.parser.assemblePartialIR(modelName, decl, sim, options, hooks, model, modelInit);
gdsge.parser.analyzeModel(ir);
end

function h = hooksFromBlocks(blocks)
h = struct('preModel','','preIter','','postIter','', ...
           'preJacCode','','postJacCode','','cxx','');
map = struct('pre_model','preModel', 'pre_iter','preIter', 'post_iter','postIter', ...
             'pre_jac_code','preJacCode', 'post_jac_code','postJacCode');
keys = fieldnames(map);
for i = 1:numel(keys)
    if isfield(blocks, keys{i})
        h.(map.(keys{i})) = blocks.(keys{i});
    end
end
end

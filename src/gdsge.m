function eq = gdsge(modelName)
% GDSGE  Backward-compatible flat orchestrator (frozen public API).
%   Generates code (gdsge_codegen), solves policies (iter_<model>) with an
%   IterRslt cache, simulates (simulate_<model>), and returns eq with fields
%   model, IterRslt, SmltRslt. The solve is skipped when the generated iter and
%   C++ code are unchanged and IterRslt_<model>.mat already exists. No v2struct.
fprintf('GDSGE: A Toolbox for Solving Global DSGE Models:\n');

[model, iterCode, cppCache, cppCode] = gdsge_codegen(modelName);

iterCacheName = ['iter_' modelName '.cache'];
iterCache = '';
if exist(iterCacheName, 'file') ~= 0
    iterCache = fileread(iterCacheName);
end

iterRsltName = ['IterRslt_' modelName '.mat'];
if strcmp(iterCache, iterCode) && strcmp(cppCache, cppCode) && exist(iterRsltName, 'file') ~= 0
    fprintf('Iter file not changed. IterRslt found. Load results:\n');
    loaded = load(iterRsltName, 'IterRslt');
    IterRslt = loaded.IterRslt;
    fprintf('ok\n');
else
    fprintf('Solve policies:\n');
    IterRslt = feval(['iter_' modelName]);
    save(iterRsltName, 'IterRslt');
    fid = fopen(iterCacheName, 'w');
    if fid == -1
        warning('gdsge:orchestrator:cacheWriteFailed', ...
            'Could not write iter cache file: %s', iterCacheName);
    else
        fprintf(fid, '%s', iterCode);
        fclose(fid);
    end
    fprintf('ok\n');
end

fprintf('Simulate:\n');
SmltRslt = feval(['simulate_' modelName], IterRslt);
fprintf('ok\n');

eq.model    = model;
eq.IterRslt = IterRslt;
eq.SmltRslt = SmltRslt;
end

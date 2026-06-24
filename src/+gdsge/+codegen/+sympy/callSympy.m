function out = callSympy(fnName, requestStruct)
% CALLSYMPY  The single MATLAB->SymPy chokepoint. jsonencode the request,
%   invoke py.gdsge_sympy.<fnName>(jsonChar), jsondecode the returned JSON.
%   Assumes gdsge.codegen.sympy.ensurePyenv() has already returned true.
reqJson = jsonencode(requestStruct);
pyFn = py.getattr(py.importlib.import_module('gdsge_sympy'), fnName);
resJson = char(pyFn(reqJson));
out = jsondecode(resJson);
end

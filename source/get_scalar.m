function outStruct = get_scalar(inStruct,state_vector)
if ~(isequal(state_vector{1},1) && length(state_vector)==2 && max(state_vector{2}==0)==1)
    outStruct = inStruct;
    return;
end
% A helper function to get the first element of the structure
varnameList = fieldnames(inStruct);
outStruct = struct;
for i=1:length(varnameList)
    varname = varnameList{i};
    if ~isempty(inStruct.(varname))
        outStruct.(varname) = inStruct.(varname)(1);
    else
        outStruct.(varname) = [];
    end
end
end
function InterpEval = convert_to_interp_eval_array(ppVec)
% CONVERT_TO_INTERP_EVAL_ARRAY(pp) converts MKL pp coefficients to interp eval
% Convert from a cell of pps, extract their common grid
% information and store coefficient in order (i_vec,coefs,i_array)
% @author: Wenlan Luo (luowenlan@gmail.com)

if isempty(ppVec)
    InterpEval = struct;
    return;
end

% Check validity of the input
if iscell(ppVec)
    sizeCoefsCommon = size(ppVec{1}.coefs);
    numVecInCell = length(ppVec);
    for i=1:numVecInCell
        if ~isequal(sizeCoefsCommon, size(ppVec{i}.coefs))
            error('coefs size is not consistent');
        end
    end
    
    % Stack coefficients
    numArray = ppVec{1}.dim;
    coefs = zeros([prod(sizeCoefsCommon),numVecInCell]);
    for i_vec=1:numVecInCell
        coefs(:,i_vec) = ppVec{i_vec}.coefs(:);
    end
    % Now coefs is in order (coefs,i_array,i_vec)
    
    % Construct dimension and permutation order of coefs
    pp = ppVec{1};
else
    % Treat it as a large vector function and num array to 1
    numArray = 1;
    sizeCoefsCommon = size(ppVec.coefs);
    numVecInCell = sizeCoefsCommon(end);
    coefs = ppVec.coefs;
    pp = ppVec;
end
dimCoefs = [];
xDim = length(pp.pieces);
for i=1:xDim
    dimCoefs = [dimCoefs,pp.order(xDim+1-i),pp.pieces(xDim+1-i)];
end

dimCoefs = [dimCoefs,numArray,numVecInCell];
coefs = reshape(coefs,dimCoefs);

% Stack coefficients in (i_vec,coefs,i_array)
newPermuteOrder = [length(dimCoefs),length(dimCoefs)-3:-2:1,length(dimCoefs)-2:-2:2,length(dimCoefs)-1];
coefs = permute(coefs,newPermuteOrder);

newCoefsSize = size(coefs);

% Now coefs is in the order
% [i_vec,order_1:-1:1,order_2:-1:1, ... , pieces_1, piece_2, ..., i_array]

% Attach other properties
dim = numVecInCell;
fullVecEvalCoefsLength = double(prod(newCoefsSize(1:xDim+1)));
singleVecEvalCoefsLength = double(fullVecEvalCoefsLength / dim);
order = pp.order;
pieces = pp.pieces;
xPts = pp.pieces+1;
breaks = pp.breaks;

arrayOffset = double(prod(newCoefsSize(1:(1+2*xDim))));

InterpEval = v2struct(coefs,fullVecEvalCoefsLength,singleVecEvalCoefsLength, ...
    xDim,order,pieces,xPts,breaks,dim,arrayOffset);

end

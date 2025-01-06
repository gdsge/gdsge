function [v,index_x] = myppual(pp,x,left,indexDim,index_x)
%MYPPUAL Construct and Evaluate splines in ppform at flexible vector-valued dimensions, table look-up index, and spline dimension reduction
% in both vectorized pure MATLAB code and CMEX implementation
%
%MYPPUAL works like the Curve Fitting Toolbox function "ppual," but improves it in several ways:
% First, "flexible vector-valued dimension selection": The evaluation can be done for any vector-valued dimension by specifying indexDim
% Second, "flexible table look-up index": the evaluation can take precomputed table lookup index of x by specifying index_x
% Third, "flexible variables dimension selection," i.e., "spline dimension reduction": partial evaluation at scalar points   
%        by specifying a cell or double vector x as an upper dimensional subset of the full variables, in which case the resulting  
%        output is a lower-dimensional pp spline
% Fourth, CMEX implementation of selected features mentioned above 
%
%Input:
%   pp: a structure of spline information, in the format of the following: 
%       MATLAB Curve Fitting Toolbox ppform; MATLAB griddedInterpolant; or myppual MKL spline forms 
%       (a) MATLAB pp form: constructed by base MATLAB (e.g., spline, pchip) or Curve Fitting Toolbox spline functions (e.g., csape)
%           as a structure from struct('form','pp','breaks',breaks,'coefs',coefs,'pieces',pieces,'order',order,'dim',dim)
%       (b) MKL pp form:
%           with all fields from MATLAB pp form, plus those from struct('Values',Values,'Method',Method), i.e.,
%           pp has the following fields, 
%           struct('form','MKLpp','breaks',break,'Values',Values,'Method',Method,'coefs',coefs,'pieces',pieces,'order',order,'dim',dim),
%           where 'forms', 'breaks',''Values', and 'Method' are required fields during the construction stage. Other fields are returned
%           values from the constructor.
%           breaks: a 1 by m cell vector with m grid vectors
%           Values: specifies a (m+1)-dimensional array, where the first dimension is the vector-valued dimension 
%           Method: a 1 by m cell vector specifying the end point conditions for splines (e.g., {'not-a-knot'} for cubic splines)
%           coefs: the coefficient array as ordered in MKL format. If coefs is empty, we call MKL constructor to construct the coefficient
%       (c) MATLAB griddedInterpolant: documentation to be added
%    x: data sites to be evaluated. x can take three forms:  
%       (a) a matrix: this form follows the same syntax as "ppual"
%       (b) a cell vector with isequal(length(x),length(pp.breaks)): this form follows the same syntax as "ppual"
%       (c) a cell vector with "length(x) < length(pp.breaks)" and "isscalar(x(i)) for all i": 
%           in this case we perform "spline dimension reduction," i.e., evaluate at length(x) outer dimension,
%           and return the resulting lower dimension spline as ppNew
%           
% left: (default: 'l') a string to specify left continuity, as standard in Curve Fitting Toolbox
%       (a detailed documentation can be found in "doc ppual")
% indexDim: (default: empty) a double matrix or array to specify the indices of dimensions to
%        be evaluated, where size(indexDim,1) indicates how many dimensions to be evaluated.
%        IndexDim can take one of three forms:
%       (a) if x is a matrix, then indexDim(:,j) represents the indices of dimensions to be 
%           evaluated at data site j. Hence, size(indexDim,2) has to equal to size(x,2).
%       (b) if x is a cell array, then 
%       (b.1)if indexDim is a multi-dimensional array, then indexDim(:,i,j,...,end)
%            specifies the dimensions to pick out corresponding to x{i},x(j),...,x(end)
%            in cell x.
%       (b.2)if indexDim is a matrix, it corresponds to indexDim in (b.1) as  
%            indexDim = reshape(indexDim,size(indexDim,1),[])
%       (c) if indexDim is empty matrix, then myppual evaluates at full dimension, i.e., the same as "ppual" 
% index_x: (default: empty) a double matrix or cell array to specify the precomputed indices of table look-up result
%
%Output: 
%    v: a double matrix or multi-dimensional array or pp spline in the following forms:
%   (a) for spline evaluation case: the same form as the output of ppual, with the exception that size(v,1) = size(indexDim,1)
%   (b) MATLAB pp form in the case of "spline dimension reduction"
%   (c) MKL pp form in the case of "spline construction" or "spline dimension reduction"
%
%   See also PPUAL, PPVAL, SPLINE, PCHIP, INTERPN, GRIDDEDINTERPOLANT
%{ 
    Authors: Jinhui Bai (Version 1.0 - current), jinhui.bai@gmail.com, Department of Economics, Georgetown University, USA 
             Wenlan Luo (Version 3.0 - current), wl239@georgetown.edu, Department of Economics, Georgetown University, USA

    TODO List for the version 3.0: C++ implementation of selected features of the current MATLAB myppual in a CMEX function.
    Essentially myppual becomes a gateway to low-level spline routine. 
    In particular, we are going to implement the following features gradually:
    (a)Interface Intel MKL C++ 1-D Spline construction, vector evaluation of values and derivatives, in a CMEX file
        Provide a new spline type "MKLpp" as a structure, so that all MKL spline construction (e.g., quadratic spline) or derivative evaluation 
        can be interfaced to MATLAB through myppual.
    (b)Interface C++ tensor product extension to Intel MKL (as of now, MKL only includes 1-D spline, which is quite limited for macro research)
    (c) Most importantly (and the unique feature of our CMEX file), implement C++ translation of MYPPUAL in a CMEX file 
       (using pointer switching to avoid memory copy; using hunt as well as bisection for table lookup; an option for OPEN MP Multithreading)
        
    Version History:
    Version 3.0 in progress: last update, 05/09/2014
    version 2.0: 11/26/2013. Major revision: introduces "spline dimension reduction" feature both for cell form and matrix form
    version 1.1: 11/23/2011. Minor revision: transform formats between base MATLAB splines and Curve Fitting Toolbox
    version 1.0: 10/15/2008. Major release with "flexible vector dimension selection" and "table look-up index" features   
%}
%{
  Input Data Check
%}
% Number of input variables and specify default values
if nargin < 5 % give default values to index_x: an empty matrix implies table lookup inside the code
    index_x = [];
    if nargin < 4 % give default values to indexDim: an empty matrix is treated as evaluating at all vector dimensions
        indexDim = [];
        if nargin < 3 % give default values: an empty value is treated as left continuity
            left = [];
            if nargin < 2
                x = [];
                if nargin < 1
                    error('SPLINES:MYPPUAL:needinput','At least one input variable is required');
                end
            end
        end
    end
end
% check the spline type: currently support pp and prp form, griddedInterpolant, and Intel MKL C Library
if isa(pp,'struct') % a structure of pp form
    % check the category of pp form    
    check_field = all(isfield(pp,{'form','breaks','coefs','order'}));
    if ~check_field,error('SPLINES:MYPPUAL:input','Missing fields in pp input'),end    
    form = pp.form;
    coefs = pp.coefs;    
    if strncmp(form,'pp',2) % MATLAB pp spline
        ispp_MATLAB = true;
        isprp_MATLAB = false;
        ispp_MKL = false;
        isobj_griddedInterpolant = false;
    elseif strncmp(form,'prp',2) % MATLAB pp spline
        ispp_MATLAB = false;
        isprp_MATLAB = true;
        ispp_MKL = false;
        isobj_griddedInterpolant = false;        
    elseif strncmp(form,'MKLpp',2) % Intel MKL pp spline
        ispp_MATLAB = false;
        isprp_MATLAB = false;        
        ispp_MKL = true;
        isobj_griddedInterpolant = false;
    elseif strncmp(form,'griddedInterpolant',2) % MATLAB griddedInterpolant Structure Interface
        ispp_MATLAB = false;
        isprp_MATLAB = false;        
        ispp_MKL = false;
        isobj_griddedInterpolant = true;
    end
    % check the variable dimensionality of the pp-form spline
    if iscell(pp.breaks) && (length(pp.breaks) > 1) % multi-variate spline
        isMultiVariateSpline = true;
    else % univariate spline: either pp.break is a vector, or pp.break is a one by one cell
        isMultiVariateSpline = false;
    end    
elseif isa(pp,'griddedInterpolant') % MATLAB griddedInterpolant object
    form = 'griddedInterpolant';
    ispp_MATLAB = false;
    isprp_MATLAB = false;
    ispp_MKL = false;
    isobj_griddedInterpolant = true;
    myppual_flag = 'eval';
    % check the variable dimensionality of the pp-form spline    
    if iscell(pp.GridVectors) && (length(pp.GridVectors) > 1) % multi-variate spline
        isMultiVariateSpline = true;
    else % univariate spline: either pp.break is a vector, or pp.break is a one by one cell
        isMultiVariateSpline = false;
    end        
else % input error
    error('SPLINES:MYPPUAL:input','The first input has to be either a structure or griddedInterpolant');
end

% check the purpose of the function calling: if it is transformation, do it and return
if isequal(nargin,1) && isa(pp,'struct')
    if isempty(coefs) % full spline construction
        myppual_flag = 'cons';
    else % spline conversion from pp form to MKLpp form
        if ispp_MATLAB
            v = pptransformation(pp,'pp2MKLpp'); return;
        elseif isprp_MATLAB
            pp.form = 'pp'; % change dim to pp-form equivalent dim
            pp.dim = 3*pp.dim; % change dim to pp-form equivalent dim
            v = pptransformation(pp,'pp2MKLpp'); return;
        elseif ispp_MKL
            v = pptransformation(pp,'MKLpp2pp'); return;
        end       
    end
end

% For MATLAB pp form, transform indexDim to make it a double precision matrix
% Reason: if indexDim is integer type, it will cause problem in calculation
if ispp_MATLAB || isprp_MATLAB
    if (~isempty(indexDim)) && (~isa(indexDim,'double'))
        indexDim = double(indexDim);
    end
end

%{
  Now the real job of pp evaluation
%}
if ispp_MATLAB && isMultiVariateSpline
    %---------------------------------------------------------
    % Case I: evaluate multivariate MATLAB pp-form spline with m variables
    %---------------------------------------------------------    
    breaks = pp.breaks; order = pp.order;
    % check pieces
    if isfield(pp,'pieces') && (~isempty(pp.pieces))
       pieces = pp.pieces; 
    else % calculate it from breaks
        pieces = cellfun('prodofsize',pp.breaks) - 1;        
    end
    % check dim
    if isfield(pp,'dim') && (~isempty(pp.dim))
       dim = pp.dim;
    else % calculate it from breaks
       dim = size(pp.Values,1);       
    end
    % give default values to input variables
    m = length(breaks); % the total number of variables in the spline
    % check the number of variables in the inquiry point
    if iscell(x) % cell form
        mx = length(x);
    else % matrix form
        mx = size(x,1);
    end
   % set up LEFT appropriately 
   if ~isempty(left) 
      if  ~iscell(left)
         temp = left; left = cell(1,mx); left(:) = {temp};
      end
   else % empty one: change it to the cell format
      left = cell(1,mx);
   end
   % set up indexDim appropriately
   if isempty(indexDim)
       mrowIndexDim = dim;       
   else
       mrowIndexDim = size(indexDim,1);
   end
   % set up index_x appropriately
   if ~isempty(index_x) % non-empty one: change it to the corresponding cell or matrix format
      if  iscell(x) 
          if ~iscell(index_x) % for reasons of backward compatibility with Version 1.0: we let index_x to be the first element of the cell
              temp = index_x; index_x = cell(1,mx); index_x{1} = temp;
          end
      end
   else % empty one: change it to either the cell or double form of empty variable 
       if iscell(x) % cell form of x: give default value as empty cell with the same dimension as x
           index_x = cell(1,mx);
       else % matrix form of x: give default value as empty matrix
           index_x = [];
       end
   end   
   %----------------------------------------------------------
   % The evaluation step for m-variate spline: case-specific treatment for eight different cases
   % (1) full tensor-product spline construction
   % (2) m-variate matrix form: combined construction and evaluation
   % (3) m-variate matrix from: general evaluation   
   % (4) m-variate cell form: combined construction and evaluation
   % (5) m-variate cell form with scalar evaluation point for variables 2:end
   % (6) m-variate cell form: general evaluation            
   % (7) mx(< m)-variate cell form with scalar point for spline dimension reduction    
   % (8) mx(< m)-variate matrix form with scalar point for spline dimension reduction    
   %----------------------------------------------------------
   % set up common tasks for spline construction   
   % check and isolate cases
   if (nargin > 1) && iscell(x) % cell form
       if isempty(coefs) && isequal(mx,m) % spline construction and evaluation in one
           myppual_flag = 'cons_eval_cell';
       elseif (~isempty(coefs)) && (mx < m) % spline reduction case
           myppual_flag = 'eval_cell_part';
           % check that all points are scalars
           is_allscalar = all(cellfun(@isscalar,x,'UniformOutput',true));
           if ~is_allscalar
               error('SPLINES:MYPPUAL:cellform','Spline reduction requires scalar inquiry points');
           end
       elseif (~isempty(coefs)) && isequal(mx,m) % normal cell-form spline evaluation case
           is_outerscalar = all(cellfun(@isscalar,x(2:end),'UniformOutput',true));  % the scalar variables for outer dimension
           if (m > 2) && is_outerscalar % with scalar point at the outer dimensions
               myppual_flag = 'eval_cell_outerscalar';
%                myppual_flag = 'eval_cell_general';
           else % general spline evaluation case
               myppual_flag = 'eval_cell_general';
           end
       else
           error('SPLINES:MYPPUAL:cellform','Please check the dimension of inquiry points');
       end       
   elseif (nargin > 1) && isnumeric(x) && ismatrix(x) % matrix form
       if isempty(coefs) && isequal(mx,m) % spline construction and evaluation in one
           myppual_flag = 'cons_eval_matrix';           
       elseif (~isempty(coefs)) && (mx < m) % spline reduction case
           myppual_flag = 'eval_matrix_part';
           % check that x is a column vector
           if ~iscolumn(x)
               error('SPLINES:MYPPUAL:matrixform','Spline reduction requires a column vector');
           end
       elseif (~isempty(coefs)) && isequal(mx,m) % normal matrix-form spline evaluation case
           myppual_flag = 'eval_matrix_general';
       else
           error('SPLINES:MYPPUAL:matrixform', ...
               ['Each X(:,j) must not have more elements than a ',num2str(m),'-vector.']);
       end
   end
   %----------------------------------------------------------
   % Start the case-specific evaluation 
   %----------------------------------------------------------   
   if strcmp(myppual_flag,'eval_cell_outerscalar')
      %---------------------------------------------------
      % cell form: scalar points for outer dimensions
      %---------------------------------------------------
      v = coefs; % initialize the value as coefficients
      sizev = [dim,pieces.*order]; % size of the coefficients
      % change outer variables to matrix form 
      dd = prod(sizev(1:2)); % the combined dimension from vector dim and the first dimension
      x_cell_outerscalar = x(2:end); % select the outer dimension from x cell
      x_matrix_outerscalar = reshape([x_cell_outerscalar{:}],[],1); % change scalar points from cell to matrix form
      index_x_cell_outerscalar = index_x(2:end); % select the outer dimension from index_x_cell
      index_x_matrix_outerscalar = reshape([index_x_cell_outerscalar{:}],[],1); % change index for scalar points from cell to matrix form
      % use the matrix form to evaluate the outer scalar point
      v_outerscalar = myppual_matrix_general(breaks(2:m),reshape(v,[dd,sizev(3:end)]),pieces(2:m),order(2:m),dd,...
                                    x_matrix_outerscalar,left(2:m),[],index_x_matrix_outerscalar);
      % use the 1-d routine to evaluate the first x variable                          
      v = myppual1(breaks{1},reshape(v_outerscalar,[dim*pieces(1),order(1)]),pieces(1),order(1),dim,...
                                    x{1},left{1},indexDim,index_x{1});
      % reshape to the required output format
      v = reshape(v,mrowIndexDim,[]);
   elseif strcmp(myppual_flag,'eval_cell_general')
      %---------------------------------------------------
      % cell form: general case for arbitrary combinations of cell
      %---------------------------------------------------
      v = myppual_cell_general(breaks,coefs,pieces,order,dim,x,left,indexDim,index_x);
   elseif strcmp(myppual_flag,'eval_matrix_general')
      %---------------------------------------------------
      % matrix case: evaluation at scattered points
      %---------------------------------------------------
      v = myppual_matrix_general(breaks,coefs,pieces,order,dim,x,left,indexDim,index_x);      
   elseif strncmp(myppual_flag,'cons',4)
      %---------------------------------------------------
      % Three Cases involing spline construction
      %---------------------------------------------------
      Values = pp.Values; % initialize the value
      % check the pp.Method field
      if isfield(pp,'Method') && (~isempty(pp.Method))
          Method = pp.Method;
          if iscellstr(Method) && isequal(length( Method),1) %
              Method = repmat( Method,1,m);
          elseif iscellstr( Method) && isequal(length(Method),m) %
              Method = pp.Method;
          elseif ischar(Method)  % a string
              Method = repmat({Method},1,m);
          else
              error('Incorrect Input Form pp.Method');
          end
      else % give it a default one as the most common usage
          Method = repmat({'not-a-knot'},1,m);
      end
      % Check the extrapolation
      if isfield(pp,'ExtrapolationOrder') && (~isempty(pp.ExtrapolationOrder))
          ExtrapolationOrder = pp.ExtrapolationOrder;
          if isscalar(ExtrapolationOrder) %
              ExtrapolationOrder = repmat(ExtrapolationOrder,1,m);
          elseif isvector(ExtrapolationOrder) && isequal(length(ExtrapolationOrder),m) %
              ExtrapolationOrder = pp.ExtrapolationOrder;
          else
              error('Incorrect Input Form pp.ExtrapolationOrder');
          end
      else % give it a default one as the most common usage
          ExtrapolationOrder = order;
      end
      % do the actual computation
      if strcmp(myppual_flag,'cons_eval_cell') % cell form combined construction and evaluation
          v = myppual_cell_cons_eval(breaks,Values,pieces,order,dim,x,left,indexDim,index_x,Method,ExtrapolationOrder);
      else % matrix
        % construction of the spline
        if isequal(order(:),4*ones(m,1))
            ppV = csape(breaks,Values,Method);
        else
            error('Currently we only support spline construction in cubic form');
        end
        % add extrapolation to the spline
        if ~isequal(ExtrapolationOrder,order),ppV = fnxtr(ppV,ExtrapolationOrder);end
        % now the following parts may change after fnxtr, so we have to redefine them.
        breaks = ppV.breaks; coefs = ppV.coefs; pieces = ppV.pieces;
        if strcmp(myppual_flag,'cons') % full construction
            v = ppV;
        elseif strcmp(myppual_flag,'cons_eval_matrix')
            % call matrix evaluation routine
            v = myppual_matrix_general(breaks,coefs,pieces,order,dim,x,left,indexDim,index_x);
        end
      end
   elseif strcmp(myppual_flag,'eval_cell_part')
      %---------------------------------------------------
      % cell form: pp spline reduction
      %---------------------------------------------------
      v = coefs; % initialize the value as coefficients
      sizev = [dim,pieces.*order]; % size of the coefficients
      outer_start = m - mx + 1; % the index for the first outer variable
      % change outer variables to matrix form
      dd = prod(sizev(1:outer_start)); % the combined dimension from vector dim and the inner variables
      % use the matrix form to evaluate the outer scalar point
      v_outerscalar = myppual_cell_general(breaks(outer_start:m),reshape(v,[dd,sizev((outer_start+1):end)]),...
                                    pieces(outer_start:m),order(outer_start:m),dd,...
                                    x,left(1:mx),[],index_x);
      % assemble the spline with the reduced dimensions: notice that we have to use another layer of cell to include the cell vector of breaks,
      % otherwise ppInner becomes a structure array (see "doc struct" for details)
      ppInner = struct('form','pp','breaks',{breaks(1:(outer_start-1))},'coefs',reshape(v_outerscalar,sizev(1:outer_start)),...
                       'pieces',pieces(1:(outer_start-1)),'order',order(1:(outer_start-1)),'dim',dim);
      v = ppInner; % return the final result as a pp spline
   elseif strcmp(myppual_flag,'eval_matrix_part')
       %---------------------------------------------------
       % matrix form: pp spline reduction
       %---------------------------------------------------
      v = coefs; % initialize the value as coefficients
      sizev = [dim,pieces.*order]; % size of the coefficients
      outer_start = m - mx + 1; % the index for the first outer variable
      % change outer variables to matrix form
      dd = prod(sizev(1:outer_start)); % the combined dimension from vector dim and the inner variables
      x_matrix_outerscalar = x; % change scalar points from cell to matrix form
      index_x_matrix_outerscalar = index_x; % change index for scalar points from cell to matrix form
      % use the matrix form to evaluate the outer scalar point
      if isequal(mx,1)
          v_outerscalar = myppual1(breaks{m},reshape(v,[dd*pieces(m),order(m)]),...
              pieces(m),order(m),dd,x_matrix_outerscalar,left(1:mx),[],index_x_matrix_outerscalar);                        
      else
          v_outerscalar = myppual_matrix_general(breaks(outer_start:m),reshape(v,[dd,sizev((outer_start+1):end)]),...
              pieces(outer_start:m),order(outer_start:m),dd,x_matrix_outerscalar,left(1:mx),[],index_x_matrix_outerscalar);              
      end
      % assemble the spline with the reduced dimensions: notice that we have to use another layer of cell to include the cell vector of breaks,
      % otherwise ppInner becomes a structure array (see "doc struct" for details)
      ppInner = struct('form','pp','breaks',{breaks(1:(outer_start-1))},'coefs',reshape(v_outerscalar,sizev(1:outer_start)),...
                       'pieces',pieces(1:(outer_start-1)),'order',order(1:(outer_start-1)),'dim',dim);
      v = ppInner; % return the final result as a pp spline       
   end   
elseif (ispp_MATLAB || isprp_MATLAB) && (~isMultiVariateSpline)
    %---------------------------------------------------------
    % Case II: evaluate univariate MATLAB pp-form or prp-form spline
    %---------------------------------------------------------
    %{
    if pp.breaks is a 1 by 1 cell, then we need to deal with two problems:
    (a)change data type of pp.break and x from cell to double vector
    (b)change the organization of pp.coefs. We have to do this because MATLAB has two conflicting ways of 
       coefficient organization in MATLAB base toolbox and Curve Fitting Toolbox.
        By default, in cell case pp.coefs is organized in "Curve Fitting Toolbox spline tradition,"  
        with [pp.dim,pp.pieces*pp.order] = size(pp.coefs).
        After reshaping it as (pp.dim*pp.pieces) by pp.order, we transform it
        to MATLAB Base Toolbox 1-dim spline storage tradition so that myppual1 works for both traditions.
    %}
    % extract common components. Do NOT extract case-specific components since they may not exist.
    breaks = pp.breaks; coefs = pp.coefs; order = pp.order;
    if isfield(pp,'Method') && (~isempty(pp.Method))
        Method = pp.Method;
        if iscell(Method),Method = Method{1};end
    else % give it a default one as the most common usage
        Method = 'not-a-knot';
    end
    if isempty(coefs) % spline construction
        % change break from cell to vector: otherwise fnxtr will throw out an error
        if iscell(breaks),breaks = breaks{1};end
        % check method field        
        % form the spline: currently only support linear or cubic spline. Will be updated in the future subject to needs
        if ispp_MATLAB
            if isequal(order,4)
                if strcmp(Method,'pchip')
                    ppV = pchip(breaks,pp.Values);
                elseif strcmp(Method,'pwch')
                    ppV = pwch(breaks,pp.Values,pp.Der1);
                else
                    ppV = csape(breaks,pp.Values,Method);
                end
            elseif isequal(order,2)
                ppV = pplin1(breaks,pp.Values);
            elseif isequal(order,3)
                if strcmp(Method,'rs3Hermite')
                    ppV = rs3Hermite(breaks,pp.Values,pp.Der1);
                end
            else
                error('The proposed spline form is not supported currently');
            end
            % deal with extrapolation
            if isfield(pp,'ExtrapolationOrder') && (~isempty(pp.ExtrapolationOrder))
                ExtrapolationOrder = pp.ExtrapolationOrder;
                if (ExtrapolationOrder < order),ppV = fnxtr(ppV,ExtrapolationOrder);end
            end            
        elseif isprp_MATLAB
            if isequal(order,3)
                if strcmp(Method,'dgHermite')
                    ppV = rs3dgHermite(breaks,pp.Values,pp.Der1);
                end
            else
                error('The proposed spline form is not supported currently');
            end
        end
        % check whether we are doing construction only
        if isequal(nargin,1)
            v = ppV; return;
        end
        % now the following parts may change after fnxtr, so we have to redefine them.
        breaks = ppV.breaks; coefs = ppV.coefs; pieces = ppV.pieces; dim = prod(ppV.dim);
    else % preconstructed spline
        if iscell(breaks)
            breaks = breaks{1};
            coefs = reshape(coefs,[],order);
        end
        pieces = pp.pieces;
        dim = prod(pp.dim);  % dim may be a vector in MATLAB spline convention: transform to a scalar        
    end
    % deal with data sites
    if iscell(x)
        if isequal(length(x),1)
            x = x{1};
        else
            error('SPLINES:MYPPUAL:wrongx','X should specify a 1-dimensional vector.');                
        end        
    end
    % deal with index_x
    if iscell(index_x), index_x = index_x{1};end
    % call myppual1 to do the actual evaluation
    if ispp_MATLAB
        v = myppual1(breaks,coefs,pieces,order,dim,x,left,indexDim,index_x);
    elseif isprp_MATLAB
        % convert pp form splines indexing to rs form
        if ~isempty(indexDim)
            indexDimRS = [indexDim;indexDim + dim;indexDim + 2*dim];
        end
        % evalute the spline as if it is a pp form
        dimRS = 3*dim; % the pp-equivalent dim for rational spline                
        vTemp = myppual1(breaks,coefs,pieces,order,dimRS,x,left,indexDimRS,index_x);
        vTemp = reshape(vTemp,[],3,size(vTemp,2));
        vTemp1 = reshape(vTemp(:,1,:),[size(vTemp,1),size(vTemp,3)]);
        vTemp2 = reshape(vTemp(:,2,:),[size(vTemp,1),size(vTemp,3)]);
        vTemp3 = reshape(vTemp(:,3,:),[size(vTemp,1),size(vTemp,3)]);
        % calculate the final value
        v = vTemp1; % initialize v as the fall back value
        idx_none_zero = (vTemp3 ~= zeros(size(vTemp3)));
        v(idx_none_zero) = vTemp1(idx_none_zero) + vTemp2(idx_none_zero)./vTemp3(idx_none_zero);
    end    
elseif ispp_MKL
    %---------------------------------------------------------
    % Case III: construct and evaluate Intel MKL pp-form spline with m variables. This part is only partially enabled.
    %---------------------------------------------------------
    if iscell(pp.breaks)
        breaks = pp.breaks;
    else % make 1-D break a one by one cell
        breaks = {reshape(pp.breaks,1,[])};
    end 
    m = length(breaks);
    breaks_MKL = breaks;
    coefs_MKL  = pp.coefs;
    order_MKL = int32(pp.order);   % typecasting to int
    pieces_MKL = int32(cellfun('prodofsize',breaks_MKL) - 1); % typecasting to int
    if isa(indexDim,'int32') % for int32, we assume that the user supplied correct C-type data
        indexDim_MKL = indexDim; % make no adjustment
    else % for other data type, we have to make the adjustment to types and zero-based indexing
        indexDim_MKL = int32(indexDim - 1); % typecasting to int and to C zero-based indexing
    end
    if iscell(index_x) && isequal(length(index_x),1)
        if isa(index_x{1},'int32') % make no adjustment
            index_x_MKL = index_x{1};
        else % typecasting
            index_x_MKL = int32(index_x{1}-1);
        end
    elseif iscell(index_x) && (length(index_x) > 1)
        % currently, CMEX only uses the first dimension, so we only typecast the first one
        index_x_MKL = index_x;
        if isa(index_x{1},'int32') % make no adjustment
            index_x_MKL{1} = index_x{1};
        else % typecasting
            index_x_MKL{1} = int32(index_x{1}-1);
        end
    elseif isnumeric(index_x) && ismatrix(index_x)
        if isa(index_x,'int32') % make no adjustment
            index_x_MKL = index_x;
        else % typecasting
            index_x_MKL = int32(index_x-1);
        end        
    end
   % set up threads
   MKL_flag = int32(1); % default value as single thread
   if isfield(pp,'thread') && (~isempty(pp.thread)),MKL_flag = int32(max(pp.thread,1)); end % MKL_flag is the number of thread to generate
   % set up method
   Method = 'not-a-knot'; % default method
   if isfield(pp,'Method') && (~isempty(pp.Method)),Method = pp.Method;end

    % check the purpose of the operation: currently we have enabled four cases
    % (1) full tensor-product spline construction; 
    % (2) matrix-form combined construction and evaluation;
    % (3) matrix-form full evaluation
    % (4) cell-form combined construction and evaluation
    if (nargin > 1) && iscell(x) % cell form
        mx = length(x);
        if isequal(mx,1)
            is_outerscalar = true;
        else
            is_outerscalar = all(cellfun(@isscalar,x(2:end),'UniformOutput',true));  % the scalar variables for outer dimension
        end
       if isempty(coefs_MKL) && isequal(mx,m) && is_outerscalar % spline construction and evaluation in one
           myppual_flag = 'cons_eval_cell';
       else
           error('SPLINES:MYPPUAL:cellform','MKL Spline does not allow the proposed evaluation format');
       end
    elseif (nargin > 1) && ismatrix(x) % matrix form
        mx = size(x,1);
        if isempty(coefs_MKL) && isequal(mx,m) % spline construction and evaluation in one
            myppual_flag = 'cons_eval_matrix';
        elseif (~isempty(coefs)) && isequal(mx,m) % normal matrix-form spline evaluation case
            myppual_flag = 'eval_matrix_general';
        else
            error('SPLINES:MYPPUAL:cellform','MKL Spline does not allow the proposed evaluation format');
        end        
    end
   %----------------------------------------------------------
   % Start the case-specific evaluation 
   %----------------------------------------------------------
   if strncmp(myppual_flag,'cons',4)
       if isfield(pp,'orient') && isequal(pp.orient,'MKLC')
           Values_MKL = pp.Values;
           dim_MKL = int32(size(pp.Values,m+1)); % the dim is the last dimension for permuted value
       elseif (~isfield(pp,'orient')) || (isfield(pp,'orient') && isequal(pp.orient,'curvefit'))
           % matlab curve fitting order: permute the variables. Fortran vs C storage transformation, sigh...
           Values_MKL = permute(pp.Values,((m+1):-1:1));
           dim_MKL = int32(size(pp.Values,1));  % the dim is the first dimension for the original value
       end
       if strcmp(myppual_flag,'cons_eval_cell') % needs a partial construction and evaluation
           if isequal(mx,1) % combined full 1-d construction and evaluation
               % first construct the coefs
               MKL_flag_cons = -MKL_flag;
               coefs_MKL = myppual_mex(MKL_flag_cons,breaks_MKL,Values_MKL,pieces_MKL,order_MKL,dim_MKL,Method,reshape(x{1},1,[]),left,indexDim_MKL,index_x_MKL);
               [v,index_x] = myppual_mex(MKL_flag,breaks_MKL,coefs_MKL,pieces_MKL,order_MKL,dim_MKL,Method,reshape(x{1},1,[]),left,indexDim_MKL,index_x_MKL);
               index_x = {index_x};
           else  % combined partial construction and evaluation
               [v,index_x] = myppual_mex(MKL_flag,breaks_MKL,Values_MKL,pieces_MKL,order_MKL,dim_MKL,Method,x,left,indexDim_MKL,index_x_MKL);
           end
       else % need a full construction first           
           % first construct the coefs
           MKL_flag_cons = -MKL_flag;

           % Dealing with extrapolation
           extrap_order_MKL = [];
           if isfield(pp,'ExtrapolationOrder') && (~isempty(pp.ExtrapolationOrder))
               if ~all(pp.order == 4)
                    error('SPLINES:MYPPUAL:MKLform','Can only apply extrapolation to cubic splines');
               end
               ExtrapolationOrder = pp.ExtrapolationOrder;
               if isscalar(ExtrapolationOrder)
                   ExtrapolationOrder = repmat(ExtrapolationOrder,[1,length(pp.breaks)]);
               end
               extrap_order_MKL = int32(ExtrapolationOrder);

               % attach new breaks
               for i_break=1:length(breaks_MKL)
                   c_break = breaks_MKL{i_break};
                   c_break = [c_break(1)-1, c_break, c_break(end)+1];
                   breaks_MKL{i_break} = c_break;
               end
               pieces_MKL = pieces_MKL + 2;
           end

           coefs_MKL = myppual_mex(MKL_flag_cons,breaks_MKL,Values_MKL,pieces_MKL,order_MKL,dim_MKL,Method,x,extrap_order_MKL,indexDim_MKL,index_x_MKL);
           % the final step
           if strcmp(myppual_flag,'cons')
               v = struct('form','MKLpp','breaks',{breaks_MKL},'Values',Values_MKL,'coefs',coefs_MKL,...
                   'pieces',pieces_MKL,'order',order_MKL,'extrap_order',extrap_order_MKL,'dim',dim_MKL,'Method',Method,'orient','MKLC','thread',-MKL_flag_cons);
           elseif strcmp(myppual_flag,'cons_eval_matrix')
               [v,index_x] = myppual_mex(MKL_flag,breaks_MKL,coefs_MKL,pieces_MKL,order_MKL,dim_MKL,Method,x,left,indexDim_MKL,index_x_MKL);
           end
       end
   elseif strcmp(myppual_flag,'eval_matrix_general')
       dim_MKL = int32(pp.dim); % directly get dim from the structure
       [v,index_x] = myppual_mex(MKL_flag,breaks_MKL,coefs_MKL,pieces_MKL,order_MKL,dim_MKL,Method,x,left,indexDim_MKL,index_x_MKL);
   end    
elseif isobj_griddedInterpolant
    %---------------------------------------------------------
    % Case IV: MATLAB griddedInterpolant spline object with m variables
    %---------------------------------------------------------
    % construct the griddedInterpolant object
    dim = size(pp.Values,1);
    if isstruct(pp)
        if isequal(dim,1) % scalar spline: no need to form a pseudo dimension
            GridVectors = pp.breaks;
        else % vector-valued spline: need to form a pseudo dimension 
            GridVectors = [{(1:dim)},pp.breaks];
        end
        order = pp.order; % the order of the spline
        if isscalar(order),order = repmat(order,[1,length(pp.breaks)]);end
        ExtrapolationOrder = order; % initialization
        if isfield(pp,'ExtrapolationOrder') && (~isempty(pp.ExtrapolationOrder))
            ExtrapolationOrder = pp.ExtrapolationOrder;
            if isscalar(ExtrapolationOrder)
                ExtrapolationOrder = repmat(ExtrapolationOrder,[1,length(pp.breaks)]);
            end
        end
        Method = 'not-a-knot'; % initialization
        if isfield(pp,'Method') && (~isempty(pp.Method))
            Method = pp.Method;
        end        
        if isequal(order(:),4*ones(length(pp.breaks),1)) && strcmp(Method,'not-a-knot')
            GIMethod = 'spline';
        elseif isequal(order(:),4*ones(length(pp.breaks),1)) && strcmp(Method,'pchip')
            % pchip does not support multi-variate problem, so redefine it to be spline
            if isequal(dim,1)
                GIMethod = 'pchip';
            else
                GIMethod = 'spline';
            end
        elseif isequal(order(:),2*ones(length(pp.breaks),1))
            GIMethod = 'linear';
        elseif isequal(order(:),ones(length(pp.breaks),1))
            GIMethod = 'nearest';
        else
            error('The Spline construction method is not supported in griddedInterpolant');
        end
        if isequal(ExtrapolationOrder(:),4*ones(length(pp.breaks),1)) && strcmp(Method,'not-a-knot')
            ExtrapolationMethod = 'spline';
        elseif isequal(ExtrapolationOrder(:),4*ones(length(pp.breaks),1)) && strcmp(Method,'pchip')            
            % pchip does not support multi-variate problem, so redefine it to be spline
            if isequal(dim,1)
                ExtrapolationMethod = 'pchip';
            else
                ExtrapolationMethod = 'spline';
            end
        elseif isequal(ExtrapolationOrder(:),2*ones(length(pp.breaks),1))
            ExtrapolationMethod = 'linear';
        elseif isequal(ExtrapolationOrder(:),ones(length(pp.breaks),1))
            ExtrapolationMethod = 'nearest';
        end                    
        pp_griddedInterpolant = griddedInterpolant(GridVectors,pp.Values,GIMethod,ExtrapolationMethod);
    elseif isa(pp,'griddedInterpolant')
        pp_griddedInterpolant = pp;
    end
    % construction
    if isequal(nargin,1)
        v = pp_griddedInterpolant;
        index_x = [];
        return;
    end
    
    % evaluation    
    if isequal(dim,1) && iscell(x) % scalar-valued spline
        sizev = [1,cellfun('prodofsize',x,'UniformOutput',true)];
        x_griddedInterpolant = x;
    elseif isequal(dim,1) && ismatrix(x) % scalar-valued spline
        sizev = [1,size(x,2)];
        x_griddedInterpolant = x.';
    elseif (dim > 1) && iscell(x) && isempty(indexDim) % use the cell form for evaluation
        sizev = [dim,cellfun('prodofsize',x,'UniformOutput',true)];
        x_griddedInterpolant = [{(1:dim)},x];        
    elseif (dim > 1) && iscell(x) && (~isempty(indexDim)) % need to use the matrix form for evaluation
        sizev = [size(indexDim,1),cellfun('prodofsize',x,'UniformOutput',true)];
        mx = length(x);
        if isequal(mx,1)
            x_mat  = reshape(x{1},1,[]);
        else
            x_ndgrid = cell(1,mx);
            [x_ndgrid{:}] = ndgrid(x{:});
            x_ndgrid = cellfun(@(xx) xx(:),x_ndgrid,'UniformOutput',false);
            x_mat = [x_ndgrid{:}].';
        end
        x_griddedInterpolant = [reshape(indexDim,[],1),...
                                (reshape(repmat(x_mat,[sizev(1) 1]),mx,[])).'];
    elseif (dim > 1) && ismatrix(x) % need to use the matrix form for evaluation
        x_mat = x;
        [mx,nx] = size(x_mat);
        if isempty(indexDim),indexDim = repmat((1:dim).',[1 nx]);end
        % in the matrix form of griddedInterpolant, it takes the transpose of matrix form in Curve Fitting Toolbox
        sizev = [size(indexDim,1),nx];
        x_griddedInterpolant = [reshape(indexDim,[],1),...
                                (reshape(repmat(x_mat,[sizev(1) 1]),mx,[])).'];
    end
    v = pp_griddedInterpolant(x_griddedInterpolant);
    v = reshape(v,sizev);    
else
    error('SPLINES:MYPPUAL:cases','The requested usage case is not currently supported.');
end

end % end function myppual

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Subfunction: myppual1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = myppual1(breaks,coefs,pieces,order,dim,x,left,indexDim,index_x)
%MYPPUAL1 Evaluate univariate function in ppform.

% error check on input dimension
if isempty(x) || (~isvector(x)) || (~isa(x,'double'))
    error('SPLINES:MYPPUAL1:wrongx: X must be a nonempty double vector.');
end

breaks = reshape(breaks,1,numel(breaks)); % reshape breaks to be a row vector 
x = reshape(x,1,numel(x)); % reshape x into a row vector
[mx,nx] = size(x);   %#ok size of x

% check indexDim
if isempty(indexDim)
    mrowIndexDim = dim;
    temp_indexDim = reshape((1:dim),dim,1);
    indexDim = temp_indexDim(:,ones(1,nx));
else    
    mrowIndexDim = size(indexDim,1);    
end

% for each data site, compute its break interval
%{
  Note: If length(breaks)==2 (or equivalently pieces == 1),"breaks(2:end-1)"
  will be empty. As a result, get_index will return index value 1 for any non NaN
  data site.
%}
if isempty(index_x)
    if pieces > 1
        [index_x,NaNx] = get_index(breaks(2:end-1),x,left);
        if ~isempty(NaNx)
            index_x(NaNx) = 1;
        end
    else
        NaNx = [];
        index_x = ones(1,nx);
    end
else
    NaNx = [];
    index_x = reshape(index_x,1,[]);
end

% now go to local coordinates ...
xs = x-breaks(index_x);
%{
  note from JB: The following line is the main difference between "myppual1" and the
  Curve Fitting Toolbox's "ppual1," where we adjust index_x in case PP is vector-valued so that it
  conforms to coefficient index arrangement in pp structure. The arrangement for
  pp structure can be seen from an example with pp.dim = 3 and pp.pieces = 99. There the location of 
  coefficients along each column is arranged according to ndgrid(1:3,1:l00), i.e.,
      for dim index, we have [1;2;3;1;2;3;1;2;3;1;2;3;...;1;2;3],
      for x index, we have [1;1;1;2;2;2;3;3;3;4;4;4;...;99;99;99],
  As a result, 
  dim*(index -1) locates the starting point according to index of x, while indexDim(:) adjusts the
  location further according to which dimension to evaluate
%}
numel_eval = mrowIndexDim*nx; % the total number of evaluation sites
if dim>1 % ... replicate XS and INDEX in case PP is vector-valued ...
    if (mrowIndexDim == 1)
        index = dim*(index_x-1)+indexDim;
    elseif (mrowIndexDim > 1)
        xs = reshape(xs(ones(mrowIndexDim,1),:),1,numel_eval); % replicate
        index = reshape(bsxfun(@plus,dim*(index_x-1),indexDim),1,numel_eval); % use bsxfun to save memory
    end
else
    index = index_x;
end

%{ 
  evaluate piecewise polynominal in an efficient way (See Numerical Recipe, 3rd Edition, Chapter 5.1, for a discussion):
  Note: (a) the local coefficients of ppform in MATLAB is organized from higher order to lower order, e.g., the order of cubic spline is  
            coefs(:,1) = cubic coef; coefs(:,2) = quadratic coef; coefs(:,3) = linear coef; coefs(:,4) = intercept
        (b) we single out cubic, linear, quadratic (in that order based on the common usage pattern) to avoid a direct loop.
%}
if isequal(order,4) % the cubic spline case
    coefs_all = coefs(index,:);
    v = xs.*(xs.*(xs.*(coefs_all(:,1).') + coefs_all(:,2).') + coefs_all(:,3).') + coefs_all(:,4).';    
%     v = xs.*(xs.*(xs.*(coefs(index,1).') + coefs(index,2).') + coefs(index,3).') + coefs(index,4).';
elseif isequal(order,2) % the linear spline case
    v = xs.*(coefs(index,1).') + coefs(index,2).';
elseif isequal(order,3) % the quadratic spline case
    v = xs.*(xs.*(coefs(index,1).') + coefs(index,2).') + coefs(index,3).';
elseif isequal(order,1) % the constant spline case
    v = coefs(index,1).'; % the constant term
    %{
        Deal with NaN: the nested multiplication takes care of NaN automatically if order > 1
        (any multiplication with NaN leads to NaN ). But we need an additional step to recover 
        NaN if order ==1, otherwise the value will be coefs(1,1) incorrectly.
    %}
    if ~isempty(NaNx), v(NaNx) = NaN; end    
else % higher-order (order >=5) spline: use nested multiplication
    v = coefs(index,1).';
    for i=2:order
        v = xs.*v + coefs(index,i).';
    end
end
% reshape the result from a row vector to the specified dimension
v = reshape(v,mrowIndexDim,nx);

end % end function myppual1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Subfunction: myppual_matrix_general
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = myppual_matrix_general(breaks,coefs,pieces,order,dim,x,left,indexDim,index_x)
%myppual_matrix_general Evaluate multivariate pp function in matrix form.

[mx,nx] = size(x);
% check indexDim
if isempty(indexDim)
    mrowIndexDim = dim;
    temp_indexDim = reshape((1:dim),dim,1);
    indexDim = temp_indexDim(:,ones(1,nx));
else    
    mrowIndexDim = size(indexDim,1);    
end

% locate the scattered data in the break sequences:
if isempty(index_x)
    index_x = zeros(mx,nx);  % index for table lookup on each dimension
    for i=1:mx
        index_x(i,:) = get_index(breaks{i}(2:end-1),x(i,:),left{i});
    end
end

% ... and now set up lockstep polynomial evaluation
% %%First,  select the relevant portion of the coefficients array.
% This has the additional pain that now there are order(i) coefficients
% for the i-th univariate interval.
% The coefficients sit in the (m+1)-dimensional array COEFS, with
% the (i+1)st dimension containing the coefficients in the i-th
% dimension, and organized to have first the highest coefficients
% for each interval, then the next-highest, etc (i.e., as if coming
% from an array of size [pieces(i),order(i)]).
% index_x(:,j) is the index vector for the lower corner of j-th point
% The goal is to extract, for the j-th point, the requisite coefficients
% from the equivalent one-dimensional array for COEFS, computing a
% base index from index_x(:,j), and adding to this the same set of offsets
% computed from the pieces(i) and order(i).
sizec = [dim,pieces.*order]; % size(coefs)
temp = pieces(1)*(0:order(1)-1)'; % starting value of each coefficent - 1
for i=2:mx
    lt = length(temp(:,1));
    temp = [repmat(temp,order(i),1), ...
        reshape(repmat(pieces(i)*(0:order(i)-1),lt,1),order(i)*lt,1)];
end
% calculate "offset" as the linear index of 1st piece for each coefs
lt = size(temp,1);

% the following vectorized code changes subscript to linear indexing, and is equivalent to the following:
% temp = num2cell([ones(lt,1),1+temp],1);
% temp1 = sub2ind(sizec,temp{:});
% admittedly, my code appears to be very cryptic (even for my future self!). Why not do simple coding, then? It turns out that the overhead 
% for calling num2cell and sub2ind is quite large, so I wrote the vectorized code purely for efficiency reasons.
temp = [ones(lt,1),1+temp];
temp_stride = [1 cumprod(sizec(1:end-1))]; % the stride for each dimension
temp1 = reshape(1 + sum(bsxfun(@times,temp - 1,temp_stride),2),1,lt); % get linear index
offset = reshape(temp1(ones(mrowIndexDim,1),:),[],1);

% calculate "base" as the linear index of different pieces. See above for the explanation to the apparently cryptic code.
temp = [ones(nx,1) index_x.'];
temp1 = reshape(1 + sum(bsxfun(@times,temp - 1,temp_stride),2),1,[]); % get linear index
% notice that we substract 2 instead of 1: one for dims index, one for subscript
base = repmat(bsxfun(@plus,temp1,indexDim) - 2,prod(order),1);

v = reshape(coefs(bsxfun(@plus,base,offset)),[mrowIndexDim,order,nx]);

% ... then do a version of local polynomial evaluation
for i=mx:-1:1
    xs = x(i,:) - breaks{i}(index_x(i,:)); % local value
    otherk = mrowIndexDim*prod(order(1:i-1));
    v = reshape(v,[otherk,order(i),nx]);
    for j=2:order(i)
        v(:,1,:) = reshape(bsxfun(@times,reshape(v(:,1,:),otherk,nx),xs),[otherk,1,nx])...
            + v(:,j,:);
    end
    v(:,2:order(i),:) = [];
end
v = reshape(v,mrowIndexDim,nx);

end % end myppual_matrix_general

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Subfunction: myppual_cell_general
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = myppual_cell_general(breaks,coefs,pieces,order,dim,x,left,indexDim,index_x)
%v = myppual_cell_general Evaluate multivariate pp function in cell form.

mx = length(x);
% check indexDim
if isempty(indexDim)
    mrowIndexDim = dim;
else    
    mrowIndexDim = size(indexDim,1);    
end

v = coefs; sizev = [dim,pieces.*order];
nsizev = zeros(1,mx);
for i=mx:-1:1
    nsizev(i) = length(x{i}(:)); dd = prod(sizev(1:mx));
    if (i > 1) || (isempty(indexDim))  % for i > 1 or empty indexDim, evaluate for all vector dimensions
        v = reshape(myppual1(breaks{i},reshape(v,dd*pieces(i),order(i)),pieces(i),order(i),dd,...
            x{i}, left{i},[],index_x{i}),...
            [sizev(1:mx),nsizev(i)]);
    elseif isequal(i,1) && (~isempty(indexDim))    % use myppual1 for the first dimension
        % convert subindex into linear index if other dimensions are not singleton
        sizeIndex2 = prod(sizev(2:mx));
        if (sizeIndex2 > 1)
            indexDim = sub2ind([sizev(i) sizeIndex2],...
                reshape(indexDim,[mrowIndexDim,nsizev(i),sizeIndex2]),...
                reshape(repmat(1:sizeIndex2,[mrowIndexDim*nsizev(i),1]),[mrowIndexDim,nsizev(i),sizeIndex2]));
            indexDim = reshape(permute(indexDim,[1 3 2]),[mrowIndexDim*sizeIndex2,nsizev(i)]);
        end
        v = reshape(myppual1(breaks{i},reshape(v,dd*pieces(i),order(i)),pieces(i),order(i),dd,...
            x{i}, left{i},indexDim,index_x{i}),...
            [mrowIndexDim,sizev(2:mx),nsizev(i)]);
    end
    sizev(mx+1) = nsizev(i);
    if mx>1
        v = permute(v,[1,mx+1,2:mx]); sizev(2:mx+1) = sizev([mx+1,2:mx]);
    end
end

end % end subfunction myppual_cell_general

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Subfunction: myppual_cell_cons_eval
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = myppual_cell_cons_eval(breaks,Values,pieces,order,dim,x,left,indexDim,index_x,Method,ExtrapolationOrder)
%v = myppual_cell_cons_eval Combined m-variate Spline Construction and Evaluation of pp splines in cell form.

mx = length(x);
% check indexDim
if isempty(indexDim)
    mrowIndexDim = dim;
else    
    mrowIndexDim = size(indexDim,1);    
end

v = Values; sizev = size(Values);
nsizev = zeros(1,mx); % the dimension after evaluation
for i=mx:-1:1
    nsizev(i) = length(x{i}(:)); dd = prod(sizev(1:mx));
    % form the spline: currently only support csape. Will be updated in the future subject to needs
    if isequal(order(i),4)
        ppV = csape(breaks{i},reshape(v,[dd,sizev(1+mx)]),Method{i});
    else
        error('Currently we only support spline construction in cubic form');
    end
    % deal with extrapolation
    if (ExtrapolationOrder(i) >= 0) && (ExtrapolationOrder(i) < order(i))
        ppV = fnxtr(ppV,ExtrapolationOrder(i));
        % now the following parts may change after fnxtr, so we redefine them.
        breaks{i} = ppV.breaks; pieces(i) = ppV.pieces;
    end
    % evaluate the spline
    if i > 1 || (isempty(indexDim))  % for i > 1, evaluate for all vector dimensions
        v = reshape(myppual1(breaks{i},reshape(ppV.coefs,dd*pieces(i),order(i)),pieces(i),order(i),dd,...
            x{i}, left{i},[],index_x{i}),...
            [sizev(1:mx),nsizev(i)]);
    elseif isequal(i,1)  && (~isempty(indexDim))  % use myppual1 for the first dimension
        % convert subindex into linear index if other dimensions are not singleton
        sizeIndex2 = prod(sizev(2:mx));
        if (sizeIndex2 > 1)
            indexDim = sub2ind([sizev(i) sizeIndex2],...
                reshape(indexDim,[mrowIndexDim,nsizev(i),sizeIndex2]),...
                reshape(repmat(1:sizeIndex2,[mrowIndexDim*nsizev(i),1]),[mrowIndexDim,nsizev(i),sizeIndex2]));
            indexDim = reshape(permute(indexDim,[1 3 2]),[mrowIndexDim*sizeIndex2,nsizev(i)]);
        end
        v = reshape(myppual1(breaks{i},reshape(ppV.coefs,dd*pieces(i),order(i)),pieces(i),order(i),dd,...
            x{i}, left{i},indexDim,index_x{i}),...
            [mrowIndexDim,sizev(2:mx),nsizev(i)]);
    end
    sizev(mx+1) = nsizev(i);
    if mx>1
        v = permute(v,[1,mx+1,2:mx]);
        sizev(2:mx+1) = sizev([mx+1,2:mx]);
    end
end

end % end subfunction myppual_cell_cons_eval

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Subfunction: pplin1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pp = pplin1(breaks,Values)
%pplin1 constructs a vector-valued 1-D linear spline of ppform

d = size(Values,1);
breaks = reshape(breaks,1,numel(breaks));
coefs2 = diff(Values,1,2)./repmat(diff(breaks,1,2),[d,1]); % slope coefficient
coefs1 = Values(:,1:(end-1));  % intercept coefficient
coefs = [coefs2(:),coefs1(:)]; % construct coefficient matrix
% construct a pp form spline
pp = struct('form','pp','breaks',breaks,'coefs',coefs,'pieces',length(breaks) - 1,'order',2,'dim',d);
end % end subfunction pplin1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Subfunction: rs3dgHermite
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pp = rs3dgHermite(breaks,Values,Der1)
%rs3DGHermite constructs a vector-valued ppform 1-D rational Hermite spline of order 3 (i.e., quadratic spline) by Delbourgo and Gregory 

d = size(Values,1); % the vector-valued dimension
breaks = reshape(breaks,1,numel(breaks));
pieces = length(breaks) - 1; % the number of pieces of the splines

% calculate the coefficients as in Cai and Judd (2012, pp162)
c1 = Values(:,1:(end-1));  % the value coefficient
diff_breaks = repmat(diff(breaks,1,2),[d,1]);
c2 = diff(Values,1,2)./diff_breaks; % slope for linear interpolation
c3 = Der1(:,1:(end-1)) - c2; % difference between left-end derivative and linear slope
c4 = Der1(:,2:end) - c2; % difference between right-end derivative and linear slope

% calculate the coefficient of the piecewise polynomial part
coefs3_pp = zeros(d,pieces); % coefficient for quadratic term
coefs2_pp = c2; % slope coefficient
coefs1_pp = c1;  % intercept coefficient

% calculate the coefficient of the numerator
coefs3_numerator = c3.*c4; % coefficient for quadratic term
coefs2_numerator = coefs3_numerator.*(-diff_breaks); % linear coefficient
coefs1_numerator = zeros(d,pieces);  % intercept coefficient

% calculate the coefficient of the denominator
coefs3_denominator = zeros(d,pieces); % coefficient for quadratic term
coefs2_denominator = c3 + c4; % linear coefficient
coefs1_denominator = c4.*(-diff_breaks);  % intercept coefficient

% construct coefficient matrix by stacking all three components together
coefs = [coefs3_pp,coefs2_pp,coefs1_pp;...
         coefs3_numerator,coefs2_numerator,coefs1_numerator;...
         coefs3_denominator,coefs2_denominator,coefs1_denominator];
coefs = reshape(coefs,[3*d*pieces,3]);     
% construct a pp form spline: notice that the vector dimension is multiplied by 3
pp = struct('form','prp','breaks',breaks,'coefs',coefs,'pieces',pieces,'order',3,'dim',d,'Method','dgHermite');
end % end subfunction rs3DGHermite

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Subfunction: pptransformation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ppNew = pptransformation(pp,conversion_flag)
%pp2MKLpp converts MATLAB pp form to MKLpp form

ppNew = pp;
breaks = pp.breaks;
pieces = pp.pieces;
order = pp.order;
dim = prod(pp.dim);
% conversion
mx = length(pieces);
idx_cell = cell(1,1+mx);
if strcmp(conversion_flag,'pp2MKLpp')
    % calculate the conversion index for each dimension
    idx_cell{1} = (1:dim); % original order in the vector-value dimension
    for ii = 1:mx % each variable dimension
        % use sub2ind to convert the different order arrangement between pp and MKL        
        idx_row = reshape(repmat((1:pieces(ii)),[order(ii) 1]),[order(ii)*pieces(ii),1]); % the piece dimension
        idx_column = repmat((order(ii):-1:1).',[pieces(ii),1]); % order dimension
        idx = sub2ind([pieces(ii),order(ii)],idx_row,idx_column); % linear index
        % record the new order in linear index
        idx_cell{ii+1} = idx;
    end
    % reshape coefs to the dimension of [dim,pieces.*order]
    coefs = reshape(pp.coefs,[dim,pieces.*order]);
    % reorder coefs array according to polynomial order of MKL and permute it from column-major to row-major storage format
    coefs_MKL = permute(coefs(idx_cell{:}),((mx+1):-1:1));
    % change fields
    ppNew.form = 'MKLpp';
    ppNew.coefs = coefs_MKL;
    ppNew.dim = int32(dim);
    ppNew.pieces = int32(pieces);
    ppNew.order = int32(order);
    ppNew.orient = 'MKLC';    
    % permute Values if present
    if isfield(pp,'Values') && (~isempty(pp.Values))
        sizValues = size(pp.Values);
        Values = reshape(pp.Values,[dim,sizValues(2:end)]);
        ppNew.Values = permute(Values,((mx+1):-1:1));
    end
    % change univariate spline to cell breaks
    if ~iscell(breaks),ppNew.breaks = {breaks};end
elseif strcmp(conversion_flag,'MKLpp2pp') % from MKLpp to pp: placeholder as of now
    error('The conversion from MKLpp to pp will be supported at a future version.');
end
end % end function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Subfunction from original "ppual": get_index
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [index,NaNx] = get_index(mesh,sites,left)
%GET_INDEX appropriate mesh intervals for given ordered data sites

if isempty(left)||left(1)~='l'
   [~,index] = histc(sites,[-inf,mesh,inf]);
   NaNx = find(index==0); 
   index = min(index,numel(mesh)+1);
else
   [~,index] = histc(-sites,[-inf,-fliplr(mesh),inf]);
   NaNx = find(index==0); 
   index = max(numel(mesh)+2-index,1);
end

end % end function get_index

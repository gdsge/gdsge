classdef asg <  handle
    % MATLAB interface to C++ AdapitveSparseGrid class
    % @author: Wenlan Luo, luowenlan@gmail.com
    properties (Access = public, Hidden = false)
        objectHandle;
        numDim;
        numVec;
        numArray;
        inputGrids;
        stateMin;
        stateMax;
        stateRange;
    end
    
    methods
        %% Constructor - Create the C++ class
        function this = asg(inputGrids, numVec, numArray)
            numDim = length(inputGrids);
            this.numDim = numDim;
            this.numVec = numVec;
            this.numArray = numArray;
            this.inputGrids = inputGrids;
            this.stateMin = zeros(numDim,1);
            this.stateMax = zeros(numDim,1);
            this.stateRange = zeros(numDim,1);
            for i=1:length(inputGrids)
                this.stateMin(i) = min(inputGrids{i});
                this.stateMax(i) = max(inputGrids{i});
                this.stateRange(i) = max(inputGrids{i}) - min(inputGrids{i});
            end
            
            this.objectHandle = asg_mex('new', numDim, numVec, numArray, this.stateMin, this.stateRange);
        end
        
        %% Destructor
        function delete(this)
            if ~isempty(this.objectHandle)
                asg_mex('delete', this.objectHandle);
            end
            this.objectHandle=[];
        end
        
        %% Assert not empty
        function assert_not_empty(this)
            if isempty(this.objectHandle)
                error('Object already deleted.');
            end
        end
        
        function [evalArrayIdx, evalGrids, evalGridsLength, evalGridsUnscaled] = get_eval_grids(this, threshhold)
            % [evalArrayIdx, evalGrids, evalGridsLength] = GET_EVAL_GRIDS(threshhold)
            % Get grids proposed to be evaluated
            assert_not_empty(this);
            [evalGrids, evalGridsLength, evalGridsUnscaled] = asg_mex('get_eval_grids',this.objectHandle,threshhold);
            % Set the vector index
            evalGridsTotalLength = sum(evalGridsLength);
            evalArrayIdx = zeros(1,evalGridsTotalLength);
            ptr = 0;
            for i=1:this.numArray
                currentGridsLength = evalGridsLength(i);
                evalArrayIdx(ptr+1:ptr+currentGridsLength) = i;
                ptr = ptr+currentGridsLength;
            end
        end
        
        function newAsg = copy(this)
            newAsg = asg.convert_from_struct(this.convert_to_struct);
        end
        
        function currentLevel = get_current_level(this)
            % currentLevel = GET_CURRENT_LEVEL()
            % Return the current level of the interpolation
            assert_not_empty(this);
            currentLevel = asg_mex('get_current_level',this.objectHandle);
        end
        
        function push_eval_results(this, value)
            % PUSH_EVAL_RESULTS(value)
            % Push evaluation results corresponding to grids returned from
            % last get_eval_grids method
            assert_not_empty(this);
            asg_mex('push_eval_results', this.objectHandle, value);
        end
        
        function push_eval_results_at_grids(this, evalArrayIdx, evalGridsUnscaled, value, currentLevel)
            % PUSH_EVAL_RESULTS_AT_GRIDS(evalArrayIdx, evalGrids, value, currentLevel)
            % Push and overwrites all evalGrids and corresponding value, setting the
            % maxLevel to currentLevel
            assert_not_empty(this);
            asg_mex('push_eval_results_at_grids', this.objectHandle, evalArrayIdx, evalGridsUnscaled, value, currentLevel);
        end
        

        
        function [grids, surplus, levels, unscaledGrids] = get_grids_info(this)
            % [grids, surplus, order, level, unscaledGrids] =
            % GET_GRIDS_INFO()
            % Get grids, surplus, order when added, level information.
            % The grids are scaled using grid range.
            assert_not_empty(this);
            [unscaledGrids, surplus, levels] = asg_mex('get_grids_info', this.objectHandle);
            grids = cell(1,this.numArray);
            for i=1:this.numArray
                if ~isempty(unscaledGrids{i})
                    grids{i} = unscaledGrids{i}.*this.stateRange(:) + this.stateMin(:);
                end
            end
        end
        
        function [grids, surplus, levels, unscaledGrids] = get_grids_info_at_level(this, level)
            % [grids, surplus, order, level, unscaledGrids] = GET_GRIDS_INFO_AT_LEVEL(level)
            % Get grids info at current level
            assert_not_empty(this);
            [grids, surplus, levels, unscaledGrids] = get_grids_info(this);
            for i=1:this.numArray
                indicatorAtLevel = sum(levels{i},1)==level;
                
                grids{i}(:,~indicatorAtLevel) = [];
                surplus{i}(:,~indicatorAtLevel) = [];
                levels{i}(:,~indicatorAtLevel) = [];
                unscaledGrids{i}(:,~indicatorAtLevel) = [];
            end
        end
        
        function evalResults = eval(this, arrayIdx, vecIdx, sites)
            % evalResults = EVAL(arrayIdx, vecIdx, sites)
            % Evaluate interpolation at sites, at arrayIdx and vecIdx
            % sites is of size [numDim,nSites]
            assert_not_empty(this);
            evalResults = asg_mex('eval', this.objectHandle, arrayIdx, vecIdx, sites);
        end
        
        function evalResults = eval_vec(this, arrayIdx, sites)
            % evalResults = EVAL(arrayIdx, sites)
            % Evaluate interpolations at all vecIdx at sites, at arrayIdx
            % sites is of size [numDim,nSites]
            assert_not_empty(this);
            evalResults = asg_mex('eval_vec', this.objectHandle, arrayIdx, sites);
        end
        
        function lowerBoundIdx = find_lower_bound_index(this, arrayIdx, sites)
            % lowerBoundIdx = FIND_LOWER_BOUND_INDEX(arrayIdx, sites)
            % Find the nearest point of sites
            % sites is of size [numDim,nSites]
            assert_not_empty(this);
            lowerBoundIdx = asg_mex('find_lower_bound_index', this.objectHandle, arrayIdx, sites);
        end
        
        function [gridsCurrentLevel, gridsNextLevel] = get_grids_current_and_next(this)
            % [gridsCurrentLevel, gridsNextLevel] =
            % GET_GRIDS_CURRENT_AND_NEXT
            % Find the internal grids at acurrent and next level
            [gridsCurrentLevel, gridsNextLevel] = asg_mex('get_grids_current_and_next', this.objectHandle);
        end
        
        function asgStruct = convert_to_struct(this)
            % asgStruct = CONVERT_TO_STRUCT()
            % Return all information in the asg object as a struct
            assert_not_empty(this);
            
            [gridsCurrentLevel, gridsNextLevel] = this.get_grids_current_and_next();
            [grids, surplus, levels, unscaledGrids] = this.get_grids_info();
            currentLevel = this.get_current_level();
            
            asgStruct.numDim = this.numDim;
            asgStruct.numVec = this.numVec;
            asgStruct.numArray = this.numArray;
            asgStruct.inputGrids = this.inputGrids;
            asgStruct.stateMin = this.stateMin;
            asgStruct.stateMax = this.stateMax;
            asgStruct.stateRange = this.stateRange;
            
            asgStruct.currentLevel = currentLevel;
            asgStruct.gridsCurrentLevel = gridsCurrentLevel;
            asgStruct.gridsNextLevel = gridsNextLevel;
            
            asgStruct.unscaledGrids = unscaledGrids;
            asgStruct.surplus = surplus;
            asgStruct.levels = levels;
        end
        
        function reconstruct_internal(this, currentLevel, unscaledGrids, surplus, levels, gridsCurrentLevel, gridsNextLevel)
            asg_mex('reconstruct_internal', this.objectHandle, ...
                currentLevel, unscaledGrids, surplus, levels, gridsCurrentLevel, gridsNextLevel);
        end
    end
    %% Static method
    methods(Static)
        function interp = construct_from_struct(asgStruct)
            % interp = CONSTRUCT_FROM_STRUCT(asgStruct)
            % Construct the asg object and return an asg wrapper
            
            % Call the constructor
            interp = asg(asgStruct.inputGrids, asgStruct.numVec, asgStruct.numArray);
            
            % Manually assign the internal objects
            interp.reconstruct_internal(asgStruct.currentLevel, asgStruct.unscaledGrids, ...
                asgStruct.surplus, asgStruct.levels, asgStruct.gridsCurrentLevel, ...
                asgStruct.gridsNextLevel);
        end
        
        function [MAX_DIM,MAX_NVEC,MAX_LEVEL] = get_mex_constants()
            [MAX_DIM,MAX_NVEC,MAX_LEVEL] = asg_mex('get_mex_constants');
        end
        
        function [maxMetric,maxMetricVec] = compute_inf_metric(interp1, interp2)
            % metric = COMPUTE_INF_METRIC(interp1, interp2)
            % eval interp1 and interp2 at interp1's grids, and
            % calculate inf norm
            if interp1.numDim ~= interp2.numDim
                error('numDim does not agree');
            end
            %{
            if interp1.numVec ~= interp2.numVec
                error('numVec does not agree');
            end
            %}
            if interp1.numArray ~= interp2.numArray
                error('numArray does not agree');
            end
            
            %{
            if any(interp1.stateMin~=interp2.stateMin)
                error('stateMin does not agree');
            end
            
            if any(interp1.stateMax~=interp2.stateMax)
                error('stateMax does not agree');
            end
            %}
            
            % get grids
            % 
            grids1 = interp1.get_grids_info;
            grids2 = interp2.get_grids_info;
            
            maxMetric = 0;
            maxMetricVec = zeros(min(interp1.numVec,interp2.numVec),1);
            for i=1:interp1.numArray
                grids = [grids1{i},grids2{i}];
                
                gridsSize = size(grids,2);
                
                fval1 = interp1.eval_vec(i*ones(1,gridsSize), grids1{i});
                fval2 = interp2.eval_vec(i*ones(1,gridsSize), grids1{i});
                % Get the minimum numVec
                minNumVec = min(size(fval2,1),size(fval1,1));
                % Compute metric
                fMetric = max(abs(fval1(1:minNumVec,:) - fval2(1:minNumVec,:)), [], 2);
                maxMetricVec = max(maxMetricVec,fMetric);
                
                maxMetric = max(max(fMetric(:)), maxMetric);
            end
        end
        
        function differenceInterp = calculate_difference(interp1, interp2, maxLevel, threshhold)
            % difference = CALCULATE_DIFFERENCE(interp1, interp2)
            % Calcualte the difference function between two asg
            % interpolations
            
            % Check valid input
            if interp1.numDim ~= interp2.numDim
                error('numDim does not agree');
            end
            %{
            if interp1.numVec ~= interp2.numVec
                error('numVec does not agree');
            end
            %}
            if interp1.numArray ~= interp2.numArray
                error('numArray does not agree');
            end
            
            minNumVec = min(interp1.numVec, interp2.numVec);
            
            differenceInterp = asg(interp1.inputGrids, minNumVec, interp1.numArray);
            
            while differenceInterp.get_current_level<maxLevel
                % Get cosntruction grids
                [evalArrayIdx,evalGrids] = differenceInterp.get_eval_grids(threshhold);
                % Evaluations
                fval1 = interp1.eval_vec(evalArrayIdx,evalGrids);
                fval2 = interp2.eval_vec(evalArrayIdx,evalGrids);
                
                difference = fval1(1:minNumVec,:)-fval2(1:minNumVec,:);
                
                differenceInterp.push_eval_results(difference);
            end
        end
        
        function [maxMetric,maxMetricVec] = compute_inf_metric_accurate(interp1, interp2, maxLevel, threshhold)
            differenceInterp = asg.calculate_difference(interp1, interp2, maxLevel, threshhold);
            
            % Evaluate at difference
            grids = differenceInterp.get_grids_info;
            maxMetricVec = zeros(min(interp1.numVec,interp2.numVec),1);
            
            for i=1:interp1.numArray
                grid = grids{i};
                difference = differenceInterp.eval_vec(i*ones(1,size(grid,2)), grid);
                maxAbsDifference = max(abs(difference),[],2);
                
                maxMetricVec = max(maxAbsDifference, maxMetricVec);
            end
            
            maxMetric = max(maxMetricVec(:));
        end
        
        function expectedInterp = integrate_expectation(interp, trans, maxLevel, threshhold)
            % expected_interp = INTEGRATE_EXPECTATION(interp, trans)
            % Calculate integration w.r.t. the array dimension using
            % transition matrix trans
            
            % Check valid input
            if size(trans,1)~=interp.numArray || size(trans,2)~=interp.numArray
                error('trans size does not agree with numArray');
            end
            
            expectedInterp = asg(interp.inputGrids, interp.numVec, interp.numArray);
            
            while expectedInterp.get_current_level<maxLevel
                % Get construction grids
                [evalArrayIdx,evalGrids] = expectedInterp.get_eval_grids(threshhold);
                expectedFval = zeros(interp.numVec,size(evalGrids,2));
                for i_current=1:interp.numArray
                    % Get the index that is within the array
                    % evalArrayIdxAtCurrentArray = evalArrayIdx(evalArrayIdx==i_current);
                    evalGridsAtCurrentArray = evalGrids(:,evalArrayIdx==i_current);
                    lenEvalGrids = size(evalGridsAtCurrentArray,2);
                    
                    % Integrate over future
                    % Replicate current grids to all future array index
                    evalGridsAtPrime = repmat(evalGridsAtCurrentArray,[1,interp.numArray]);
                    arrayIndexAtPrime = repmat([1:interp.numArray],[lenEvalGrids,1]);
                    % Reshape array index to a row vector
                    arrayIndexAtPrime = reshape(arrayIndexAtPrime,1,[]);
                    
                    % Do a vector evaluation
                    fvalAtPrime = interp.eval_vec(arrayIndexAtPrime, evalGridsAtPrime);
                    
                    % Now fvalAtFuture has the order (vec, sites, array),
                    % applying the transition matrix carefully
                    expectedFvalAtCurrentArray = reshape(fvalAtPrime,[],interp.numArray) * trans(i_current,:)';
                    
                    % Now expectedFval has shape (vec, sites, 1), Shape
                    % expectedFval to the right shape and push to
                    % construction
                    expectedFval(:,evalArrayIdx==i_current) = reshape(expectedFvalAtCurrentArray, [interp.numVec, lenEvalGrids]);
                end
                expectedInterp.push_eval_results(expectedFval);
            end
        end
    end
end

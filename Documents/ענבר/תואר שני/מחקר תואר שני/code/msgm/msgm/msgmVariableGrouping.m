function [vg, mapFineToCoarse] = msgmVariableGrouping(G, param, bInitialized)
% msgmVariableGrouping(G, bInitialized) select a variable grouping for the coarsening stage
%
% input:
%
%       - G : graphical model
%       - bInitialized : boolean flag for whether labels are initialized
%
% output:
%       - vg, struct array
%           vg(i).seed  : the i-th seed
%           vg(i).vars  : variables in the i-th group (including seed)
%           vg(i).edges : (some of the...) edges in the i-th group
%           vg(i).brev  : boolean flag specifying if the edge is 'reversed'
%           vg(i).map   : interpolation rule, set later
%       
%       - mapFineToCoarse : mapping of variables from a fine scale to a
%                           coarse scale, i.e. mapFineToCoarse(v1) = vc
%                           means that fine-variable indexed by v1 is
%                           mapped to a coarse-variable indexed by vc

    % output data structures
    vg = cell(size(G.u, 1), 1);
    mapFineToCoarse = zeros(size(G.u,1),1);

    % VARS and SEEDS
    vbVars = true(size(G.u,1),1);
    vbSeeds = false(size(G.u,1),1);
    
    % local-conditional-entropy scores
    % vbReverseEdge defines whether the pairwise is stored
    % in G.p as (v1,v2) or as (v2,v1)
    [orderedEdgeList, vbReverseEdge] = msgmScoreEdges(G, param, bInitialized);
    
    % assign SEED variables and their respective group
    iEdge = 0;
    iGroup = 0;
    while (any(vbVars) &&  (iEdge < numel(orderedEdgeList)))
       
        % get the next edge
        iEdge = iEdge + 1;
        v1 = G.adj(orderedEdgeList(iEdge), 1);
        v2 = G.adj(orderedEdgeList(iEdge), 2);
        if (vbReverseEdge(iEdge))
            % the relevant direction is (v2,v1)
           
            v_ = v1;
            v1 = v2;
            v2 = v_;
        end
        v_ = v1;
        v1 = v2;
        v2 = v_;
        
        % verify that v2 has not been assigned
        % ..and that v1 is either seed or not assigned
        if (vbVars(v2) && (vbVars(v1) || vbSeeds(v1)))
            
            % set v1 to be v2's seed
            if (vbSeeds(v1))
                % v1 is already a seed variable
                
                v1c = mapFineToCoarse(v1);
                
                % update v1c's group
                vg{v1c}.vars = cat(1, vg{v1c}.vars, v2);
                vg{v1c}.edges = cat(1, vg{v1c}.edges, orderedEdgeList(iEdge));
                vg{v1c}.brev = cat(1, vg{v1c}.brev, vbReverseEdge(iEdge));
                
                % update v2's coarse representative
                mapFineToCoarse(v2) = v1c;
                
            else
                % v1 is a new seed variable
                
                iGroup = iGroup + 1;
                
                % construct v1c's group              
                v1group.seed = v1;
                v1group.vars = [v1; v2];
                v1group.edges = orderedEdgeList(iEdge);
                v1group.brev = vbReverseEdge(iEdge);
                vg{iGroup} = v1group;
                
                % update v1,v2 coarse representative
                mapFineToCoarse(v1) = iGroup;
                mapFineToCoarse(v2) = iGroup;

            end
            
            % update seeds and vars
            vbVars([v1, v2]) = false;
            vbSeeds(v1) = true;
        end
    end
    
    % collect "leftover" variables
    leftoverVars = find(vbVars);
    for i = 1 : numel(leftoverVars)
        
        iGroup = iGroup + 1;
       
        % the group is singleton
        v1group.seed = leftoverVars(i);
        v1group.vars = leftoverVars(i);
        v1group.edges = [];
        v1group.brev = [];
        vg{iGroup} = v1group;
        mapFineToCoarse(leftoverVars(i)) = iGroup;
    end
    
    % clear pre-allocated space, make vg a struct array
    vg(iGroup + 1 : end) = [];
    vg = cat(1, vg{:});
end


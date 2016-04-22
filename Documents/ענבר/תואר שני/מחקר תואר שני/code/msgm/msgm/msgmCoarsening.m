function [Gc, xc, vg] = msgmCoarsening(G, param, x)
% msgmCoarsening(G, x) apply coarsening by variable-grouping
%

    % select a variable-grouping
    [vg, mapFineToCoarse] = msgmVariableGrouping(G, param, any(x));

    % set an interpolation rule
    [vg, mapInterpolation] = msgmSetInterpolationRule(G, x, vg);
    
    % set the coarse potentials
    Gc = msgmSetCoarsePotentials(G, vg, mapFineToCoarse, mapInterpolation);
    
    % intialize coarse scale's labeling
    xc = msgmInherit(x, [vg.seed]);
    
    msgmEnergyAssert(G, x, Gc, xc);
end


%% Coarsening subroutines

function [vg, mapInterpolation] = msgmSetInterpolationRule(G, x, vg)
% msgmSetInterpolationRule set an interpolation-rule, from labeling of a
% coarse scale to labeling of a fine scale
%
% output:
%
%   - vg : appends a 'map' field for the interpolation rule
%
%   - mapInterpolation : (optinal, for internal use) interpolation table,
%                        used for slightly improving the running time
%

    % intialize the interpolation-rule table
    mapInterpolation = zeros(size(G.u,1), G.numLabels);    

    % itearte seed variables
    for i = 1 : numel(vg)
               
        % allocate space for the interpolation rule of the i-th group
        map = zeros(numel(vg(i).vars), G.numLabels);
        
        % iterate variables in the seed's group
        for j = 1 : numel(vg(i).vars)
           
            if (vg(i).vars(j) == vg(i).seed)
                % seed variable gets the label
                % of the coarse representative (Eq. (2))
                
                map(j,:) = 1 : G.numLabels;
                mapInterpolation(vg(i).vars(j),:) = 1 : G.numLabels;
                
            else
                % find the minimizer (Eq. (3))

                pairwise = squeeze(G.p(:,:,vg(i).edges(j-1)));
                if (~vg(i).brev(j-1))
                    % transpose the pairwise s.t. seed is on 2nd dim

                    pairwise = pairwise';
                end
                pairwise = bsxfun(@plus, pairwise, G.u(vg(i).vars(j),:)');
                [~, map_] = min(pairwise,[],1);            

                if any(x)
                    % labels are intiailized,
                    % reset the interpolatin rule (Eq. (6))

                    map_(x(vg(i).seed)) = x(vg(i).vars(j));
                end
                
                map(j,:) = map_;
                mapInterpolation(vg(i).vars(j),:) = map_;
            end
        end
        
        vg(i).map = map;
    end
end

function Gc = msgmSetCoarsePotentials(G, vg, mapFineToCoarse, mapInterpolation)

    % keep track which edges have been accounted for
    vbTouch = false(size(G.p,3), 1);

    % set the coarse unary terms (Eq. (4))
    uc = zeros(numel(vg), G.numLabels);
    for i = 1 : numel(vg)
        
        uGroup = zeros(1, G.numLabels);
        
        % sum j's unary term, considering the interpolation
        % this is the first term in Eq. (4)        
        for j = 1 : numel(vg(i).vars)
            
            vj = vg(i).vars(j);
            uGroup = uGroup + G.u(vj, vg(i).map(j,:));
        end
        
        % sum the pairwise energy, considering the interpolation
        % this is (part of) the second term in Eq. (4)
        for j = 1 : numel(vg(i).edges)
            
            vbTouch(vg(i).edges(j)) = true;
            pairwise = G.p(:,:,vg(i).edges(j));
            
            % interpolation rule applied to v1, v2
            v1 = G.adj(vg(i).edges(j), 1);
            v2 = G.adj(vg(i).edges(j), 2);
            map1 = mapInterpolation(v1,:);
            map2 = mapInterpolation(v2,:);
%             map1 = vg(i).map(vg(i).vars == v1,:);
%             map2 = vg(i).map(vg(i).vars == v2,:);

            uGroup = uGroup + diag(pairwise(map1, map2))';
        end
        
        uc(i,:) = uGroup;
    end
    
    % set the coarse pairwise terms (Eq. (5))
    adjMat = zeros(numel(vg));    % coarse scale's adjacency matrix
    pc = zeros(G.numLabels, G.numLabels, nnz(~vbTouch));
    adjc = zeros(nnz(~vbTouch),2);
    nEdgeCounter = 1;
    vNoTouch = find(~vbTouch);
    for i = 1 : numel(vNoTouch)

        iEdge = vNoTouch(i);
              
        % coarse representatives of edge's variables
        v1 = G.adj(iEdge,1);
        v2 = G.adj(iEdge,2);
        v1c = mapFineToCoarse(v1);
        v2c = mapFineToCoarse(v2);
        map1 = mapInterpolation(v1,:);
        map2 = mapInterpolation(v2,:);
%         map1 = vg(v1c).map(vg(v1c).vars == v1,:);
%         map2 = vg(v2c).map(vg(v2c).vars == v2,:);

        % the pairwise term, considering the interpolation
        pairwise = G.p(:,:,iEdge);
        pairwise = pairwise(map1, map2);
        
        % check if the edge defines a self-loop
        if (v1c ~= v2c)

            % check if (v1c,v2c) is already defined on the coarse scale
            if (~adjMat(v1c,v2c))
                % append as a new edge
                
                % track changes to adjacency
                adjMat(v1c,v2c) = nEdgeCounter;     % pairwise stored as (v1c,v2c)
                adjMat(v2c,v1c) = -nEdgeCounter;	% flag that pairwise is transposed
                
                % update adjacency and pairwise
                adjc(nEdgeCounter,:) = [v1c, v2c];
                pc(:,:,nEdgeCounter) = pairwise;

                nEdgeCounter = nEdgeCounter + 1;
            else
                % add to an existing edge
                
                idx = adjMat(v1c,v2c);
                if (idx < 0)
                    pairwise = pairwise';
                end
                pc(:,:,abs(idx)) = pc(:,:,abs(idx)) + pairwise;
            end
            
        else
            % self-loop, add the edge to respective coarse unary term
            % this is the remainder of the second term in Eq. (4)
            
            uc(v1c,:) = uc(v1c,:) + diag(pairwise)';
        end
    end
    
    % resize to clear reserved space
    pc(:,:,nEdgeCounter:end) = [];
    adjc(nEdgeCounter:end,:) = [];
    
    % the coarse graph
    Gc.u = uc;
    Gc.p = pc;
    Gc.adj = adjc;
    Gc.numLabels = G.numLabels;
end
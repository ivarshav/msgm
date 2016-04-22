function [Gcond, xcond] = msgmConditionalDist(G, x, vb)
% msgmConditionalDist(G, x, vb) compute the condional distribution for
% variables in G given the labels of variables for which (vb == true)
%
% terminology
%   - conditioned variables: variables that are conditioned upon
%   - conditional graph: graph *after* conditioning upon the cond. vars.
%
% note that the energy resulting from the conditioned variables is a
% constant that does not affect the optimization of the conditional graph,
% and hence its calculation is dropped.
%

    % maps: (var inds in cond. graph) <-> (var inds in orig. graph)
    mapCondToOrig = find(~vb);
    mapOrigToCond = zeros(1, size(G.u, 1));
    mapOrigToCond(mapCondToOrig) = 1 : numel(mapCondToOrig);

    % account for the unary terms in the conditional graph
    ucond = G.u(mapCondToOrig,:);   

    % account for pairwise terms
    % some edges are fully contained, and some edges
    % have only one endpoint in the cond. graph
    adjcondInds = [];
    for i = 1 : size(G.p, 3)

        v1 = G.adj(i,1);
        v2 = G.adj(i,2);

        if (~vb(v1) && ~vb(v2))
            % both variables are in the conditional graph

            % append the i-th edge
            adjcondInds = cat(1, adjcondInds, i);

        elseif (~vb(v1))
            % only the first variable is in the conditional graph

            pairwise = G.p(:,:,i);
            v1cond = mapOrigToCond(v1);
            ucond(v1cond,:) = ucond(v1cond,:) + pairwise(:,x(v2))';

        elseif (~vb(v2))
            % only the second variable is in the conditional graph

            pairwise = G.p(:,:,i);
            v2cond = mapOrigToCond(v2);
            ucond(v2cond,:) = ucond(v2cond,:) + pairwise(x(v1),:);

        end

    end
    
    % collect all components
    Gcond.u = ucond;
    Gcond.p = G.p(:,:,adjcondInds);
    Gcond.adj = mapOrigToCond(G.adj(adjcondInds,:));
    Gcond.numLabels = G.numLabels;
    Gcond.bProcessed = true;
    
    % initial guess for the conditional graph
    xcond = msgmInherit(x, find(~vb));
    
end





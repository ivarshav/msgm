function x = msgmInterpolate(G, vg, xc, param)
% msgmInterpolate(G, vg, xc, mapFineToCoarse, mapInterpoation)
% given a labeling 'xc' of the coarse scale, interpolate the assignment
% according to the interpolation rule 'mapInterpolation'
%
    
    x = zeros(size(G.u, 1), 1);

    % interpolate the i-th group
    for i = 1 : numel(vg)

        x(vg(i).vars) = vg(i).map(:,xc(i));
    end
    
    if (param.bSoftInterpolation)
        % 'soft' interpolation:
        % fix labels of seed variables and optimize labels of all the rest

        % condition upon the seed variables
        vb = false(size(G.u, 1), 1);
        vb([vg.seed]) = true;
        [Gcond, ~] = msgmConditionalDist(G, x, vb);
        
        % optimize labels of the conditional graph 
        [~, xcond] = min(Gcond.u, [], 2);
		xASSERT = x;
		xASSERT(~vb) = xcond;
        xcond = msgmOptimizeScale(Gcond, xcond, param);
        
        % map the optimized labels to the original variables
        x(~vb) = xcond;
		
		assert(msgmEnergy(G, x) <= msgmEnergy(G, xASSERT));
    
    end
end
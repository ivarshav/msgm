function x = msgmVcycle(G, x, param)
% msgmVcycle(G, x, param) optimization & coarsening followed by
% interpolation & optimization 
%
% input:
%
%   G   -   graphical model (see msgm())
%
%   x   -   an initial labeling assignment to the variables (column vector,
%           empty if an initial guess is not available)
%
%   param  -   set of parameters, see msgmParams().
%
%
% output:
%
%   x   -   a labeling assignment of the variables
%


    % check stopping condition
    if (size(G.u,1) <= param.numMinVars)
        
        x = msgmOptimizeScale(G, x, param);       
        return;
    end

    
    % reparameterization of the energy potentials
    G = msgmReparam(G);
	
    
    % run inference on the current scale
    if (any(x))
    
        x = msgmOptimizeScale(G, x, param);
    end
  

    % coarsen the graph
    [Gc, xc, vg] = msgmCoarsening(G, param, x);

    
    % recursive call
    xc = msgmVcycle(Gc, xc, param);

    
    % interpolate solution
    x = msgmInterpolate(G, vg, xc, param);
    

    % run inference on the current scale
    x = msgmOptimizeScale(G, x, param);


end
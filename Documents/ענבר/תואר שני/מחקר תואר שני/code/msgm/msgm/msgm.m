function [x, e, t] = msgm(G, x, param)
% msgm(G, x, param) a multiscale framework for computing approximate
% minimum energy assignments of pairwise MRFs
% 
% Main entrypoint to multiscale optimization of graphical models.
% See msgmDemo() for an example.
% 
%
% input:
%
%   G   -   graphical model
%           G.u     [num variables x num labels]            unary term
%           G.adj   [num edges x 2]                         adjacency relations
%           G.p     [num labels x num labels x num edges]   pairwise terms
%           G.numLabels  (for internal use) number of labels
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

    tStart = tic;

    % initialize internal values
    G.numLabels = size(G.u, 2);
    
    if (isempty(param))
        % initialize default params
        
        param = msgmParams();
    end

    % V-cycles
    for i = 1 : param.numVcycles

        x = msgmVcycle(G, x, param);
    end

	% last optimization
    x = msgmOptimizeScale(G, x, param);

    % energy & time
    t = toc(tStart);
    e = msgmEnergy(G, x);

end
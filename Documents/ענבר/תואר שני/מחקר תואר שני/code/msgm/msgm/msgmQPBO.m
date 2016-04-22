function x = msgmQPBO(G, x)            
% msgmQPBO(G, x, param) wrapper for QPBO
%

    % prepare pairwise potentials for QPBO mex
    p = cat(1, G.adj(:,1)', ...
        G.adj(:,2)', ...
        squeeze(G.p(1,1,:))', ...
        squeeze(G.p(1,2,:))', ...
        squeeze(G.p(2,1,:))', ...
        squeeze(G.p(2,2,:))');
    
    x = QPBO_wrapper_mex(G.u', p, int32(x == 2), 'i');
    
    % map 0,1 labels to 1,2 labels
    x = x + 1;
end
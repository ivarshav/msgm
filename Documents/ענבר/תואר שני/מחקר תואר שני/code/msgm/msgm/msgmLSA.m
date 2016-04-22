function x = msgmLSA(G, x, param)
% msgmLSA(G, x, param) wrapper for LSA-TR
%

    % reparameterization
    p = cat(2, G.adj(:,1), ...
            G.adj(:,2), ...
            squeeze(G.p(1,1,:)), ...
            squeeze(G.p(1,2,:)), ...
            squeeze(G.p(2,1,:)), ...
            squeeze(G.p(2,2,:)));   
    [newUE, newSubPE, newSuperPE, newConst] = reparamEnergy(G.u', p);
    eng.UE = newUE;
    eng.subPE = newSubPE;
    eng.superPE = newSuperPE;
    eng.constTerm = newConst;

    % used for LSA-euc mode, for optimizing grids
    if (numel(x) == prod(param.imSz))

        x = reshape(x, param.imSz);
    end
    x = LSA_TR(eng, 0, x);
    x = x(:);

end


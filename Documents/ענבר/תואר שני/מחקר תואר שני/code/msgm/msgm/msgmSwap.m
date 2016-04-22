function x = msgmSwap(G, x, param)
% msgmSwap(G, x, param) ab-Swap for "hard" energies with general pairwise potentials
%
    % possible a/b combinations
    ab_combs = combnk(1:G.numLabels, 2);

    % do param.numSwapIterations iterations
    for i = 1 : param.numSwapIterations

        % random ordering of ab_combs
        perm = randperm(size(ab_combs, 1));
        for j = 1 : size(ab_combs, 1)

            a = ab_combs(perm(j),1);
            b = ab_combs(perm(j),2);
            x = swap(G, x, a, b, param);
        end
    end
end


%% Swap-QPBO

function xnew = swap(G, x, a, b, param)
% optimize by allowing only swap move a <-> b

    % find variables whose label is a or b
    vb_ab = (x == a) | (x == b);
    xnew = x;
    
    if (any(vb_ab))
        
        % fix labels of variables whose labels is not a,b
        [Gab, xab] = msgmConditionalDist(G, x, ~vb_ab);
        
        % only (a <--> b) swap move is allowed
        Gab.u = Gab.u(:,[a,b]);
        Gab.p = Gab.p([a,b],[a,b],:);
        xab = 1 * (xab == a) + 2 * (xab == b);
        
        % improve a,b initial guess
        switch (param.optimization)

            case 'QPBO'                  
                xab_opt = msgmQPBO(Gab, xab);

            case 'LSA'
                xab_opt = msgmLSA(Gab, xab, param);
        end   
        
        % insert optimized a,b move (xab_opt) into the initial labeling x
        xab_opt = a * (xab_opt == 1) + b * (xab_opt == 2);
        xnew(vb_ab) = xab_opt;

        % assert that energy does not increase
        assert(msgmEnergy(G, xnew) <= msgmEnergy(G, x));
    end
end

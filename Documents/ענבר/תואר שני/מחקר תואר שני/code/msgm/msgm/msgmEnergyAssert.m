function msgmEnergyAssert(G, x, Gc, xc)
% msgmEnergyAssert(G, x, Gc, xc) assertion of the consistency criterion

    if (any(x) && any(xc))
        
        assert(isequal(msgmEnergy(G, x), msgmEnergy(Gc, xc)));
    end
end
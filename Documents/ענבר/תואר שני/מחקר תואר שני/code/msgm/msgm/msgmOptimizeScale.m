function x = msgmOptimizeScale(G, x, param)
% msgmOptimizeScale(G, x, param) single-scale optimization of a graphical
% model 'G', given an initial guess 'x'

    if (isempty(x))
        % labels not initialized,
        % currently, move-making methods require initialization

        % initialize according to unary term
        [~, x] = min(G.u, [], 2);
    end

    if (~strcmp(param.optimization,'NONE'))

        % keep current labeling assignment to assert improvement
        x_ = x;

        % apply a single-scale optimization method
        if (G.numLabels == 2)
            % binary model, apply selected binary optimization method
            
            switch (param.optimization)
                
                case 'QPBO'                  
                    x = msgmQPBO(G, x);
                    
                case 'LSA'
                    x = msgmLSA(G, x, param);
    
            end            
        else
            % model with >2 variables, apply ab-swap move-making method
            
            x = msgmSwap(G, x, param);    
        end
        disp('optimizition energy:');
        disp(msgmEnergy(G, x));
        assert(msgmEnergy(G, x) <= msgmEnergy(G, x_));
    end
end
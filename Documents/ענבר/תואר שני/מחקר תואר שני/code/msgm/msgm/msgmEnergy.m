function e = msgmEnergy(G, x)
% msgmEnergy(G, x) compute the energy of a labeling assignment 'x' on
% graphical model 'G'
%

    % N-th significant digit for round-off
    N = 6;

    e = 0;
    
    % pairwise term
    if ~isempty(G.p)
        ind = sub2ind(size(G.p), ...
            x(G.adj(:,1)), ...  % label of v1
            x(G.adj(:,2)), ...  % label of v2
            (1 : size(G.p, 3))');	% index of edge
        e = sum(G.p(ind));
    end

    % unary term
    ind = sub2ind(size(G.u),...
        (1 : size(G.u, 1))', ...	% index of varialbe
        x);                         % label of variable
    e = e + sum(G.u(ind));
    
    % rouding to N-th significant digit
    e = round(e, N);
end

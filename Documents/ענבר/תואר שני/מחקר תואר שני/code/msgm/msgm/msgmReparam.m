function G = msgmReparam(G)
% msgmReparam(U,E,P) reparameterization of the energy potentials of G.
% This is necessary as a pre-processing step in the computation of
% the local conditional entropy score.
%
% A pairwise potential 'p12' is rewritten as
%
%       p12(x1,x2) = p12_(x1,x2) + u1_(x1) + u2_(x2)    s.t.
%       sum(p12_(x1,x2).^2) is minimized.
%
% The above optimization problem can be solved analytically, up to an
% additive constant which does not affect our framework.
% The solution is given by:
%
%       u1_(x1) = (1 / numLabels) * sum_x2(p12(x1,x2))
%
% The new pairwise becomes p12_, and we add the residuals u1_,u2_
% to the unary terms of the respective variables v1,v2,
% i.e. u1 = u1 + u1_, etc.

    % compute edge means u1_,u2_
    u1_ = mean(G.p, 2);
    u2_ = mean(G.p, 1);

    % update the pairwise:
    % p12_(x1,x2) = p12(x1,x2) - u1_(x1) - u2_(x2)
    G.p = bsxfun(@minus, G.p, u1_);
    G.p = bsxfun(@minus, G.p, u2_);

    % update the unary term u1 = u1 + u1_
    % ...for each label separately
    u1_ = reshape(u1_, G.numLabels, size(G.p, 3))';
    u2_ = reshape(u2_, G.numLabels, size(G.p, 3))';
    for i = 1 : G.numLabels

        G.u(:,i) = G.u(:,i) + ...
            accumarray([G.adj(:,1); G.adj(:,2)], ...
            [u1_(:,i); u2_(:,i)], [size(G.u, 1), 1]);
    end
end
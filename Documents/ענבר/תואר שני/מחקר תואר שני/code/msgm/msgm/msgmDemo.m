function [eMS, tMS, eSS, tSS] = msgmDemo()
% msgmDemo() toy problem demo for msmgm
%
% sample random potentials for a 4-conncected grid and optimize the energy
% 
% parameters
%   - GRID_SIZE : generate grid of size [GRID_SIZE x GRID_SIZE]
%   - N_LABELS  : size of the label set
%   - N_REPETITIONS : num of test-repetitions
%   - COUPLING : coupling parameter, values >1 correspond to "harder" models
%

    % parameters
    GRID_SIZE = 100;
    N_LABELS = 2;
    N_REPETITIONS = 10;
    COUPLING = 1;

    % generate the adjacency relations for [GRID_SIZE x GRID_SIZE] grid
    sz = [GRID_SIZE, GRID_SIZE];
    [ii, jj] = sparse_adj_matrix(sz, 1, 1);
    sel = ii<jj;
    G.adj = [ii(sel), jj(sel)];

    % initialize output data variables
    eMS = zeros(N_REPETITIONS,1);
    eSS = zeros(N_REPETITIONS,1);
    tMS = zeros(N_REPETITIONS,1);
    tSS = zeros(N_REPETITIONS,1);

    % do N_REPETITIONS iterations
    for i = 1 : N_REPETITIONS

        % fix random seed, for reproducibility
        rng(i);
        disp(strcat('iteration: ',num2str(i)));

        % generate the energy potentials by sampling from a random distribution
        G.u = round(randn(GRID_SIZE^2, N_LABELS), 1);
        G.p = COUPLING * round(randn(N_LABELS, N_LABELS, size(G.adj, 1)), 1);


        % set parameters for multiscale optimization
        param = msgmParams;
        param.imSz = [GRID_SIZE, GRID_SIZE];
        param.optimization = 'LSA';
        param.numSwapIterations = 1;
        param.bSoftInterpolation = false;
        param.numVcycles = 1;

        % multiscale
        [~, eMS(i), tMS(i)] = msgm(G, [], param);

        % single scale
        G.numLabels = size(G.u, 2);
        tSS_ = tic;
        x = msgmOptimizeScale(G, [], param);
        tSS(i) = toc(tSS_);
        eSS(i) = msgmEnergy(G, x);
    end
end


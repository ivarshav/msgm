function param = msgmParams()
% msgmParams() setting parameters for msgm()

    % V-cycle
    param.numVcycles = 1;               % num of V-cycles
    param.numMinVars = 10;              % num of variables for V-cycle stopping criterion
  
    % optimization
    param.optimization = 'LSA';         % 'QPBO' or 'LSA', 'NONE' for skipping
    param.numSwapIterations = 1;        % num of 'SWAP' iterations
    param.imSz = [];                    % used for LSA-euc mode, for optimizing grids    

    % coarsening
    param.numEntropyBins = 20;          % num of bins for conditional entorpy scores
    
    % interpolation
    param.bSoftInterpolation = true;    % boolean flag for using soft interpolation
end







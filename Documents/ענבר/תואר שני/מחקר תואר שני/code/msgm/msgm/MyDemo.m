function [eMS, tMS, eSS, tSS] = MyDemo()
%
% sample random potentials for a 4-connceted grid and optimize the energy
% 
% parameters
%   - GRID_SIZE : generate grid of size [GRID_SIZE x GRID_SIZE]
%   - N_LABELS  : size of the label set
%   - N_REPETITIONS : num of test-repetitions
%   - COUPLING : coupling parameter, values >1 correspond to "harder" models
%
    % experiment number
    j = 27; %not real exp 
    
    % parameters
    GRID_SIZE = 70;
    N_LABELS = 2;
    N_REPETITIONS = 1;
    COUPLING = 1;
    
    % generate the adjacency relations for [GRID_SIZE x GRID_SIZE] grid
    sz = [GRID_SIZE, GRID_SIZE];
    [ii, jj] = sparse_adj_matrix(sz, 1, 'inf');
    sel = ii<jj;
    G.adj = [ii(sel), jj(sel)];
    
    % set parameters for multiscale optimization
    param = msgmParams;
    param.imSz = [GRID_SIZE, GRID_SIZE];
    param.optimization = 'QPBO';
    param.numSwapIterations = 1;
    param.bSoftInterpolation = false;
    param.numVcycles = 15;
    
    % initialize output data variables
    eMS = zeros(N_REPETITIONS, param.numVcycles);
    eSS = zeros(N_REPETITIONS, param.numVcycles);
    tMS = zeros(N_REPETITIONS, 1);
    tSS = zeros(N_REPETITIONS, 1);
    
    header = {'exp', 'grid size', 'random initial assignment', 'num repetitions', 'optimization', 'numVcycles', 'bSoftInterpolation', ...
        'Variable grouping','num labels','coupling', 'eMS', 'tMS', 'eSS', 'tSS'};
    xlswrite('exp.xls', header);
    
    % random initial assignment
    % fix random seed, for reproducibility
    rng(j);
    y = ones(GRID_SIZE^2, 1) + round(rand(GRID_SIZE^2, 1));
    
    %graph
    fig = figure('Name', strcat('exp. ', num2str(j)));
    title(sprintf('Grid size: %d, optimization: %s', GRID_SIZE, param.optimization));
    xlabel('Vcycle');
    ylabel('Energy');
    graph_colors = num2cell(jet(N_REPETITIONS * 2), 2);
    legendInfo = '';
    leg_info = 1;
    
    % do N_REPETITIONS iterations
    for i = 1 : N_REPETITIONS

        % fix random seed, for reproducibility
        rng(i);
        disp(strcat('iteration: ',num2str(i)));
        
        % generate the energy potentials by sampling from a random distribution
        % unary potential of every variable is 0
        G.u = zeros(GRID_SIZE^2, N_LABELS);           
        % pairwise potentials between [-1, 1]
        G.p = COUPLING * round(-1 + (1 - (-1)) * rand(N_LABELS, N_LABELS, size(G.adj, 1)) * 10^1) / 10^1; 
        
        
        % multiscale
        [~, eMS(i, 1:param.numVcycles), tMS(i)] = msgm(G, y, param);

        % single scale
        G.numLabels = size(G.u, 2); 
        x = y;
        tSS_ = tic;
        for k = 1: param.numVcycles
            x = msgmOptimizeScale(G, x, param);
            eSS(i, k) = msgmEnergy(G, x);
        end
        tSS(i) = toc(tSS_);
        eSS(i) = msgmEnergy(G, x);
        hold on
        plot([1:param.numVcycles], eMS(i, :), '-o', 'Color', cell2mat(graph_colors(i))); %, 'LineWidth', 2);
        legendInfo{leg_info}= [strcat('eMS Rep', num2str(i))]; 
        plot([1:param.numVcycles], eSS(i,:),  '-.', 'Color', cell2mat(graph_colors(end - i + 1))); %, 'LineWidth', 2);
        legendInfo{leg_info + 1}= [strcat('eSS Rep', num2str(i))]; 
        leg_info = leg_info + 2;
    end
    % plot avg of all repetitions, if more than one
    if N_REPETITIONS > 1
        plot([1:param.numVcycles], mean(eMS, 1), '-*', 'Color', 'm', 'LineWidth', 2);
        legendInfo{N_REPETITIONS * 2 + 1}= ['avg. eMS ']; 
        plot([1:param.numVcycles], mean(eSS, 1), '-*', 'Color', [1,0.4,0.6], 'LineWidth', 2);
        legendInfo{N_REPETITIONS * 2 + 2}= ['avg. eSS ']; 
    end
    % add legend and more info to plot
    l = legend(legendInfo);
    set(l, 'Location', 'bestoutside'); 
    descr = {strcat('numVcycles: ', num2str(param.numVcycles));
        strcat('numLabels: ', num2str(N_LABELS));
        strcat('Coupling: ', num2str(COUPLING));
        strcat('Variable grouping: ', '8-connected');
        };
    ax1 = axes('Position',[0 0 1 1],'Visible','off');
    axes(ax1) % sets ax1 to current axes
    text(0.7,0.35,descr);
    hold off
    
    print(fig, strcat('results/exp', num2str(j)), '-djpeg');

    xlswrite('exp.xls', ...
        {j, GRID_SIZE, 'too long', N_REPETITIONS, param.optimization, param.numVcycles, param.bSoftInterpolation, '8-connected', N_LABELS, COUPLING, ReprVector(eMS(:, param.numVcycles)), ReprVector(tMS),... 
        ReprVector(eSS), ReprVector(tSS);}, 1, sprintf('A%d' ,(j+1)));
    
    
end

function [x_str] = ReprVector(x)
%
% parameters
%   - x  : vector to print
%
    x_str = '[';
    for k = 1 : size(x, 1)
        x_str = strcat(x_str, num2str(x(k)), ','); 
    end
    x_str(end) = ']';
end

% function [] = PrintVector(fileID, x)
% %
% % parameters
% %   - fileID : file descriptor (open for writing)
% %   - x  : vector to print
% %
%     fprintf(fileID, '[');
%     for k = 1 : size(x, 1)
%         fprintf(fileID,'%2.4f ', x(k));
%     end
%     fprintf(fileID, ']\r\n');
% end
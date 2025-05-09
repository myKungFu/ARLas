function [x0,k,ci,cix,rmse] = audiometer_piFit(stimLvl,resp,xmin,xmax)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [x0,k,ci,cix,rmse] = audiometer_piFit(stimLvl,resp);
%
% Author: Shawn Goodman
% Date: July 28, 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %xmin = 0;
    %xmax = 100;
    
    % pin down the ends
    frontPadx = ones(10,1) * xmin;
    frontPady = ones(10,1) * 0;
    backPadx = ones(10,1) * xmax;
    backPady = ones(10,1) * 1;
    stimLvl = [frontPadx;stimLvl;backPadx];
    resp = [frontPady;resp;backPady];

   [fitresult,gof] = createFit(stimLvl,resp);
    x0 = fitresult.x0;
    k = fitresult.k;
    % calculate Simultaneous Functional Bounds':
    step = 0.1;
    cix = (xmin:step:xmax)'; % vector of presentation levels
    ci = predint(fitresult,cix,0.95,'functional','off');
    % see the following link for more info:
    % https://www.mathworks.com/help/curvefit/confidence-and-prediction-bounds.html
    rmse = gof.rmse; % root mean square error
end

% Internal Functions ------------------------------------------------------
function [fitresult,gof] = createFit(x,y)
    [xData, yData] = prepareCurveData( x, y );
    % Set up fittype and options.
    ft = fittype( '1./(1+exp(-(k*(x-x0))))', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Robust = 'Bisquare';
    opts.StartPoint = [0.0971317812358475 0.823457828327293];

    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );

%     % Plot fit with data.
%     figure(10);
%     h = plot( fitresult, xData, yData );
%     legend( h, 'y vs. x', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
%     % Label axes
%     xlabel( 'x', 'Interpreter', 'none' );
%     ylabel( 'y', 'Interpreter', 'none' );
%     grid on
%     xlim([-10 100])
end


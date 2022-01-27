function [fitresult, gof] = gaussN(lag1, acf1,npeaks)
% by Xiaoyu Shi @ Bo Huang Lab, UCSF

%CREATEFIT(LAG1,ACF1)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : lag1
%      Y Output: acf1
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 05-Jan-2016 14:39:12


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( lag1, acf1 );
a = 'gauss';
b = num2str(npeaks);
ftype =[a b];
% Set up fittype and options.
ft = fittype( ftype );

opts = fitoptions( ft );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0];
opts.StartPoint = [1.57894736842105 1.46031158962866 0.0629411730092827 1.28333333333333 -1.58318570396634 0.054292145424026 0.841004184100418 -0.439511255130954 0.0410337534521606 0.821875 -1.20511150600423 0.0699431830273753 0.75 -0.0708889121178959 0.0488231359433184 0.741935483870929 -0.817585453093066 0.0817507249200038 0.666666666666667 1.05388182681939 0.0567584826701476 0.658385093167702 0.297733430895163 0.0258126833816592];
opts.Upper = [Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf];
opts.Normalize = 'on';

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'untitled fit 1' );
h = plot( fitresult, xData, yData );
legend( h, 'acf1 vs. lag1', 'untitled fit 1', 'Location', 'NorthEast' );
% Label axes
xlabel( 'lag1' );
ylabel( 'acf1' );
grid on


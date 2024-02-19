function [distributionType, secondMean] = histogramGaussFit(fName, x, y, targetX, type2DistFactor)
%     figure;
%     normalizedRasGtpEverySecond = load(fName).normalizedRasGtpEverySecond;
%     h = histogram(normalizedRasGtpEverySecond(end,:),'BinWidth',0.05,'Normalization','probability', 'FaceColor', '#ffc0cb', 'EdgeColor', 'white');


    scatter(x,y, 30, 'filled');
%     txt = ['p=',num2str(p,3)];
%     text(0.75,0.45,txt,'FontSize',16,'FontWeight','bold');
%     xlim ([0 1])
%     ylim ([0 max(ceil(max(y)*10)/10, .5)])
%     ylim ([0 .7])
%     set(gca, 'CLim',[0 0.2],'FontSize',12,'FontWeight','bold','LineWidth',2);
%     title(sprintf('Histogram at t=%.4f', t))
%     xlabel('RasGTP/Ras')
%     ylabel('Frequency')
    
    workspace;  % Make sure the workspace panel is showing.
    format long g;
    format compact;
    fontSize = 20;

    % First specify how many Gaussians there will be.
    numGaussians = 2;

    % Initial Gaussian Parameters
    initialGuesses = [[0 .5 ]', [.068 .1]'];
    % Add a little noise so that our first guess is not dead on accurate.
    startingGuesses = reshape(initialGuesses', 1, []);

    global c NumTrials TrialError
    % 	warning off

    % Initializations
    NumTrials = 0;  % Track trials
    TrialError = 0; % Track errors

    
    
    tFit = reshape(x, 1, []);
    y = reshape(y, 1, []);

    %-------------------------------------------------------------------------------------------------------------------------------------------
    % Perform an iterative fit using the FMINSEARCH function to optimize the height, width and center of the multiple Gaussians.
    options = optimset('TolX', 1e-4, 'MaxFunEvals', 10^12);  % Determines how close the model must fit the data
    % First, set some options for fminsearch().
    options.TolFun = 1e-4;
    options.TolX = 1e-4;
    options.MaxIter = 100000;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % HEAVY LIFTING DONE RIGHT HERE:
    % Run optimization
    [parameter, fval, flag, output] = fminsearch(@(lambda)(fitgauss(lambda, tFit, y)), startingGuesses, options);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %----------------------------------------------------------------------------------------------------------------
    % Now plot results.
    yhat = PlotComponentCurves(x, y, tFit, c, parameter);
    % Compute the residuals between the actual y and the estimated y and put that into the graph's title.
%     meanResidual = mean(abs(y - yhat));
%     fprintf('The mean of the absolute value of the residuals is %f.\n', meanResidual);
%     caption = sprintf('Estimation of %d Gaussian Curves that will fit data.  Mean Residual = %f.', numGaussians, meanResidual);
%     title(caption, 'FontSize', fontSize, 'Interpreter', 'none');
%     drawnow;

    % Make table for the fitted, estimated results.
    % First make numGaussians row by 3 column matrix: Column 1 = amplitude, column 2 = mean, column 3 = width.
    % 	parameter % Print to command window.
    estimatedMuSigma = reshape(parameter, 2, [])';
    gaussianParameters = [c, estimatedMuSigma];
    % Now sort parameters in order of increasing mean
    gaussianParameters = sortrows(gaussianParameters, 2);
%     tActual % Display actual table in the command window.
    % Create table of the output parameters and display it below the actual, true parameters.
%     tEstimate = table((1:numGaussians)', c(:), estimatedMuSigma(:, 1), estimatedMuSigma(:, 2), 'VariableNames', {'Number', 'Amplitude', 'Mean', 'Width'})

    secondMean = gaussianParameters(2, 2);
    if (gaussianParameters(2, 2) - targetX)*targetX > type2DistFactor % g = 0.02/t + t
        distributionType = 2;
    else
        distributionType = 1;
    end
    
end

%=======================================================================================================================================================
function yhat = PlotComponentCurves(x, y, t, c, parameter)
try
	fontSize = 20;
	% Get the means and widths.
	means = parameter(1 : 2 : end);
	widths = parameter(2 : 2 : end);
	% Now plot results.
% 	hFig2 = figure;
% 	hFig2.Name = 'Fitted Component Curves';
	% 	plot(x, y, '--', 'LineWidth', 2)
	hold on;
	yhat = zeros(1, length(t));
	numGaussians = length(c);
% 	legendStrings = cell(numGaussians + 2, 1);
    
    colorArray = ['#00E3CB'; '#67349C'];
    
	for k = 1 : numGaussians
		% Get each component curve.
        xxx = linspace(0,1,100);
		thisEstimatedCurve = c(k) .* gaussian(xxx, means(k), widths(k));
%         thisEstimatedCurve = c(k) .* gaussian(t, means(k), widths(k));
		% Plot component curves.
		p = plot(xxx, thisEstimatedCurve,  'LineWidth', 2);
        p.Color = colorArray(k,:);
		hold on;
		% Overall curve estimate is the sum of the component curves.
% 		yhat = yhat + thisEstimatedCurve;
% 		legendStrings{k} = sprintf('Estimated Gaussian %d', k);
	end
	% Plot original summation curve, that is the actual curve.
% 	plot(x, y, 'r-', 'LineWidth', 1)
	% Plot estimated summation curve, that is the estimate of the curve.
% 	plot(x, yhat, 'k--', 'LineWidth', 2)
	grid on;
	xlabel('RasGTP/Ras', 'FontSize', fontSize)
	ylabel('Frequency', 'FontSize', fontSize)
% 	caption = sprintf('Estimation of %d Gaussian Curves that will fit data.', numGaussians);
% 	title(caption, 'FontSize', fontSize, 'Interpreter', 'none');
	grid on
% 	legendStrings{numGaussians+1} = sprintf('Actual original signal');
% 	legendStrings{numGaussians+2} = sprintf('Sum of all %d Gaussians', numGaussians);
% 	legend('RasGTP distribution', 'Values', 'Gaussian1', 'Gaussian2', 'FontSize', 10);
% 	xlim(sort([x(1) x(end)]));
	hFig2.WindowState = 'maximized';
	drawnow;
	
catch ME
	% Some error happened if you get here.
	callStackString = GetCallStack(ME);
	errorMessage = sprintf('Error in program %s.\nTraceback (most recent at top):\n%s\nError Message:\n%s', ...
		mfilename, callStackString, ME.message);
	WarnUser(errorMessage);
end
end % of PlotComponentCurves

%=======================================================================================================================================================
function theError = fitgauss(lambda, t, y)
% Fitting function for multiple overlapping Gaussians, with statements
% added (lines 18 and 19) to slow the progress and plot each step along the
% way, for educational purposes.
% Author: T. C. O'Haver, 2006

global c NumTrials TrialError
try
	
	A = zeros(length(t), round(length(lambda) / 2));
	for j = 1 : length(lambda) / 2
		A(:,j) = gaussian(t, lambda(2 * j - 1), lambda(2 * j))';
	end
	
	c = A \ y';
	z = A * c;
	theError = norm(z - y');
	
	% Penalty so that heights don't become negative.
	if sum(c < 0) > 0
		theError = theError + 1000000;
	end
	
	NumTrials = NumTrials + 1;
	TrialError(NumTrials) = theError;
catch ME
	% Some error happened if you get here.
	callStackString = GetCallStack(ME);
	errorMessage = sprintf('Error in program %s.\nTraceback (most recent at top):\n%s\nError Message:\n%s', ...
		mfilename, callStackString, ME.message);
	WarnUser(errorMessage);
end
end % of fitgauss()

%=======================================================================================================================================================
% Gets a string describing the call stack where each line is the filename, function name, and line number in that file.
% Sample usage
% try
% 	% Some code that might throw an error......
% catch ME
% 	callStackString = GetCallStack(ME);
% 	errorMessage = sprintf('Error in program %s.\nTraceback (most recent at top):\n%s\nError Message:\n%s', ...
% 		mfilename, callStackString, ME.message);
% 	WarnUser(errorMessage);
% end
function callStackString = GetCallStack(errorObject)
try
	theStack = errorObject.stack;
	callStackString = '';
	stackLength = length(theStack);
	% Get the date of the main, top level function:
	% 	d = dir(theStack(1).file);
	% 	fileDateTime = d.date(1:end-3);
	if stackLength <= 3
		% Some problem in the OpeningFcn
		% Only the first item is useful, so just alert on that.
		[folder, baseFileName, ext] = fileparts(theStack(1).file);
		baseFileName = sprintf('%s%s', baseFileName, ext);	% Tack on extension.
		callStackString = sprintf('%s in file %s, in the function %s, at line %d\n', callStackString, baseFileName, theStack(1).name, theStack(1).line);
	else
		% Got past the OpeningFcn and had a problem in some other function.
		for k = 1 : length(theStack)-3
			[folder, baseFileName, ext] = fileparts(theStack(k).file);
			baseFileName = sprintf('%s%s', baseFileName, ext);	% Tack on extension.
			callStackString = sprintf('%s in file %s, in the function %s, at line %d\n', callStackString, baseFileName, theStack(k).name, theStack(k).line);
		end
	end
catch ME
	errorMessage = sprintf('Error in program %s.\nTraceback (most recent at top):\nError Message:\n%s', ...
		mfilename, ME.message);
	WarnUser(errorMessage);
end
end % from callStackString

%==========================================================================================================================
% Pops up a warning message, and prints the error to the command window.
function WarnUser(warningMessage)
if nargin == 0
	return; % Bail out if they called it without any arguments.
end
try
	fprintf('%s\n', warningMessage);
	uiwait(warndlg(warningMessage));
	% Write the warning message to the log file
	folder = 'C:\Users\Public\Documents\MATLAB Settings';
	if ~exist(folder, 'dir')
		mkdir(folder);
	end
	fullFileName = fullfile(folder, 'Error Log.txt');
	fid = fopen(fullFileName, 'at');
	if fid >= 0
		fprintf(fid, '\nThe error below occurred on %s.\n%s\n', datestr(now), warningMessage);
		fprintf(fid, '-------------------------------------------------------------------------------\n');
		fclose(fid);
	end
catch ME
	message = sprintf('Error in WarnUser():\n%s', ME.message);
	fprintf('%s\n', message);
	uiwait(warndlg(message));
end
end % from WarnUser()

%=======================================================================================================================================================
function g = gaussian(x, peakPosition, width)
%  gaussian(x,pos,wid) = gaussian peak centered on pos, half-width=wid
%  x may be scalar, vector, or matrix, pos and wid both scalar
%  T. C. O'Haver, 1988
% Examples: gaussian([0 1 2],1,2) gives result [0.5000    1.0000    0.5000]
% plot(gaussian([1:100],50,20)) displays gaussian band centered at 50 with width 20.
g = exp(-((x - peakPosition) ./ (0.60056120439323 .* width)) .^ 2);
end % of gaussian()

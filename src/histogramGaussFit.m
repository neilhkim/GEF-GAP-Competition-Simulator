function [distributionType, higherGaussianMean] = histogramGaussFit(x, y, rateLawXVal, typeOffset)

workspace;  % Make sure the workspace panel is showing.
format long g;
format compact;

% Adjusted Initial Gaussian Parameters
initialGuesses = [[0.0 0.5]', [0.5 0.1]'];  % Slightly adjusted for a broader starting point

% Add a little noise so that our first guess is not on a potential singular point.
noiseToAdd = rand(size(initialGuesses)) * 0.01;  % Adding small noise
startingGuesses = reshape(initialGuesses + noiseToAdd', 1, []);

global c NumTrials TrialError

% Initializations
NumTrials = 0;  % Track trials
TrialError = 0; % Track errors

tFit = reshape(x, 1, []);
y = reshape(y, 1, []);

% Perform an iterative fit using the FMINSEARCH function to optimize the height, width, and center of the multiple Gaussians.
options = optimset('TolX', 1e-4, 'MaxFunEvals', 10^12);  % Determines how close the model must fit the data
% Set fminsearch() options.
options.TolFun = 1e-4;
options.TolX = 1e-4;
options.MaxIter = 100000;

% Run optimization (fitting)
warning('off', 'all')
[parameter, ~, ~, ~] = fminsearch(@(lambda)(fitgauss(lambda, tFit, y)), startingGuesses, options);
warning('on', 'all')

% Handle results.
plotComponentCurves(tFit, c, parameter);
estimatedMuSigma = reshape(parameter, 2, [])';
gaussianParameters = [c, estimatedMuSigma];
% Sort parameters in order of increasing mean
gaussianParameters = sortrows(gaussianParameters, 2);
higherGaussianMean = gaussianParameters(2, 2);
if higherGaussianMean > rateLawXVal + typeOffset
    distributionType = 2;
else
    distributionType = 1;
end
    
end


%=======================================================================================================================================================
function yhat = plotComponentCurves(t, c, parameter)
try
	% Get the means and widths.
	means = parameter(1 : 2 : end);
	widths = parameter(2 : 2 : end);
	% Now plot results.
	hold on;
	yhat = zeros(1, length(t));
	numGaussians = length(c);
    
    colorArray = ['#00E3CB'; '#67349C'];
    
	for k = 1 : numGaussians
		% Get each component curve.
        xxx = linspace(0,1,100);
		thisEstimatedCurve = c(k) .* gaussian(xxx, means(k), widths(k));
		% Plot component curves.
		p = plot(xxx, thisEstimatedCurve,  'LineWidth', 2);
        p.Color = colorArray(k,:);
		hold on;
	end
	drawnow;
	
catch ME
	% Some error happened if you get here.
	callStackString = GetCallStack(ME);
	errorMessage = sprintf('Error in program %s.\nTraceback (most recent at top):\n%s\nError Message:\n%s', ...
		mfilename, callStackString, ME.message);
	WarnUser(errorMessage);
end
end % of plotComponentCurves

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
function g = gaussian(x, peakPosition, width)
%  gaussian(x,pos,wid) = gaussian peak centered on pos, half-width=wid
%  x may be scalar, vector, or matrix, pos and wid both scalar
%  T. C. O'Haver, 1988
% Examples: gaussian([0 1 2],1,2) gives result [0.5000    1.0000    0.5000]
% plot(gaussian([1:100],50,20)) displays gaussian band centered at 50 with width 20.
g = exp(-((x - peakPosition) ./ (0.60056120439323 .* width)) .^ 2);
end % of gaussian()


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
% This function is adapted from 
% Mark Hayworth (2024). MAGIC - MATLAB Generic Imaging Component 
% (https://www.mathworks.com/matlabcentral/fileexchange/24224-magic-matlab-generic-imaging-component), 
% MATLAB Central File Exchange. Retrieved March 7, 2024.
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
		[~, baseFileName, ext] = fileparts(theStack(1).file);
		baseFileName = sprintf('%s%s', baseFileName, ext);	% Tack on extension.
		callStackString = sprintf('%s in file %s, in the function %s, at line %d\n', callStackString, baseFileName, theStack(1).name, theStack(1).line);
	else
		% Got past the OpeningFcn and had a problem in some other function.
		for k = 1 : length(theStack)-3
			[~, baseFileName, ext] = fileparts(theStack(k).file);
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
% This function is adapted from 
% Mark Hayworth (2024). MAGIC - MATLAB Generic Imaging Component 
% (https://www.mathworks.com/matlabcentral/fileexchange/24224-magic-matlab-generic-imaging-component), 
% MATLAB Central File Exchange. Retrieved March 7, 2024.
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

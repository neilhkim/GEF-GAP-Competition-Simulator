% Histogram Gaussian Fit for resulting x-values from the GEF-GAP competition simulation.
% This function applies Gaussian fitting to the histogram of simulated x-values derived from GEF-GAP competition simulations.
% It assesses the resulting distribution to determine whether it is unimodal or bimodal by comparing the position of the higher Gaussian mean
% to a predefined threshold. This threshold is determined by the sum of rateLawXVal and typeOffset.
% Neil H. Kim, 2024.
function [distributionType, higherGaussianMean] = histogramGaussFit(x, y, rateLawXVal, typeOffset)

workspace;  % Ensures that the MATLAB workspace panel is open for easy access to variables during debugging.
format long g; % Sets the numeric output format to 'long' for more precision in the display.
format compact; % Sets the display format to 'compact' to reduce the amount of vertical space used for outputs.

% Defines initial guesses for the parameters of the Gaussian functions. These parameters are slightly adjusted
% to provide a broader starting point for the optimization process, aiming to improve the chances of finding
% a global rather than a local optimum.
initialGuesses = [[0.0 0.5]', [0.5 0.1]'];

% Introduces a small amount of random noise to the initial parameter guesses. This helps prevent the optimization
% algorithm from getting stuck in singular points by starting the search from a slightly perturbed position.
noiseToAdd = rand(size(initialGuesses)) * 0.01;
startingGuesses = reshape(initialGuesses + noiseToAdd', 1, []);

global c NumTrials TrialError % Declares global variables for tracking the fitting process and any errors.

% Initializes the counters for the number of trials (attempts at fitting) and the cumulative error across these trials.
NumTrials = 0;
TrialError = 0;

% Reshapes the input data vectors 'x' and 'y' into row vectors to ensure compatibility with MATLAB's optimization functions.
tFit = reshape(x, 1, []);
y = reshape(y, 1, []);

% Configures the options for MATLAB's fminsearch optimization function. These options specify the tolerances for the solution
% (TolX, TolFun) and set a maximum number of function evaluations (MaxFunEvals) and iterations (MaxIter) to prevent infinite loops.
options = optimset('TolX', 1e-4, 'MaxFunEvals', 10^12, 'TolFun', 1e-4, 'MaxIter', 100000);

% Temporarily disables all MATLAB warnings to prevent them from interfering with the user interface. This is particularly
% useful for optimization processes that might generate non-critical warnings.
warning('off', 'all')
% Executes the optimization process by fitting the Gaussian model to the data. The fitgauss function is called with the current
% guesses for the Gaussian parameters and the data to fit.
[parameter, ~, ~, ~] = fminsearch(@(lambda)(fitgauss(lambda, tFit, y)), startingGuesses, options);
% Re-enables MATLAB warnings after the optimization is complete.
warning('on', 'all')

% Processes the results of the fitting. This includes plotting the fitted Gaussians over the data and determining the
% parameters of the fitted Gaussians (means, widths, and amplitudes).
plotComponentCurves(tFit, c, parameter);
estimatedMuSigma = reshape(parameter, 2, [])'; % Reshapes the fitted parameters for easier interpretation.
gaussianParameters = [c, estimatedMuSigma]; % Combines coefficients (amplitudes) with the means and widths.
% Sorts the Gaussian parameters by their mean values to identify the Gaussian with the highest mean.
gaussianParameters = sortrows(gaussianParameters, 2);
higherGaussianMean = gaussianParameters(2, 2); % The mean of the Gaussian with the higher mean value.

% Determines the type of distribution (unimodal or bimodal) based on the position of the higher Gaussian mean relative
% to the predefined threshold (rateLawXVal + typeOffset).
if higherGaussianMean > rateLawXVal + typeOffset
    distributionType = 2; % Indicates a bimodal distribution.
else
    distributionType = 1; % Indicates a unimodal distribution.
end
end % of histogramGaussFit()


%=======================================================================================================================================================
function yhat = plotComponentCurves(t, c, parameter)
try
    % Extracts the means and widths for Gaussian components from the parameter array.
    % Means are located at odd indices, and widths are at even indices.
    means = parameter(1 : 2 : end);
    widths = parameter(2 : 2 : end);

    % Prepares for plotting by holding the current figure.
    hold on;
    
    % Initializes a vector to store the combined Gaussian curve data. This implementation focuses on plotting,
    % so 'yhat' remains unused but is prepared for potential future use or data analysis.
    yhat = zeros(1, length(t));
    
    % Determines the number of Gaussian components based on the length of the coefficients array.
    numGaussians = length(c);

    % Defines a color array for distinguishing between different Gaussian components on the plot.
    colorArray = ['#00E3CB'; '#67349C'];

    for k = 1 : numGaussians
        % Generates a dense set of x-values for plotting each Gaussian component smoothly.
        xxx = linspace(0, 1, 100);
        
        % Calculates the y-values for the current Gaussian component based on its mean, width, and coefficient.
        thisEstimatedCurve = c(k) .* gaussian(xxx, means(k), widths(k));
        
        % Plots the current Gaussian component with a specified line width and color.
        p = plot(xxx, thisEstimatedCurve, 'LineWidth', 2);
        p.Color = colorArray(k,:);
    end
    
    % Forces MATLAB to update the figure window to display the newly plotted components.
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

% --------------------------------------------------------------------------------------------------------------------
% This script conducts a simulation of Ras signaling dynamics on an compartmentalized (corraled) in-vitro membrane, 
% focusing on the interactions between GEF (Sos), GAP, and Ras in various states (GDP-bound and GTP-bound). 
% It utilizes a stochastic approach to model the molecular interactions and changes over time within a corraled membrane. 
% Key aspects of the simulation include initializing parameters that define the biochemical reaction rates, setting up a simulation environment, 
% and dynamically updating the system using Gillespie method. 
% Neil H. Kim and Albert A. Lee, 2024.


close all; % Close all open figures
% clear all; % Uncomment to clear all variables from the workspace
addpath('../lib') % Add the 'lib' directory to the path for accessing external functions

% Simulation parameters
kCatGef = 1; % Set 1/kCatGEF as unit time
kCatGap = 2^-2; % Set kcat of GAP
kOnGtp2KOnGdp = 2^0; % Feedback strength: GEF's propensity to bind to RasGTP compared to RasGDP
kOffGdp = 2^-2; % GEF unbinding rate from RasGDP
kOffGtp2KOffGdp = 1; % Ratio of GEF's unbinding rates for RasGTP to RasGDP, set to 1 by default
kOffGtp = kOffGdp * kOffGtp2KOffGdp; % Calculate GEF's unbinding rate from RasGTP

rateLawXVal = 0.6; % x-value used to determine the binding rates by solving ODEs later
                   % x-value is defined as (# of RasGTP) / (# of RasGTP + # RasGDP)

% Simulation grid and entity settings
gridLen = 1; % Length of the grid
RAS_DENSITY = 1000; % Density of Ras molecules
N_CORRALS = 500; % Number of corrals
N_MIN_DURATION = 30; % Minimum duration for the simulation
LOWEST_AVG_X_REQD_FOR_DIP_TEST = 0.1; % Minimum average X required for the DIP test
N_TRACKED_VARIABLES = 5; % Number of variables tracked in the simulation
TYPE_OFFSET = 0; % Offset used to determine the distribution type in the DIP test

% Create output directory if it doesn't exist
folderIndex = 1; % Start with index 1 for naming output folders
while exist(sprintf('.\\output%d', folderIndex), 'dir')
    folderIndex = folderIndex + 1; % Increment index if folder already exists
end
mkdir(sprintf('.\\output%d', folderIndex)); % Create the new output directory

% Show initial parameter values
fprintf('kCatGap = %f \n', kCatGap);
fprintf('kOnGtp2KOnGdp = %f \n',kOnGtp2KOnGdp);
fprintf('kOffGdp = %f \n',kOffGdp);
fprintf('kOffGtp = %f \n',kOffGtp);

nRasTotal = RAS_DENSITY * gridLen^2; % Total Ras molecules
nTargetRasGtp = rateLawXVal * nRasTotal; % Target RasGTP based on rate law value

% Solve for initial conditions using analytical solution
[SosRasGdp,SosRasGtp,RasGdp,kOnGdp] = getKOnGtp(kCatGap, kCatGef, kOnGtp2KOnGdp, kOffGdp, kOffGtp, nTargetRasGtp, nRasTotal); % Solves for binding rate and concentrations. The binding rate kOnGdp is k_on (conventional on-rate) times the supposed concentration of GEF.
str = sprintf('\n  Analytical solution for x=%f\n-> [kOnGdp = %f, kOnGtp = %f, RasGdp = %f, RasGtp = %f, SosRasGdp = %f, SosRasGtp = %f]\n', rateLawXVal, kOnGdp, kOnGtp2KOnGdp * kOnGdp, RasGdp, nRasTotal-RasGdp-SosRasGdp-SosRasGtp, SosRasGdp, SosRasGtp);
disp(str); % Display the solution

if length(kOnGdp) > 1
    disp('Error in finding deterministic solution.'); % Check for error in solution
    return;
end

kOnGtp = kOnGtp2KOnGdp * kOnGdp; % Calculate binding rate of GEF (Sos) to Ras-GTP

% Simulation duration setup
duration = 1 + nTargetRasGtp / 2 / (kOffGtp + nTargetRasGtp * kOnGdp) * kOffGdp;
duration = duration * 2; % Adjust duration based on simulation needs
duration = round(duration); % Round to nearest whole number
if duration > 20000
    duration = 20000; % Cap the duration to avoid excessively long simulations
end

% Initialize random number generator for simulation
rng('shuffle');

nTrackedVars = 5; % tracked variables: (t,nRasGtp,nRasGdp,nSosRasGtp,nSosRasGdp)

% Prepare matrix to track time and variables across all corrals. Give me much
% time points (10 * duration (in unit times))
timeTraceAllCorrals = zeros(10 * double(duration), nTrackedVars * N_CORRALS);

% Simulation loop for each corral, iterating through the total number of corrals
for iCorral = 1:N_CORRALS
    % Periodically output progress for every half of the total corrals processed
    if mod(iCorral, N_CORRALS/2) == 0
        str = sprintf(' Results: iCorral = %d out of %d', iCorral, N_CORRALS);
        disp(str) % Display the progress
    end

    % Setting up initial conditions for the simulation within a corral
    t = 0; % Initialize time
    x_initial = 0; % Initial x-value ((# of RasGTP)/(# of RasGDP + # of RasGTP), assumed to be 0 for start)
    nRasGtp = nRasTotal * x_initial; % Calculate initial number of RasGTP based on total Ras and initial concentration
    nRasGdp = nRasTotal * (1-x_initial); % Remaining Ras are RasGDP
    nSosRasGtp = 0; % Initial number of Sos bound to RasGTP
    nSosRasGdp = 0; % Initial number of Sos bound to RasGDP

    % Record the initial state for this corral
    record = [t, nRasGtp, nRasGdp, nSosRasGtp, nSosRasGdp];

    % Initialize record index for storing simulation results
    recIdx = 0;

    % Store the initial state record into the timeTraceAllCorrals array
    timeTraceAllCorrals(recIdx+1, ((iCorral-1)*nTrackedVars + 1) : iCorral*nTrackedVars) = record;

    % Use the Gillespie algorithm for stochastic simulation, iterating until the simulation duration is reached
    while (t < duration) || (t < N_MIN_DURATION)
        r1 = kOnGtp * nRasGtp; % Rate of SOS (GEF) binding to RasGtp
        r2 = kOnGdp * nRasGdp; % Rate of SOS binding to RasGdp
        r3 = kOffGtp * nSosRasGtp; % Rate of SOS unbinding from RasGtp
        r4 = kOffGdp * nSosRasGdp; % Rate of SOS unbinding from RasGdp
        r5 = kCatGef * (nSosRasGtp + nSosRasGdp) / gridLen^2 * nRasGdp / gridLen^2; % Rate of SOS converting RasGdp to RasGtp
        r6 = kCatGap / gridLen^2 * nRasGtp / gridLen^2; % Rate of GAP converting RasGtp to RasGdp

        % Total rate of all reactions
        rT = r1 + r2 + r3 + r4 + r5 + r6;

        % Draw random time interval based on total rate, simulating the time to the next reaction
        drawT = rand;
        tau = 1/rT * log(1/drawT);

        % Advance simulation time by the drawn interval
        t = t + tau;
        t = double(t);

        % Determine the reaction that occurs based on their probabilities
        p1 = r1 / rT;
        p2 = r2 / rT;
        p3 = r3 / rT;
        p4 = r4 / rT;
        p5 = r5 / rT;
        drawReaction = rand; % Draw a random number to decide which reaction occurs

        % Update the system state based on the reaction that occurred
        if drawReaction < p1
            nRasGtp = nRasGtp - 1;
            nSosRasGtp = nSosRasGtp + 1;
        elseif drawReaction < (p1+p2)
            nSosRasGdp = nSosRasGdp + 1;
            nRasGdp = nRasGdp - 1;
        elseif drawReaction < (p1+p2+p3)
            nSosRasGtp = nSosRasGtp - 1;
            nRasGtp = nRasGtp + 1;
        elseif drawReaction < (p1+p2+p3+p4)
            nSosRasGdp = nSosRasGdp - 1;
            nRasGdp = nRasGdp + 1;
        elseif drawReaction < (p1+p2+p3+p4+p5)
            nRasGtp=nRasGtp+1;
            nRasGdp=nRasGdp-1;
        else
            nRasGtp=nRasGtp-1;
            nRasGdp=nRasGdp+1;
        end

        % Record the updated state at specified time intervals
        if floor(t*10) > recIdx % Condition to record data at fixed intervals (every 0.1 unit time)
            recIdx = recIdx + 1; % Increment the record index
            % Update the timeTraceAllCorrals array with the new state
            timeTraceAllCorrals(recIdx+1, (iCorral-1)*nTrackedVars + 1 : iCorral * nTrackedVars) = [t, nRasGtp, nRasGdp, nSosRasGtp, nSosRasGdp];
        end
    end
end % Full runtime for a corral end

% Extracting specific time trace data from the comprehensive simulation results for further analysis
rasGtpTimeTrace = timeTraceAllCorrals(:,2:N_TRACKED_VARIABLES:N_CORRALS*N_TRACKED_VARIABLES); % Extract RasGTP time trace data
sosGtpTimeTrace = timeTraceAllCorrals(:,4:N_TRACKED_VARIABLES:N_CORRALS*N_TRACKED_VARIABLES); % Extract Sos-RasGTP complex time trace data
sosGdpTimeTrace = timeTraceAllCorrals(:,5:N_TRACKED_VARIABLES:N_CORRALS*N_TRACKED_VARIABLES); % Extract Sos-RasGDP complex time trace data
recordedTimeArray = timeTraceAllCorrals(:,1:N_TRACKED_VARIABLES:N_CORRALS*N_TRACKED_VARIABLES); % Extract time points of the recordings
dataLength = max(duration, N_MIN_DURATION); % Determine the length of data based on the maximum of either the simulation duration or the minimum duration


% Reorganizing extracted time trace data for analysis at each unit time
nUnitTimes = 0:1:dataLength; % Create an array representing each unit time interval up to the data length
rasGtpEveryUnitTime = reorganizeForEveryUnitTime(recordedTimeArray, rasGtpTimeTrace, dataLength, N_CORRALS); % Reorganize RasGTP data to align with unit times
sosGtpEveryUnitTime = reorganizeForEveryUnitTime(recordedTimeArray, sosGtpTimeTrace, dataLength, N_CORRALS); % Reorganize Sos-RasGTP data to align with unit times
sosGdpEveryUnitTime = reorganizeForEveryUnitTime(recordedTimeArray, sosGdpTimeTrace, dataLength, N_CORRALS); % Reorganize Sos-RasGDP data to align with unit times


% Calculating mean values of SOS-GTP and SOS-GDP across all corrals for each unit time
meanOfSosGtpEveryUnitTime = mean(sosGtpEveryUnitTime,2); % Calculate mean SOS-GTP for each unit time
meanOfSosGdpEveryUnitTime = mean(sosGdpEveryUnitTime,2); % Calculate mean SOS-GDP for each unit time

% Calculating and recording the mean of normalized RasGTP 
normalizedRasGtpEveryUnitTime = rasGtpEveryUnitTime/nRasTotal; % Normalize RasGTP counts by total Ras to get a proportion (== x-value)
meanOfNormalizedRasGtpEveryUnitTime = mean(normalizedRasGtpEveryUnitTime,2); % Calculate the mean normalized RasGTP for each unit time

% Saving processed data to .mat files for further analysis or visualization

% Save norm RasGTP every unit time as .mat file
if ~exist(sprintf('.\\output%d\\variables', folderIndex), 'dir')
    mkdir(sprintf('.\\output%d\\variables', folderIndex));
end
filename = sprintf('.\\output%d\\variables\\normalizedRasGtpEveryUnitTime.mat', folderIndex);
save(filename,'normalizedRasGtpEveryUnitTime');

% Save the last timepoint's norm RasGTPs as .mat file
lastNormRasGtpArray = normalizedRasGtpEveryUnitTime(end,:);
filename = sprintf('.\\output%d\\variables\\endPointNormXDistribution.mat', folderIndex);
save(filename,'lastNormRasGtpArray');

% Performing a DIP test to analyze the distribution of the last timepoint's normalized RasGTP values
if sum(lastNormRasGtpArray) / length(lastNormRasGtpArray) > LOWEST_AVG_X_REQD_FOR_DIP_TEST % Check if the average X-value is sufficient for the DIP test
    bootStrapSampleSize = 100; % Set the bootstrap sample size for the DIP test
    [~, p] = HartigansDipSignifTest(normalizedRasGtpEveryUnitTime(end,:), bootStrapSampleSize); % Perform the DIP test and store the p-value, where low p-value indicates a non-unimodal distribution
end

% Plot the result
figure;
plot(nUnitTimes,meanOfNormalizedRasGtpEveryUnitTime,'LineWidth',2, 'Color', '#004C99');
set(gca, 'CLim',[0 0.2],'FontSize',12,'FontWeight','bold','LineWidth',2);
xlim ([0 size(nUnitTimes,2)])
ylim ([0 1])
title('The mean of normalized RasGTP every second')
xlabel('Time (s)')
ylabel('RasGTP/Ras')

figure;
h = histogram(lastNormRasGtpArray,'BinWidth',0.05,'Normalization','probability','FaceColor', '#ffc0cb', 'EdgeColor', 'white');
histXArray = h.BinEdges(1:end-1);
histYArray = h.Values;
xlim ([0 1])
ylim ([0 max(histYArray)*1.1]);
set(gca, 'CLim',[0 0.2],'FontSize',12,'FontWeight','bold','LineWidth',2);
title(sprintf('Histogram at t=%.2f', t))
xlabel('RasGTP/Ras')
ylabel('Frequency')
hold on;
higherGaussMean = nan;

% Make sure that the average X-value is bigger than the required for the DIP test (threshold set arbitrarly).
if sum(lastNormRasGtpArray) / length(lastNormRasGtpArray) > LOWEST_AVG_X_REQD_FOR_DIP_TEST
    plot(rateLawXVal*[1 1], [0 1],'--', 'LineWidth', 2, 'color', '#D90000');
    txt = ['p = ',num2str(p,3),','];
    text(0.1,max(histYArray),txt,'FontSize',12,'FontWeight','bold','HorizontalAlignment', 'left','VerticalAlignment', 'top');
    % If the p-value is less than 0.05, further
    % characterize the bimodal distribution.
    if p < 0.05
        % fit two gaussians and retrieve distribution type and the higher gaussian mean of the two
        % If the higher gauss mean is larger than rate law x-value + the
        % type determining offset, it is type 2. (else type 1)
        [distribType, higherGaussMean] = histogramGaussFit(histXArray,histYArray,rateLawXVal,TYPE_OFFSET);
        plot([1 1]*(rateLawXVal + TYPE_OFFSET), [0,1], '-.', 'LineWidth', 2, 'color', '#e03fd8');
        plot([1 1]*higherGaussMean, [0,1], '-', 'LineWidth', 2, 'color', '#67349C');
        text(higherGaussMean*1.03, max(histYArray)*0.5, sprintf('%.3f',higherGaussMean),'FontSize',12,'FontWeight','bold','HorizontalAlignment', 'left','VerticalAlignment', 'bottom', 'color', '#67349C');
        text((rateLawXVal + TYPE_OFFSET)*1.03, max(histYArray)*0.2, 'rate law x-val','FontSize',12,'FontWeight','bold','HorizontalAlignment', 'left','VerticalAlignment', 'bottom', 'color', '#e03fd8');
    else
        % distribType == 0 means Unimodal distribution
        distribType = 0;
        text(0.5, max(histYArray)*0.2, 'Unimodal','FontSize',12,'FontWeight','bold','HorizontalAlignment', 'center','VerticalAlignment', 'bottom');
    end
    txt2 = ['Stochastic Bistability Type = ', num2str(distribType)];
    text(0.5,max(histYArray),txt2,'FontSize',12,'FontWeight','bold','HorizontalAlignment', 'center','VerticalAlignment', 'top');
else
    txt2 = 'x-value too low to draw conclusion';
    text(0.5,max(histYArray),txt2,'FontSize',12,'FontWeight','bold','HorizontalAlignment', 'center','VerticalAlignment', 'top');
end
hold off;

if ~exist(sprintf('.\\output%d\\png', folderIndex), 'dir')
    mkdir(sprintf('.\\output%d\\png', folderIndex));
    mkdir(sprintf('.\\output%d\\fig', folderIndex));
end
saveas(gcf,[sprintf('.\\output%d\\fig\\', folderIndex),'histogram.fig'])
saveas(gcf,[sprintf('.\\output%d\\png\\', folderIndex),'histogram.png'])

disp('Process Finished.')

% The main script ends here.
% ---------------------------------------------------------------------------------------


% Analytical Solution for Binding Rates and Concentrations in the GEF(Sos)-GAP competition model.
% This function analytically solves a set of equations representing the dynamics of Sos-RasGDP, Sos-RasGTP, RasGDP, and RasGTP
% given the catalytic rates of GAP and GEF, feedback strengths, unbinding rates, target number of RasGTP, and total Ras.
% The function returns the concentrations of Sos-RasGDP, Sos-RasGTP complexes, free RasGDP, and the binding rate of SOS (GEF) to RasGDP (kOnGdp).
function [SosRasGdp,SosRasGtp,RasGdp,kOnGdp] = getKOnGtp(kCatGap, kCatGef, kOnGtp2KOnGdp, kOffGdp, kOffGtp, nTargetRasGtp, nRasTotal)

%setting all symbols involved in the equation
syms kOnGdp SosRasGdp SosRasGtp RasGdp real

% eqn1 = d(SosRasGdp)/dt == 0
eqn1 = kOnGdp * RasGdp - kOffGdp * SosRasGdp == 0;
% eqn2 = d(SosRasGdp)/dt = 0
eqn2 = kOnGtp2KOnGdp * kOnGdp * nTargetRasGtp - kOffGtp * SosRasGtp == 0;
% eqn3 = d(RasGtp)/dt == 0
eqn3 = kCatGef * (SosRasGdp + SosRasGtp) * RasGdp - kCatGap * nTargetRasGtp ...
    - kOnGtp2KOnGdp * kOnGdp * nTargetRasGtp + kOffGtp * SosRasGtp == 0;
% eqn4 = d(RasGdp)/dt == 0
eqn4 = -kCatGef * (SosRasGdp + SosRasGtp) * RasGdp + kCatGap * nTargetRasGtp ...
    - kOnGdp * RasGdp + kOffGdp * SosRasGdp == 0;
% eqn5 = total Ras conserve
eqn5 =  SosRasGdp + SosRasGtp + RasGdp + nTargetRasGtp == nRasTotal;

%Solve analytically
[SosRasGdp,SosRasGtp,RasGdp,kOnGdp] = solve([eqn1,eqn2,eqn3,eqn4,eqn5, ...
    SosRasGdp >= 0, 100 > SosRasGtp, SosRasGtp >= 0,RasGdp >= 0,kOnGdp >= 0], ...
    [SosRasGdp,SosRasGtp,RasGdp,kOnGdp]);

kOnGdp = double(kOnGdp);

end % from getKOnGtp()


% This function reorganizes recorded time trace data so that it aligns with each unit time, making it easier to analyze and visualize.
% It takes in a recorded time array, the input time trace data, the total data length, and the number of corrals as inputs.
function traceEveryUnitTime = reorganizeForEveryUnitTime(recordedTimeArray, inputTimeTrace, dataLength, nCorrals)
traceEveryUnitTime = zeros(dataLength,nCorrals);
for j = 0:1:dataLength % record for every (integer) unit time
    for i = 1:nCorrals
        jUnitTimesMinusRecordedTimes = j - recordedTimeArray(:,i);
        jUnitTimesMinusRecordedTimes(jUnitTimesMinusRecordedTimes<=0) = nan;
        [~, idxes] = min(jUnitTimesMinusRecordedTimes);
        traceEveryUnitTime(j+1,i) = inputTimeTrace(idxes,i);
    end
end
end % from reorganizeForEveryUnitTime()

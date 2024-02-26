close all;
% clear all;

% Set 1/kCatGEF as the unit time.
kCatGef = 1;
% Set the kcat of GAP
kCatGap = 2^-2; %
% Set the feedback strengths (propensity to bind to RasGtp than RasGdp)
kOnGtp2KOnGdp = 2^0;
% Set the SOS unbinding rate
kOffGdp = 2^-2;
% Set the koff of GEF from GTP and GDP the same (as default)
kOffGtp2KOffGdp = 1;
kOffGtp = kOffGdp * kOffGtp2KOffGdp;

rateLawXVal = 0.6;  % Based on this value, k_on's will be determined by solving the ODEs later.

% Set the grid length
gridLen = 1;

RAS_DENSITY = 1000;
N_CORRALS = 200;
N_MIN_DURATION = 30;
LOWEST_AVG_X_REQD_FOR_DIP_TEST = 0.1; % If avg x is too low, we cannot draw conclusions about unimodal/bimodality
N_TRACKED_VARIABLES = 5;
TYPE_OFFSET = 0; % Offset to determine type2 distribution when the higher gaussian mean is higher than the rate law x-value by the offset.

% Automatically name the output folder as "outputN", where N is the lowest
% number for such non-existing folder name.
folderIndex = 1;
while exist(sprintf('.\\output%d', folderIndex), 'dir')
    folderIndex = folderIndex + 1;
end
mkdir(sprintf('.\\output%d', folderIndex));

fprintf('\n kCatGap = %f \n', kCatGap);
fprintf('\n kOnGtp2KOnGdp = %f \n',kOnGtp2KOnGdp);
fprintf('\n kOffGdp = %f \n',kOffGdp);
fprintf('\n kOffGtp = %f \n',kOffGtp);

nRasTotal = RAS_DENSITY * gridLen^2;
nTargetRasGtp = rateLawXVal * nRasTotal;
[SosRasGdp,SosRasGtp,RasGdp,kOnGdp] = getKOnGtp(kCatGap, kCatGef, kOnGtp2KOnGdp, kOffGdp, kOffGtp, nTargetRasGtp, nRasTotal); %kon*[K], binding rate to RasGDP
str = sprintf('\n  Analytical solution for x=%f\n-> [kOnGdp = %f, kOnGtp = %f, RasGdp = %f, RasGtp = %f, SosRasGdp = %f, SosRasGtp = %f]\n', ...
    rateLawXVal, kOnGdp, kOnGtp2KOnGdp * kOnGdp, RasGdp, nRasTotal-RasGdp-SosRasGdp-SosRasGtp, SosRasGdp, SosRasGtp);
disp(str);
if length(kOnGdp) > 1
    disp('Error in finding deterministic solution.');
    return;
end
kOnGtp = kOnGtp2KOnGdp * kOnGdp; %binding rate to RasGTP

% Set the appropriate rxn time to run
duration = 1 + nTargetRasGtp / 2 / (kOffGtp + nTargetRasGtp * kOnGdp) * kOffGdp;
duration = duration * 2;
duration = round(duration);
if duration > 20000
    duration = 20000;
end

% Initialize random number generator
rng('shuffle')

nTrackedVars = 5; % tracked variables: (t,nRasGtp,nRasGdp,nSosRasGtp,nSosRasGdp)

% Set the time trace matrix for all tracked variabled and give me much
% time points (10 * duration (in unit times))
timeTraceAllCorrals = zeros(10 * double(duration), nTrackedVars * N_CORRALS);

% Iterate over each corral
for iCorral = 1:N_CORRALS
    % Print output message only under below condition
    if mod(iCorral, N_CORRALS/2) == 0
        str = sprintf(' Results: iCorral = %d out of %d', iCorral, N_CORRALS);
        disp(str)
    end

    % Initial condition setup
    t = 0;
    x_initial = 0;
    nRasGtp = nRasTotal * x_initial;
    nRasGdp = nRasTotal * (1-x_initial);
    nSosRasGtp = 0;
    nSosRasGdp = 0;

    % Record the initial state
    record = [t,nRasGtp,nRasGdp,nSosRasGtp,nSosRasGdp];

    % Record index
    recIdx = 0;

    % Insert the record to the dedicated place for the current corral
    timeTraceAllCorrals(recIdx+1,((iCorral-1)*nTrackedVars + 1) : iCorral*nTrackedVars) = record;

    % Gillespie algorithm
    while (t < duration) || (t < N_MIN_DURATION)
        r1 = kOnGtp * nRasGtp; % Rate of SOS binding to RasGtp
        r2 = kOnGdp * nRasGdp; % Rate of SOS binding to RasGdp
        r3 = kOffGtp * nSosRasGtp; % Rate of SOS unbinding from RasGtp
        r4 = kOffGdp * nSosRasGdp; % Rate of SOS unbinding from RasGdp
        r5 = kCatGef * (nSosRasGtp + nSosRasGdp) / gridLen^2 * nRasGdp / gridLen^2; % Rate of SOS converting RasGdp to RasGtp
        r6 = kCatGap / gridLen^2 * nRasGtp / gridLen^2; % Rate of GAP converting RasGtp to RasGdp

        % Combined rate
        rT=r1+r2+r3+r4+r5+r6;

        % Draw random time interval following exponential distribution
        drawT = rand;
        tau = 1/rT * log(1/drawT);

        % Advance current time
        t = t + tau;
        t = double(t);

        % Probability of each reaction for this time interval
        p1=r1/rT;
        p2=r2/rT;
        p3=r3/rT;
        p4=r4/rT;
        p5=r5/rT;

        % Randomly determine the occurred reaction
        drawReaction = rand;

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

        % Record the trajectory under the below condition
        if floor(t*10) > recIdx % recording every 0.1 unit time
            % Increment record index
            recIdx = recIdx + 1;

            % Record
            timeTraceAllCorrals(recIdx+1,(iCorral-1)*nTrackedVars + 1 : iCorral * nTrackedVars) ...
                = [t,nRasGtp,nRasGdp,nSosRasGtp,nSosRasGdp];
        end
    end
end % Full runtime for a corral end

% Extract time trace of each value from the master record
rasGtpTimeTrace = timeTraceAllCorrals(:,2:N_TRACKED_VARIABLES:N_CORRALS*N_TRACKED_VARIABLES);
sosGtpTimeTrace = timeTraceAllCorrals(:,4:N_TRACKED_VARIABLES:N_CORRALS*N_TRACKED_VARIABLES);
sosGdpTimeTrace = timeTraceAllCorrals(:,5:N_TRACKED_VARIABLES:N_CORRALS*N_TRACKED_VARIABLES);
recordedTimeArray = timeTraceAllCorrals(:,1:N_TRACKED_VARIABLES:N_CORRALS*N_TRACKED_VARIABLES);
dataLength = max(duration, N_MIN_DURATION);

% Reorganize for every unit time
nUnitTimes = 0:1:dataLength;
rasGtpEveryUnitTime = reorganizeForEveryUnitTime(recordedTimeArray, rasGtpTimeTrace, dataLength, N_CORRALS);
sosGtpEveryUnitTime = reorganizeForEveryUnitTime(recordedTimeArray, sosGtpTimeTrace, dataLength, N_CORRALS);
sosGdpEveryUnitTime = reorganizeForEveryUnitTime(recordedTimeArray, sosGdpTimeTrace, dataLength, N_CORRALS);

% Record the mean SOS-GTP and SOS-GDP numbers from all corrals
meanOfSosGtpEveryUnitTime = mean(sosGtpEveryUnitTime,2);
meanOfSosGdpEveryUnitTime = mean(sosGdpEveryUnitTime,2);

% Record the mean of normalized RasGTP (1 means all Ras are RasGTP, 0 means all Ras are RasGDP).
normalizedRasGtpEveryUnitTime = rasGtpEveryUnitTime/nRasTotal;
meanOfNormalizedRasGtpEveryUnitTime = mean(normalizedRasGtpEveryUnitTime,2);

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

% Check whether avg x-value is too low. If so, DIP test is not possible.
if sum(lastNormRasGtpArray) / length(lastNormRasGtpArray) > LOWEST_AVG_X_REQD_FOR_DIP_TEST
    bootStrapSampleSize = 100;
    % Run the DIP test storing the result to p-value. Low p-value indicates
    % non-unimodal distribution
    [~, p] = HartigansDipSignifTest(normalizedRasGtpEveryUnitTime(end,:), bootStrapSampleSize);
end

% Plot the result
figure;
plot(nSeconds,meanOfNormalizedRasGtpEverySecond,'LineWidth',2, 'Color', '#004C99');
set(gca, 'CLim',[0 0.2],'FontSize',12,'FontWeight','bold','LineWidth',2);
xlim ([0 size(nSeconds,2)])
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


function [SosRasGdp,SosRasGtp,RasGdp,kOnGdp] = getKOnGtp(kCatGap, kCatGef, kOnGtp2KOnGdp, kOffGdp, kOffGtp, nTargetRasGtp, nRasTotal)

%setting all symbols involved in the equation
syms kOnGdp SosRasGdp SosRasGtp RasGdp real

%eqn1=d(SosRasGdp)/dt
eqn1 = kOnGdp * RasGdp - kOffGdp * SosRasGdp == 0;
%eqn2=d(SosRasGdp)/dt
eqn2 = kOnGtp2KOnGdp * kOnGdp * nTargetRasGtp - kOffGtp * SosRasGtp == 0;
%eqn3=d(RasGtp)/dt
eqn3 = kCatGef * (SosRasGdp + SosRasGtp) * RasGdp - kCatGap * nTargetRasGtp ...
    - kOnGtp2KOnGdp * kOnGdp * nTargetRasGtp + kOffGtp * SosRasGtp == 0;
%eqn4=d(RasGdp)/dt
eqn4 = -kCatGef * (SosRasGdp + SosRasGtp) * RasGdp + kCatGap * nTargetRasGtp ...
    - kOnGdp * RasGdp + kOffGdp * SosRasGdp == 0;
%eqn5=total Ras conserve
eqn5 =  SosRasGdp + SosRasGtp + RasGdp + nTargetRasGtp == nRasTotal;

%Solve analytically
[SosRasGdp,SosRasGtp,RasGdp,kOnGdp] = solve([eqn1,eqn2,eqn3,eqn4,eqn5, ...
    SosRasGdp >= 0, 100 > SosRasGtp, SosRasGtp >= 0,RasGdp >= 0,kOnGdp >= 0], ...
    [SosRasGdp,SosRasGtp,RasGdp,kOnGdp]);

kOnGdp = double(kOnGdp);

end

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
end

% This script simulates Ras signaling dynamics incorporating the competition between GEF (Sos) and GAP on Ras, with additional consideration for inhibitor kinetics.
% Refer to the "main" script, which is devoid of inhibitor.
% The simulation utilizes a stochastic approach to model the molecular interactions and changes over time within a corraled membrane.
% Using Gillespie method, it dynamically tracks the system's evolution, specifically focusing on RasGTP and RasGDP concentrations and their interaction with Sos in the presence of inhibitors.
% Neil H. Kim and Albert A. Lee, 2024.

close all;  % Close all open figure windows.
clear all;  % Clear all variables from the workspace.
disp('-- Stochastic Simulation Initiating --');  % Notify start of simulation.

% Setup for output directory creation.
folderIndex = 1;
while exist(sprintf('.\\output%d', folderIndex), 'dir')
    folderIndex = folderIndex + 1;  % Increment folder index if it already exists.
end

% Simulation parameters initialization.
nSimulations = 200;  % Number of simulation runs.
gridLength = 1;  % Grid length for simulation, A=gridLength^2 represents area.

% Inhibitor kinetics parameters.
ki_on = 0.0001;  % On-rate for inhibitor binding.
ki_off = 0.0001;  % Off-rate for inhibitor unbinding.

% Binding/unbinding rate parameters
kOnGdp = (ki_off/(ki_off+ki_on))*[2.5e-06]; % binding rate of GEF to RasGDP as a function of ki_on and ki_off.
kOnGtp = 100 * kOnGdp; % binding rate of GEF to RasGTP, presumed significantly faster than to RasGDP
kOffGdp = 0.1; % Unbinding rate of SOS from SOS-RasGDP
kOffGtp = 0.1; % Unbinding rate of SOS from SOS-RasGTP

% Catalytic rate parameters.
nGapPerGrid = 1 * gridLength^2;  % GAP concentration.
kCatGef = 0.01 / gridLength^2;  % Catalytic rate of GEF.
kCatGap = 0.01 * 0.5 / gridLength^2;  % Catalytic rate of GAP, adjusted per grid area.

% Initial system conditions.
x_initial = 0;  % Initial proportion of RasGTP (among RasGDP and RasGTP combined).
nRasTotal = 1000 * gridLength^2;  % Total number of Ras molecules.

% Time-related parameters for the simulation.
duration = 5000;  % Total duration in unit times.
samplefactor = 100;  % Sampling frequency for recording the system state.
duration = round(duration);  % Ensure duration is an integer.
N_TRACKED_VARIABLES = 7;  % Number of variables being tracked in the simulation.

% Initialize array to record the time trace of the system.
recordedTimeTraceArray = zeros(10 * duration, N_TRACKED_VARIABLES * nSimulations);

rng('shuffle');  % Initialize random number generator for stochastic elements.

% Main simulation loop.
for iSimulation = 1:nSimulations
    % Periodic display of simulation progress.
    if mod(iSimulation, nSimulations / 10) == 0
        str = sprintf('nSimulations = %d out of %d', iSimulation, nSimulations);
        disp(str);
    end

    % Initialize state variables for each simulation run.
    t = 0;  % Initialize time in unit times.
    nRasGtp = nRasTotal * x_initial;  % Initial RasGTP count.
    nRasGdp = nRasTotal * (1 - x_initial);  % Initial RasGDP count.
    nSosRasGtp = 0;  % Initial count of Sos-RasGTP complexes.
    nSosRasGdp = 0;  % Initial count of Sos-RasGDP complexes.
    nSOSRti = 0;  % Initial count of inhibitor-bound Sos-RasGTP.
    nSOSRdi = 0;  % Initial count of inhibitor-bound Sos-RasGDP.
    nReactionNumber = 0;  % Track number of reactions occurred.

    % Record initial state.
    recordedTimeTraceArray(1, ((iSimulation - 1) * N_TRACKED_VARIABLES + 1):iSimulation * N_TRACKED_VARIABLES) = [t, nRasGtp, nRasGdp, nSosRasGtp, nSosRasGdp, nSOSRti, nSOSRdi];

    % Initialize record index for storing simulation results
    recIdx = 0;

    % Use the Gillespie algorithm for stochastic simulation, iterating until the simulation duration is reached
    while t < duration
        % Calculate rates for all reactions based on current state.
        r1 = kOnGtp * nRasGtp;
        r2 = kOnGdp * nRasGdp;
        r3 = kOffGtp * nSosRasGtp;
        r4 = kOffGdp * nSosRasGdp;
        r5 = kCatGef * (nSosRasGtp + nSosRasGdp) / gridLength^2 * nRasGdp;  % GEF catalyzed RasGDP to RasGTP conversion.
        r6 = kCatGap * nGapPerGrid / gridLength^2 * nRasGtp;  % GAP catalyzed RasGTP to RasGDP conversion.
        % Additional reactions for inhibitor binding/unbinding to/from Sos-Ras complexes and direct inhibitor effects on binding rates.
        r7 = ki_on * nSosRasGtp;
        r8 = ki_off * nSOSRti;
        r9 = ki_on * nSosRasGdp;
        r10 = ki_off * nSOSRdi;
        r11 = (ki_on / ki_off) * kOnGtp * nRasGtp;
        r12 = (ki_on / ki_off) * kOnGdp * nRasGdp;
        r13 = kOffGtp * nSOSRti;
        r14 = kOffGdp * nSOSRdi;
        rT = r1 + r2 + r3 + r4 + r5 + r6 + r7 + r8 + r9 + r10 + r11 + r12 + r13 + r14;  % Total rate sum.

        % Draw next reaction time and update system state based on drawn reaction.
        drawT = rand;
        tau = 1 / rT * log(1 / drawT);
        t = t + tau;
        nReactionNumber = nReactionNumber + 1;

        % Calculate probabilities for each reaction based on current rates and total rate.
        p1 = r1 / rT;
        p2 = r2 / rT;
        p3 = r3 / rT;
        p4 = r4 / rT;
        p5 = r5 / rT;
        p6 = r6 / rT;
        p7 = r7 / rT;
        p8 = r8 / rT;
        p9 = r9 / rT;
        p10 = r10 / rT;
        p11 = r11 / rT;
        p12 = r12 / rT;
        p13 = r13 / rT;

        drawReaction = rand;

        if drawReaction < p1
            %             rxn_case=1;
            nRasGtp = nRasGtp - 1;
            nSosRasGtp = nSosRasGtp + 1;
        elseif drawReaction < (p1+p2)
            %             rxn_case=2;
            nSosRasGdp = nSosRasGdp + 1;
            nRasGdp = nRasGdp - 1;
        elseif drawReaction < (p1+p2+p3)
            %             rxn_case=3;
            nSosRasGtp = nSosRasGtp - 1;
            nRasGtp = nRasGtp + 1;
        elseif drawReaction < (p1+p2+p3+p4)
            %             rxn_case=4;
            nSosRasGdp = nSosRasGdp - 1;
            nRasGdp = nRasGdp + 1;
        elseif drawReaction < (p1+p2+p3+p4+p5)
            %             rxn_case=5;
            nRasGtp=nRasGtp+1;
            nRasGdp=nRasGdp-1;
        elseif drawReaction < (p1+p2+p3+p4+p5+p6)
            %             rxn_case=6;
            nRasGtp=nRasGtp-1;
            nRasGdp=nRasGdp+1;
        elseif drawReaction < (p1+p2+p3+p4+p5+p6+p7)
            %             rxn_case=7;
            nSosRasGtp=nSosRasGtp-1;
            nSOSRti=nSOSRti+1;
        elseif drawReaction < (p1+p2+p3+p4+p5+p6+p7+p8)
            %             rxn_case=8;
            nSOSRti=nSOSRti-1;
            nSosRasGtp=nSosRasGtp+1;
        elseif drawReaction < (p1+p2+p3+p4+p5+p6+p7+p8+p9)
            %             rxn_case=9;
            nSosRasGdp=nSosRasGdp-1;
            nSOSRdi=nSOSRdi+1;
        elseif drawReaction < (p1+p2+p3+p4+p5+p6+p7+p8+p9+p10)
            %             rxn_case=10;
            nSOSRdi=nSOSRdi-1;
            nSosRasGdp=nSosRasGdp+1;
        elseif drawReaction < (p1+p2+p3+p4+p5+p6+p7+p8+p9+p10+p11)
            %             rxn_case=11;
            nRasGtp = nRasGtp - 1;
            nSOSRti = nSOSRti + 1;
        elseif drawReaction < (p1+p2+p3+p4+p5+p6+p7+p8+p9+p10+p11+p12)
            %             rxn_case=12;
            nSOSRdi = nSOSRdi + 1;
            nRasGdp = nRasGdp - 1;
        elseif drawReaction < (p1+p2+p3+p4+p5+p6+p7+p8+p9+p10+p11+p12+p13)
            %             rxn_case=13;
            nSOSRti = nSOSRti - 1;
            nRasGtp = nRasGtp + 1;
        else
            %             rxn_case=14;
            nSOSRdi = nSOSRdi - 1;
            nRasGdp = nRasGdp + 1;
        end

        % Record the updated state at specified time intervals
        if floor(t*10) > recIdx % Condition to record data at fixed intervals (every 0.1 unit time)
            % Update the timeTraceAllCorrals array with the new state
            recordedTimeTraceArray(recIdx+2,((iSimulation-1)*N_TRACKED_VARIABLES + 1):iSimulation*N_TRACKED_VARIABLES) ...
                = [t nRasGtp nRasGdp nSosRasGtp nSosRasGdp nSOSRti nSOSRdi];
            recIdx = recIdx + 1; % Increment the record index
        end
    end
end % of one reaction duration.

% Extracting specific time trace data from the comprehensive simulation results for further analysis
rasGtpTimeTrace = recordedTimeTraceArray(:,2:N_TRACKED_VARIABLES:nSimulations*N_TRACKED_VARIABLES);
sosGtpTimeTrace = recordedTimeTraceArray(:,4:N_TRACKED_VARIABLES:nSimulations*N_TRACKED_VARIABLES);
sosGdpTimeTrace = recordedTimeTraceArray(:,5:N_TRACKED_VARIABLES:nSimulations*N_TRACKED_VARIABLES);
recordedTimeArray = recordedTimeTraceArray(:,1:N_TRACKED_VARIABLES:nSimulations*N_TRACKED_VARIABLES);

% Reorganizing extracted time trace data for analysis at each unit time
nUnitTimes = 0:samplefactor:duration;
rasGtpEveryUnitTime = zeros(duration/samplefactor,nSimulations);
sosGtpEveryUnitTime = zeros(duration/samplefactor,nSimulations);
sosGdpEveryUnitTime = zeros(duration/samplefactor,nSimulations);

% Record for every (integer) UnitTime
for j = 0:samplefactor:duration 
    for i = 1:1:nSimulations
        jUnitTimesMinusRecordedTimes = j - recordedTimeArray(:,i);
        jUnitTimesMinusRecordedTimes(jUnitTimesMinusRecordedTimes<=0) = nan;
        [mins, idxes] = min(jUnitTimesMinusRecordedTimes);
        rasGtpEveryUnitTime(j/samplefactor+1,i) = rasGtpTimeTrace(idxes,i);
        sosGtpEveryUnitTime(j/samplefactor+1,i) = sosGtpTimeTrace(idxes,i);
        sosGdpEveryUnitTime(j/samplefactor+1,i) = sosGdpTimeTrace(idxes,i);
    end
end

% Calculating mean values of SOS-GTP and SOS-GDP across all corrals for each unit time
meanOfSosGtpEveryUnitTime = mean(sosGtpEveryUnitTime,2);
meanOfSosGdpEveryUnitTime = mean(sosGdpEveryUnitTime,2);

% Calculating and recording the mean of normalized RasGTP 
normalizedRasGtpEveryUnitTime = rasGtpEveryUnitTime/nRasTotal;
meanOfNormalizedRasGtpEveryUnitTime = mean(normalizedRasGtpEveryUnitTime,2);


% Saving processed data to .mat files for further analysis or visualization
% Save RasGTP every unit time as .mat file
filename = 'rasGtpEveryUnitTime.mat';
save(filename,'rasGtpEveryUnitTime')
[~, p] = HartigansDipSignifTest(normalizedRasGtpEveryUnitTime(end,:), 100); % The 2nd parameter is sample size of boot-strap

% Plot 1) the mean of normalized RasGTP, 2) all norm RasGTP traces, and 3) histogram of RasGTP/Ras values.
figure(1);

% 1) the mean of normalized RasGTP
subplot(3,1,1);
plot(nUnitTimes,meanOfNormalizedRasGtpEveryUnitTime,'LineWidth',2)
set(gca, 'CLim',[0 0.2],'FontSize',12,'FontWeight','bold','LineWidth',2);
xlim ([0 samplefactor*size(nUnitTimes,2)])
ylim ([0 1])
title('The mean of normalized RasGTP every UnitTime')
xlabel('Time (s)')
ylabel('RasGTP/Ras')

% 2) all norm RasGTP traces
subplot(3,1,2);
plot(nUnitTimes,normalizedRasGtpEveryUnitTime,'LineWidth',2);
lgnd = legend('Rt');
set(lgnd, 'FontSize', 8);
set(gca, 'CLim',[0 0.2],'FontSize',12,'FontWeight','bold','LineWidth',2);
xlim ([0 samplefactor*size(nUnitTimes,2)])
ylim ([0 1])
title('RasGTP every UnitTime')
xlabel('Time (s)')
ylabel('RasGTP/Ras')

% 3) histogram of RasGTP/Ras values
subplot(3,1,3);
histogram(normalizedRasGtpEveryUnitTime(end,:),'BinWidth',0.05,'Normalization','probability')
txt = ['p=',num2str(p,3)];
text(0.75,0.45,txt,'FontSize',16,'FontWeight','bold');
xlim ([0 1])
ylim ([0 1])
set(gca, 'CLim',[0 0.2],'FontSize',12,'FontWeight','bold','LineWidth',2);
title(sprintf('Histogram at t=%f', t))
xlabel('RasGTP/Ras')
ylabel('Frequency')

% Save plot
saveas(gcf,['Hist_end',num2str(kOnGdp),'.png'])

% Define symbolic variables for the concentrations of different species.
syms SRd SRt Rd Rt SRdi SRti % Sos-RasGDP, Sos-RasGTP, RasGDP, RasGTP, Sos-RasGDP-inhibitor, Sos-RasGTP-inhibitor

% Define equations based on the kinetics of the system, assuming steady state (derivative with respect to time equals zero).

% eqn1=dSRd/dt
eqn1 = kOnGdp*Rd - kOffGdp*SRd - ki_on*SRd + ki_off*SRdi== 0;
% eqn2=dSRt/dt
eqn2 = kOnGtp*Rt - kOffGtp*SRt - ki_on*SRt + ki_off*SRti== 0;
% eqn3=dRt/dt
eqn3 = kCatGef*(SRd+SRt)*Rd - kCatGap*Rt - kOnGtp*Rt - (ki_on/ki_off)*kOnGtp*Rt + kOffGtp*SRt + kOffGtp*SRti == 0;
% eqn4=dRd/dt
eqn4 = -kCatGef*(SRd+SRt)*Rd + kCatGap*Rt - kOnGdp*Rd - (ki_on/ki_off)*kOnGdp*Rd + kOffGdp*SRd + kOffGdp*SRdi == 0;
% eqn5=dSRti/dt
eqn5 = (ki_on/ki_off)*kOnGtp*Rt - kOffGtp*SRti + ki_on*SRt - ki_off*SRti== 0;
% eqn6=dSRdi/dt
eqn6 = (ki_on/ki_off)*kOnGdp*Rd - kOffGdp*SRdi + ki_on*SRd - ki_off*SRdi== 0;
% eqn7=total Ras conserve
eqn7 =  SRd + SRt + Rd + Rt + SRdi + SRti == nRasTotal;
% Solve analytically
[SRd,SRt,Rd,Rt,SRdi,SRti] = solve([eqn1,eqn2,eqn3,eqn4,eqn5,eqn6,eqn7,SRd>=0,SRt>=0,Rd>=0,Rt>0,SRdi>=0,SRti>=0],[SRd,SRt,Rd,Rt,SRdi,SRti]);

%record the solution in double precision
SRd_value(1:length(SRd)) = double(SRd);
SRt_value(1:length(SRt)) = double(SRt);
Rd_value(1:length(Rd)) = double(Rd);
Rt_value(1:length(Rt)) = double(Rt);
SRdi_value(1:length(SRd)) = double(SRdi);
SRti_value(1:length(SRt)) = double(SRti);

% Visualize the distribution of normalized RasGTP at the simulation's end.
figure(2)
histogram(normalizedRasGtpEveryUnitTime(end,:),'BinWidth',0.1,'Normalization','probability','EdgeColor','none','FaceColor','[1 0.5 0.5]')
xlim ([0 1])
ylim ([0 1])
set(gca, 'CLim',[0 0.2],'FontSize',30,'FontWeight','bold','LineWidth',4);
set(gcf,'Position',[100 100 500 500]);
yticks ([0 0.5 1]);
yticklabels ([0 0.5 1]);
hold on
% Highlight the theoretical equilibrium concentration of RasGTP.
plot([1 1]*(Rt_value/nRasTotal), [0,1],'--', 'LineWidth', 3, 'color', '#D90000');
xlabel('RasGTP/totalRas')
saveas(gcf,'histogram.tif')

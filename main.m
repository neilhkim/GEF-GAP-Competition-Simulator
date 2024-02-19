close all;
clear all;

% Main controls
kCatGef = 1; % Setting 1/kCatGEF as the unit time
kOffGtp2KOffGdp = 1; % Default setting

listLog2_kCatGap = -3:-1;%-2;
nKCatGaps = length(listLog2_kCatGap);

listLog2_kOnGtp2KOnGdp = 0;%0;%5:-2:1; % --> Feedback Strength (propensity to bind to RasGtp than RasGdp)
nFeedBackLevels = length(listLog2_kOnGtp2KOnGdp);

listLog2_kOffGdp = -3:-1;%-2;%-2:2:2; % --> Mean SOS unbinding rate
nKOffs = length(listLog2_kOffGdp);



TARGET_X_TO_MATCH_BY_ITERATION = 0.6;

IS_KON_ADJ_PERFORMED = false;%true;
ALPHA = 2; % Jump size of the target_x matching iteration
MAX_DELTA = 0.02;
N_MAX_KON_ADJ_ITERATION = 10;
N_CORRALS_FOR_KON_ADJ = 100;

RAS_DENSITY = 1000;
N_CORRALS = 200;
N_MIN_DURATION = 30;
LOWEST_AVG_X_REQD_FOR_DIP_TEST = 0.1; % If avg x is too low, we cannot draw conclusions about unimodal/bimodality
N_TRACKED_VARIABLES = 5;
N_TRACKED_VARS_IN_KON_ADJ = 2;
TYPE2_DISTRIB_FACTOR = 0; % Where to draw the line for TYPE2 criterion
% GRID_LENGTH = 2; % Now it's no longer a constant but an iteration variable

SHOW_PVAL_MATRIX = false;
SHOW_DISTRIB_TYPE_MATRIX = false;

SUBPLOT_MEAN_NORM_RASGTP = true;
SUBPLOT_MEAN_SOS = true;

% Automatically name the output folder as "outputN".
% Start from "output1", and if it already exists try, "output2", and so on.
folderIndex = 1;
while exist(sprintf('.\\output%d', folderIndex), 'dir')
    folderIndex = folderIndex + 1;
end
mkdir(sprintf('.\\output%d', folderIndex));

for Log2KCatGap = listLog2_kCatGap % [-2, -1, 0, 1, 2]
    kCatGap = 2^Log2KCatGap;
    fprintf('\n~~~~~~~~~~~~~~~ Log2KCatGap = %f ~~~~~~~~~~~~~~~~~~\n', Log2KCatGap);
    FB_kOff_pValues = NaN(nFeedBackLevels * nKOffs, 3); % Think of this as a collection of vectors (feedBack, kOff, pValue).
    FB_kOff_distributionType = NaN(nFeedBackLevels * nKOffs, 3); % Think of this as a collection of vectors (feedBack, kOff, distributionType).
    pValueIndex = 1;

    for Log2KOnGtp2KOnGdp = listLog2_kOnGtp2KOnGdp
        fprintf('\n .............. Log2KOnGtp2KOnGdp = %f ..............\n',Log2KOnGtp2KOnGdp);
        kOnGtp2KOnGdp = 2^Log2KOnGtp2KOnGdp;

        for Log2KOffGdp = listLog2_kOffGdp
            fprintf('\n `````````````` Log2KOffGdp = %f ``````````````\n',Log2KOffGdp);
            kOffGdp = 2^Log2KOffGdp;
            kOffGtp = kOffGdp * kOffGtp2KOffGdp;

            

            for gridlen = 1
                fprintf('target x = %.4f\n', TARGET_X_TO_MATCH_BY_ITERATION);
                x_initial = 0;
                nRasTotal = RAS_DENSITY * gridlen^2;
                nTargetRasGtp = TARGET_X_TO_MATCH_BY_ITERATION * nRasTotal;

%                 disp([kCatGap, kCatGef, kOnGtp2KOnGdp, kOffGdp, kOffGtp, nTargetRasGtp, nRasTotal]);


                [SosRasGdp,SosRasGtp,RasGdp,kOnGdp] = getKOnGtp(kCatGap, kCatGef, kOnGtp2KOnGdp, kOffGdp, kOffGtp, nTargetRasGtp, nRasTotal); %kon*[K], binding rate to RasGDP
                str = sprintf('\n  Analytical solution for x=%f\n-> [kOnGdp = %f, kOnGtp = %f, RasGdp = %f, RasGtp = %f, SosRasGdp = %f, SosRasGtp = %f]\n', ...
                    TARGET_X_TO_MATCH_BY_ITERATION, kOnGdp, kOnGtp2KOnGdp * kOnGdp, RasGdp, nRasTotal-RasGdp-SosRasGdp-SosRasGtp, SosRasGdp, SosRasGtp);
                disp(str);
                if length(kOnGdp) > 1
                    disp('Error in finding deterministic solution.');
                    return;
                end
                kOnGtp = kOnGtp2KOnGdp * kOnGdp; %binding rate to RasGTP

                % Estimate rxn time
                duration = 1 + nTargetRasGtp / 2 / (kOffGtp + nTargetRasGtp * kOnGdp) * kOffGdp;
                duration = duration * 2; % Arbitrary
                duration = round(duration);
                if duration > 20000
                    duration = 20000;
                end
                rng('shuffle')

                predictXPostIteration = TARGET_X_TO_MATCH_BY_ITERATION;

                if IS_KON_ADJ_PERFORMED == true

                    if gridlen == 1 % Just to do kOn-adjustment only once.
                        % START ----------- Test run for kOn adjustment ------------- %
                        if IS_KON_ADJ_PERFORMED == true
                            if ~exist(sprintf('.\\output%d\\adjusting', folderIndex), 'dir')
                                mkdir(sprintf('.\\output%d\\adjusting', folderIndex));
                            end
                            if ~exist(sprintf('.\\output%d\\adjusting\\fig', folderIndex), 'dir')
                                mkdir(sprintf('.\\output%d\\adjusting\\fig', folderIndex));
                            end
                            if ~exist(sprintf('.\\output%d\\adjusting\\png', folderIndex), 'dir')
                                mkdir(sprintf('.\\output%d\\adjusting\\png', folderIndex));
                            end
                            fprintf('-- Starting kOn Adjustment. kOnGtp is now %.8f\n', kOnGtp);
                            disp('- Testing to see what xIS_KON_ADJ_PERFORMED we get with the current kOnGdp & kOnGtp');

                            kOnGtpArray = zeros(N_CORRALS_FOR_KON_ADJ,1);
                            deltaArray = zeros(N_CORRALS_FOR_KON_ADJ,1);

                            for kOnAdjIdx = 1:N_MAX_KON_ADJ_ITERATION

                                kOnGtpArray(kOnAdjIdx) = kOnGtp;
                                
                                % I tried making the below (about 100 lines)
                                % routine as a local function call, but for
                                % some reason it slowed the run to an
                                % enormous degree... [2022-04-06 NKIM]
                                % Now for the actual result

                nMinDuration = N_MIN_DURATION;

                if IS_KON_ADJ_PERFORMED == true
                    nTrackedVars = 2; % (t, nRasGtp)
                    nCorrals = N_CORRALS_FOR_KON_ADJ;
                else
                
                    nTrackedVars = 5; % (t,nRasGtp,nRasGdp,nSosRasGtp,nSosRasGdp)
                    nCorrals = N_CORRALS;
                end

                timeTraceAllCorrals = zeros(10 * double(duration), nTrackedVars * nCorrals);

                for iCorral = 1:nCorrals
                    if mod(iCorral, nCorrals/2) == 0

                        if IS_KON_ADJ_PERFORMED == true
                            str = sprintf('=== Iteration %d, kOn Testing = %d out of %d === ', kOnAdjIdx, iCorral, nCorrals);
                        else
                            str = sprintf(' Results: iCorral = %d out of %d  |||  FB %d kOff %d kCatGap %.5f', iCorral, N_CORRALS, Log2KOnGtp2KOnGdp, Log2KOffGdp, Log2KCatGap);
                        end
                        disp(str)
                    end

                    t = 0;
                    nRasGtp = nRasTotal * x_initial; % rasGTP
                    nRasGdp = nRasTotal * (1-x_initial);
                    nSosRasGtp = 0; % SosRasGtp
                    nSosRasGdp = 0; % SosRasGdp

                    % Writing the initial state
                    if IS_KON_ADJ_PERFORMED == true
                        input = [t,nRasGtp];
                    else
                        input = [t,nRasGtp,nRasGdp,nSosRasGtp,nSosRasGdp];
                    end
                    timeTraceAllCorrals(1,((iCorral-1)*nTrackedVars + 1) : iCorral*nTrackedVars) = input;

                    recIdx = 0;
                    while (t < duration) || (t < nMinDuration)
                        r1 = kOnGtp * nRasGtp;
                        r2 = kOnGdp * nRasGdp;
                        r3 = kOffGtp * nSosRasGtp;
                        r4 = kOffGdp * nSosRasGdp;
                        r5 = kCatGef * (nSosRasGtp + nSosRasGdp) / gridlen^2 * nRasGdp / gridlen^2;
                        r6 = kCatGap / gridlen^2 * nRasGtp / gridlen^2;
                        rT=r1+r2+r3+r4+r5+r6;
                        drawT = rand;
                        tau = 1/rT * log(1/drawT);
                        t = t + tau;
                        p1=r1/rT;
                        p2=r2/rT;
                        p3=r3/rT;
                        p4=r4/rT;
                        p5=r5/rT;
                        
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

                        if floor(t*10) > recIdx % recording every 0.1 second

                            if IS_KON_ADJ_PERFORMED == true
                                 timeTraceAllCorrals(recIdx+1,(iCorral-1)*nTrackedVars + 1: iCorral * nTrackedVars) = [t,nRasGtp];
                            else
                                 timeTraceAllCorrals(recIdx+1,(iCorral-1)*nTrackedVars + 1 : iCorral * nTrackedVars) = [t,nRasGtp,nRasGdp,nSosRasGtp,nSosRasGdp];
                            end

%                             if nSosRasGtp + nSosRasGdp > 0
%                                 disp(double([t,nRasGtp,nRasGdp,nSosRasGtp,nSosRasGdp]));
%                             end


                            recIdx = recIdx + 1;
%                             if recIdx > 20
%                                 disp('stop')
%                             end
%                             timeTraceAllCorrals(recIdx+1,((iCorral-1)*nTrackedVars + 1) : iCorral*nTrackedVars) = input;
                        end
                    end
                end


                            
                                rasGtpTimeTrace = timeTraceAllCorrals(:, 2 : N_TRACKED_VARS_IN_KON_ADJ : N_CORRALS_FOR_KON_ADJ*N_TRACKED_VARS_IN_KON_ADJ);
                                recordedTimeArray = timeTraceAllCorrals(:, 1 : N_TRACKED_VARS_IN_KON_ADJ : N_CORRALS_FOR_KON_ADJ*N_TRACKED_VARS_IN_KON_ADJ);
                                dataLength = max(duration, N_MIN_DURATION);
                                nSeconds = 0:1:dataLength;
                                rasGtpEverySecond = reorganizeForEverySecond(recordedTimeArray, rasGtpTimeTrace, dataLength, N_CORRALS_FOR_KON_ADJ);

                                normalizedRasGtpEverySecond = rasGtpEverySecond / nRasTotal;
                                meanOfNormalizedRasGtpEverySecond = mean(normalizedRasGtpEverySecond,2);
                                meanTestX = mean(mean(normalizedRasGtpEverySecond(round(end/2):end,:)));
                                %                     f1 = figure(1);
                                figure(999);
                                hold on;
                                ccmap = turbo(N_MAX_KON_ADJ_ITERATION);
                                ccmap = ccmap(end:-1:1,:);
                                plot(nSeconds,meanOfNormalizedRasGtpEverySecond,'LineWidth',2, 'Color', ccmap(kOnAdjIdx,:));

                                if kOnAdjIdx == 1 % Starting plotting the kOn adjustment figure
                                    strSgtitle = sprintf('kOnGTP = %.4f kOnGDP\n kOffGdp = %.4f \n kOffGtp = %.4f \n kCatGAP = %.4f \n Unit Time: 1/kCatGEF', kOnGtp2KOnGdp, kOffGdp, kOffGtp, kCatGap);
                                    sgtitle(strSgtitle, 'FontSize',12);
                                    plot([0 dataLength], [TARGET_X_TO_MATCH_BY_ITERATION, TARGET_X_TO_MATCH_BY_ITERATION], '--', 'LineWidth',2,'color', '#66B2FF');
                                    set(gca, 'CLim',[0 0.2],'FontSize',12,'FontWeight','bold','LineWidth',2);
                                    xlim ([0 size(nSeconds,2)])
                                    ylim ([0 1])
                                    xlabel('Time (s)')
                                    ylabel('RasGTP/Ras')
                                    txt = sprintf('target x = %.4f', TARGET_X_TO_MATCH_BY_ITERATION);
                                    text(dataLength*0.05,TARGET_X_TO_MATCH_BY_ITERATION*1.1,txt,'FontSize',12,'FontWeight','bold', 'horizontalalignment', 'left', 'verticalalignment', 'bottom', 'color','#66B2FF');
                                    %                         txt2 = sprintf('kOnGtp = %.8f', kOnGtp);
                                    %                         text(dataLength*0.95,TARGET_X_TO_MATCH_BY_ITERATION*1.1,txt2,'FontSize',12,'FontWeight','bold', 'horizontalalignment', 'right', 'verticalalignment', 'bottom');
                                    %                         colorbar;
                                end
                                hold off;
                                delta = TARGET_X_TO_MATCH_BY_ITERATION - meanTestX;
                                fprintf('---- delta = TARGET_X - meanTestX = %.6f\n', delta);
                                deltaArray(kOnAdjIdx) = delta;

                                if abs(delta) < MAX_DELTA %/ 2^onRateRatioExponent
                                    fprintf('---- abs(delta) is smaller than MAX_DELTA = %f\n',MAX_DELTA);
                                    if kOnAdjIdx > 1
                                        % kOnGtpArray
                                        figure
                                        plot(0:kOnAdjIdx-1,kOnGtpArray(1:kOnAdjIdx),'LineWidth',2, 'Color', 'b', 'marker', '*');
                                        title('kOnGtp Adjustment Trace')
                                        xlim([0 max(kOnAdjIdx-1,1)]);
                                        xlabel('Iteration')
                                        ylabel('kOnGtp')

                                        saveas(gcf,[sprintf('.\\output%d\\adjusting\\png\\kOnGtpArray-Log2KCatGap=%.4f_Log2KOnGtp2Gdp=%d_Log2KOffGdp=%.2f', folderIndex, Log2KCatGap, Log2KOnGtp2KOnGdp, Log2KOffGdp),'.png'])
                                        saveas(gcf,[sprintf('.\\output%d\\adjusting\\fig\\kOnGtpArray-Log2KCatGap=%.4f_Log2KOnGtp2Gdp=%d_Log2KOffGdp=%.2f', folderIndex, Log2KCatGap, Log2KOnGtp2KOnGdp, Log2KOffGdp),'.fig'])
                                        % deltaArray
                                        figure
                                        plot(1:kOnAdjIdx,deltaArray(1:kOnAdjIdx),'LineWidth',2, 'Color', 'g', 'marker', '*');
                                        title('delta Trace')
                                        xlim([1 max(kOnAdjIdx,1)]);
                                        xlabel('Iteration')
                                        ylabel('delta')
                                        saveas(gcf,[sprintf('.\\output%d\\adjusting\\deltaArray-Log2KCatGap=%.4f_Log2KOnGtp2Gdp=%d_Log2KOffGdp=%.2f', folderIndex, Log2KCatGap, Log2KOnGtp2KOnGdp, Log2KOffGdp),'.png'])
                                        saveas(gcf,[sprintf('.\\output%d\\adjusting\\deltaArray-Log2KCatGap=%.4f_Log2KOnGtp2Gdp=%d_Log2KOffGdp=%.2f', folderIndex, Log2KCatGap, Log2KOnGtp2KOnGdp, Log2KOffGdp),'.fig'])
                                    end
                                    figure(999);
                                    saveas(gcf,[sprintf('.\\output%d\\adjusting\\Traces-Log2KCatGap=%.4f_Log2KOnGtp2Gdp=%d_Log2KOffGdp=%.2f_kOffRatio=%.2f', folderIndex, Log2KCatGap, Log2KOnGtp2KOnGdp, Log2KOffGdp, kOffGtp2KOffGdp),'.png'])
                                    break;
                                end
                                %                     close(f1);
                                kOnGtp = kOnGtp * (1 + ALPHA * delta);
                                kOnGdp = kOnGtp / kOnGtp2KOnGdp;
                                fprintf('[[ --- kOn Adjusted. kOnGtp is now %.8f-- ]]\n', kOnGtp);
                            end
                        end
                        fprintf('******* kOnGdp Adjustment Complete. ****************************************************************\n******* Beginning simulation with kOnGdp = %.8f, kOnGtp = %.8f *************************\n', kOnGdp, kOnGtp);

                        nTargetRasGtp = getNTargetRasGtp(kCatGap, kCatGef, kOnGdp, kOnGtp2KOnGdp, kOffGdp, kOffGtp, nRasTotal);
                        predictXPostIteration = double(nTargetRasGtp/nRasTotal); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                        disp(['The rate equation predicted x value as result of kOnGdp iteration is ', num2str(predictXPostIteration)]);
                        fprintf('Analytical solution for kOnGdp=%f\n-> [newly predicted X = %f, kOnGtp = %f, \n    RasGdp = %f, RasGtp = %f, \n    SosRasGdp = %f, SosRasGtp = %f,\n    kCatGap=%f, Unit Time: 1/kCatGEF', ...
                            kOnGdp, predictXPostIteration, kOnGtp,                RasGdp, nRasTotal-RasGdp-SosRasGdp-SosRasGtp, SosRasGdp, SosRasGtp, kCatGap);
                        % END ----------- Test run for kOn adjustment ------------- %
                        close all;
                    end
                end
                

                % Now for the actual result

                nMinDuration = N_MIN_DURATION;

                IS_KON_ADJ_PERFORMED = IS_KON_ADJ_PERFORMED;
                if IS_KON_ADJ_PERFORMED == true
                    nTrackedVars = 2; % (t, nRasGtp)
                    nCorrals = N_CORRALS_FOR_KON_ADJ;
                else
                
                    nTrackedVars = 5; % (t,nRasGtp,nRasGdp,nSosRasGtp,nSosRasGdp)
                    nCorrals = N_CORRALS;
                end

                timeTraceAllCorrals = zeros(10 * double(duration), nTrackedVars * nCorrals);

                for iCorral = 1:nCorrals
                    if mod(iCorral, nCorrals/2) == 0

                        if IS_KON_ADJ_PERFORMED == true
                            str = sprintf('=== Iteration %d, kOn Testing = %d out of %d === ', kOnAdjIdx, iCorral, nCorrals);
                        else
                            str = sprintf(' Results: iCorral = %d out of %d  |||  FB %d kOff %d kCatGap %.5f', iCorral, N_CORRALS, Log2KOnGtp2KOnGdp, Log2KOffGdp, Log2KCatGap);
                        end
                        disp(str)
                    end

                    t = 0;
                    nRasGtp = nRasTotal * x_initial; % rasGTP
                    nRasGdp = nRasTotal * (1-x_initial);
                    nSosRasGtp = 0; % SosRasGtp
                    nSosRasGdp = 0; % SosRasGdp

                    % Writing the initial state
                    if IS_KON_ADJ_PERFORMED == true
                        input = [t,nRasGtp];
                    else
                        input = [t,nRasGtp,nRasGdp,nSosRasGtp,nSosRasGdp];
                    end

                    recIdx = 0;
                    timeTraceAllCorrals(recIdx+1,((iCorral-1)*nTrackedVars + 1) : iCorral*nTrackedVars) = input;
                    while (t < duration) || (t < nMinDuration)
                        r1 = kOnGtp * nRasGtp; % Rate of SOS binding to RasGtp
                        r2 = kOnGdp * nRasGdp; % Rate of SOS binding to RasGdp
                        r3 = kOffGtp * nSosRasGtp; % Rate of SOS unbinding from RasGtp
                        r4 = kOffGdp * nSosRasGdp; % Rate of SOS unbinding from RasGdp
                        r5 = kCatGef * (nSosRasGtp + nSosRasGdp) / gridlen^2 * nRasGdp / gridlen^2; % Rate of SOS converting RasGdp to RasGtp
                        r6 = kCatGap / gridlen^2 * nRasGtp / gridlen^2; % Rate of GAP converting RasGtp to RasGdp

                        rT=r1+r2+r3+r4+r5+r6;
                        drawT = rand;
                        tau = 1/rT * log(1/drawT);
                        t = t + tau;
                        t = double(t);
                        p1=r1/rT;
                        p2=r2/rT;
                        p3=r3/rT;
                        p4=r4/rT;
                        p5=r5/rT;

%                         if nSosRasGtp > 0
%                             disp('waitup')
%                         end

                        
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

                        if floor(t*10) > recIdx % recording every 0.1 second

                             recIdx = recIdx + 1;
                            
                            if IS_KON_ADJ_PERFORMED == true
                                 timeTraceAllCorrals(recIdx+1,(iCorral-1)*nTrackedVars + 1: iCorral * nTrackedVars) ...
                                     = [t,nRasGtp];
                            else
                                 timeTraceAllCorrals(recIdx+1,(iCorral-1)*nTrackedVars + 1 : iCorral * nTrackedVars) ...
                                     = [t,nRasGtp,nRasGdp,nSosRasGtp,nSosRasGdp];
                            end

%                             if nSosRasGtp + nSosRasGdp > 0
%                                 disp(double([t,nRasGtp,nRasGdp,nSosRasGtp,nSosRasGdp]));
%                             end


                           

%                             if recIdx > 20
%                                 disp('stop')
%                             end
%                             timeTraceAllCorrals(recIdx+1,((iCorral-1)*nTrackedVars + 1) : iCorral*nTrackedVars) = input;
                        end
                    end
                end
                rasGtpTimeTrace = timeTraceAllCorrals(:,2:N_TRACKED_VARIABLES:N_CORRALS*N_TRACKED_VARIABLES);
                sosGtpTimeTrace = timeTraceAllCorrals(:,4:N_TRACKED_VARIABLES:N_CORRALS*N_TRACKED_VARIABLES);
                sosGdpTimeTrace = timeTraceAllCorrals(:,5:N_TRACKED_VARIABLES:N_CORRALS*N_TRACKED_VARIABLES);
                recordedTimeArray = timeTraceAllCorrals(:,1:N_TRACKED_VARIABLES:N_CORRALS*N_TRACKED_VARIABLES);
                dataLength = max(duration, N_MIN_DURATION);

%                 losers = rasGtpTimeTrace(:,rasGtpTimeTrace(end,:)==0);

                nSeconds = 0:1:dataLength;
                rasGtpEverySecond = reorganizeForEverySecond(recordedTimeArray, rasGtpTimeTrace, dataLength, N_CORRALS);
                sosGtpEverySecond = reorganizeForEverySecond(recordedTimeArray, sosGtpTimeTrace, dataLength, N_CORRALS);
                sosGdpEverySecond = reorganizeForEverySecond(recordedTimeArray, sosGdpTimeTrace, dataLength, N_CORRALS);
                
                meanOfSosGtpEverySecond = mean(sosGtpEverySecond,2);
                meanOfSosGdpEverySecond = mean(sosGdpEverySecond,2);

                normalizedRasGtpEverySecond = rasGtpEverySecond/nRasTotal;
                meanOfNormalizedRasGtpEverySecond = mean(normalizedRasGtpEverySecond,2);

                if ~exist(sprintf('.\\output%d\\variables', folderIndex), 'dir')
                    mkdir(sprintf('.\\output%d\\variables', folderIndex));
                end
                filename = sprintf('.\\output%d\\variables\\normalizedRasGtpEverySecond.mat', folderIndex);
                save(filename,'normalizedRasGtpEverySecond');
                

                % Check if it was close to all rasGTP=0
                lastNormRasGtpArray = normalizedRasGtpEverySecond(end,:);
                filename = sprintf('.\\output%d\\variables\\endPointNormXDistribution_log2_of_kOn%d_kOff%d_kCatGAP%.3f.mat', folderIndex, Log2KOnGtp2KOnGdp, Log2KOffGdp, Log2KCatGap);
                save(filename,'lastNormRasGtpArray');
                if sum(lastNormRasGtpArray) / length(lastNormRasGtpArray) > LOWEST_AVG_X_REQD_FOR_DIP_TEST
                    [~, p] = HartigansDipSignifTest(normalizedRasGtpEverySecond(end,:), 100); % second parameter is sample size of boot-strap
                    FB_kOff_pValues(pValueIndex,:) = [Log2KOnGtp2KOnGdp Log2KOffGdp p];
                else
                    FB_kOff_pValues(pValueIndex,:) = [Log2KOnGtp2KOnGdp Log2KOffGdp NaN];
                end

                hh = figure(pValueIndex);
                hh.Position = [0 0 600 1000];
                strSgtitle = sprintf('kOnRatiosGTP = %.4f kOnRatiosGDP\n kOffRasGtp = %.4f kCatGEF\n kOffGtp = %.4f kOffGdp\n kCatGAP = %.4f kCatGEF', kOnGtp2KOnGdp, kOffGdp, kOffGtp2KOffGdp, kCatGap);
                sgtitle(strSgtitle, 'FontSize',12);

                nSubPlots = 1;
                if SUBPLOT_MEAN_NORM_RASGTP == true
                    nSubPlots = nSubPlots + 1;
                end
                if SUBPLOT_MEAN_SOS == true
                   nSubPlots = nSubPlots + 1;
                end 

                subplotIdx = 1;

                randomCorralIdx = randi(nCorrals);

                if SUBPLOT_MEAN_NORM_RASGTP == true
                    subplot(nSubPlots,1,subplotIdx);
                    subplotIdx = subplotIdx + 1;
                    
                    list = rasGtpEverySecond(:,randomCorralIdx);
                    %  plot(nSeconds,meanOfNormalizedRasGtpEverySecond,'LineWidth',2, 'Color', '#004C99');
                    plot(nSeconds,list,'LineWidth',2, 'Color', '#004C99');
                    hold on;
                    plot([0 dataLength], [TARGET_X_TO_MATCH_BY_ITERATION, TARGET_X_TO_MATCH_BY_ITERATION], '--', 'LineWidth',2,'color', '#66B2FF');
                    hold off;
                    set(gca, 'CLim',[0 0.2],'FontSize',12,'FontWeight','bold','LineWidth',2);
                    xlim ([0 size(nSeconds,2)])
%                     ylim ([0 1])
                    ylim ([0 nRasTotal])
%                     title('The mean of normalized RasGTP every second')
title('Example RasGTP time trace');
                    xlabel('Time (s)')
                    ylabel('RasGTP/Ras')
                    txt = sprintf('target x = %.4f', TARGET_X_TO_MATCH_BY_ITERATION);
%                     text(4,TARGET_X_TO_MATCH_BY_ITERATION*1.1,txt,'FontSize',12,'FontWeight','bold', 'verticalalignment', 'bottom', 'color','#66B2FF');
                end

                if SUBPLOT_MEAN_SOS == true
                    subplot(nSubPlots,1,subplotIdx);
                    subplotIdx = subplotIdx + 1;
%                     plot(nSeconds,meanOfSosGtpEverySecond,'LineWidth',2, 'color', '#FF7373');
                    list = sosGtpEverySecond(:,randomCorralIdx);
                    plot(nSeconds,list,'LineWidth',2, 'color', '#FF7373');
                    hold on;
                    list = sosGdpEverySecond(:,randomCorralIdx);
%                     plot(nSeconds,meanOfSosGdpEverySecond,'LineWidth',2, 'color', '#bada55');
                    plot(nSeconds,list,'LineWidth',2, 'color', '#bada55');
                    hold off;
                    lgnd = legend('SRt', 'SRd');
                    set(lgnd, 'FontSize', 8);
                    set(gca, 'CLim',[0 0.2],'FontSize',12,'FontWeight','bold','LineWidth',2);
                    xlim ([0 size(nSeconds,2)])
%                     ylim ([0 max(ceil(max([meanOfSosGtpEverySecond; meanOfSosGdpEverySecond])),1)])
                    ylim ([0 max([max(max([sosGtpEverySecond(:,randomCorralIdx); sosGdpEverySecond(:,randomCorralIdx)])),1])])
%                     title('The mean of recruited SOS number every second')
                    title('Example SOS number time trace')
                    xlabel('Time (s)')
                    ylabel('# of SOS')
                end

                subplot(nSubPlots,1,subplotIdx);
%                 figure;
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
                secondMean = nan;
                % --------------------
                % var distribType
                % --------------------
                % 0: Unimodal
                % 1: Zero & Mean
                % 2: Zero & 2Mean
                if sum(lastNormRasGtpArray) / length(lastNormRasGtpArray) > LOWEST_AVG_X_REQD_FOR_DIP_TEST
                    plot(predictXPostIteration*[1 1], [0 1],'--', 'LineWidth', 2, 'color', '#D90000');
                    txt = ['p = ',num2str(p,3),','];
                    text(0.1,max(histYArray),txt,'FontSize',12,'FontWeight','bold','HorizontalAlignment', 'left','VerticalAlignment', 'top');
                    if p < 0.05
                        [distribType, secondMean] = histogramGaussFit(sprintf('.\\output%d\\normalizedRasGtpEverySecond.mat', folderIndex), histXArray,histYArray,predictXPostIteration,TYPE2_DISTRIB_FACTOR);
                        plot([1 1]*(predictXPostIteration + TYPE2_DISTRIB_FACTOR/predictXPostIteration), [0,1], '-.', 'LineWidth', 2, 'color', '#e03fd8');
                        plot([1 1]*secondMean, [0,1], '-', 'LineWidth', 2, 'color', '#67349C');
                        text(secondMean*1.03, max(histYArray)*0.5, sprintf('%.3f',secondMean),'FontSize',12,'FontWeight','bold','HorizontalAlignment', 'left','VerticalAlignment', 'bottom', 'color', '#67349C');
                        text((predictXPostIteration + TYPE2_DISTRIB_FACTOR/predictXPostIteration)*1.03, max(histYArray)*0.2, 'new x','FontSize',12,'FontWeight','bold','HorizontalAlignment', 'left','VerticalAlignment', 'bottom', 'color', '#e03fd8');
                    else
                        distribType = 0;
                        text(0.5, max(histYArray)*0.2, 'Unimodal','FontSize',12,'FontWeight','bold','HorizontalAlignment', 'center','VerticalAlignment', 'bottom');
                    end
                    FB_kOff_distributionType(pValueIndex,:) = [Log2KOnGtp2KOnGdp Log2KOffGdp distribType];
                    txt2 = ['Stochastic Bistability Type = ',num2str(FB_kOff_distributionType(pValueIndex, 3))];
                    text(0.5,max(histYArray),txt2,'FontSize',12,'FontWeight','bold','HorizontalAlignment', 'center','VerticalAlignment', 'top');
                else
                    FB_kOff_distributionType(pValueIndex,:) = [Log2KOnGtp2KOnGdp Log2KOffGdp NaN];
                    txt2 = 'x-value too low to draw conclusion';
                    text(0.5,max(histYArray),txt2,'FontSize',12,'FontWeight','bold','HorizontalAlignment', 'center','VerticalAlignment', 'top');
                end
                hold off;

                filename = sprintf('.\\output%d\\variables\\secondMean_kOn%d_kOff%d_kCatGAP%d.mat', folderIndex, Log2KOnGtp2KOnGdp, Log2KOffGdp, Log2KCatGap);
                save(filename,'secondMean', 'predictXPostIteration');

                if ~exist(sprintf('.\\output%d\\png', folderIndex), 'dir')
                    mkdir(sprintf('.\\output%d\\png', folderIndex));
                    %                 mkdir(sprintf('.\\output%d\\pdf', folderIndex));
                    mkdir(sprintf('.\\output%d\\fig', folderIndex));
                end
                saveas(gcf,[sprintf('.\\output%d\\fig\\kCatGap2Gef=%.4f_kOnRatioExpo=%d_kOffGtp2kCat=%.2f', folderIndex, Log2KCatGap, Log2KOnGtp2KOnGdp, Log2KOffGdp),'.fig'])
                saveas(gcf,[sprintf('.\\output%d\\png\\kCatGap2Gef=%.4f_kOnRatioExpo=%d_kOffGtp2kCat=%.2f', folderIndex, Log2KCatGap, Log2KOnGtp2KOnGdp, Log2KOffGdp),'.png'])
                %             saveas(gcf,[sprintf('.\\output%d\\pdf\\kOnRatioExpo=%d_kOffGtp2kCat=%.2f_kCatGap2Gef=%.4f', folderIndex, onRateRatioExponent, offRateLevel, kCatGap2Gef),'.pdf'])

                pValueIndex = pValueIndex + 1;
            end
        end
    end
    %     for i=1:pValueIndex-1
    %         close(i);
    %     end
    %     close all;
    pValueFileName = sprintf('.\\output%d\\pValues_kOnRatioT2D_kOffD2KCat_kCatGap2Gef%.4f.csv', folderIndex,Log2KCatGap);
    writematrix(FB_kOff_pValues, pValueFileName);
    distribTypeFileName = sprintf('.\\output%d\\distribType_kOnRatioT2D_kOffD2KCat_kCatGap2Gef%.4f.csv', folderIndex,Log2KCatGap);
    writematrix(FB_kOff_distributionType, distribTypeFileName);

    if SHOW_PVAL_MATRIX == true
        fid = fopen(pValueFileName,'rt');
        C = textscan(fid, '%f %f %f', 'Delimiter', ',');
        OnRateRatioExponentAxis = reshape(C{1},nKOffs,[]); OnRateRatioExponentAxis = OnRateRatioExponentAxis(1,:);
        offRateLevelAxis = reshape(C{2},nKOffs,[]); offRateLevelAxis = offRateLevelAxis(:,1);
        pValuesHeight = reshape(C{3},nKOffs,[]);
        pValuesHeight = pValuesHeight';
        figure;
        pValuesHeight(pValuesHeight==0) = 0.001;
        customColorMap = [ones(1,3);summer(1000+1)];
        h = bar3(pValuesHeight);
        colormap(customColorMap);
        caxis([0 1]);
        strTitle = sprintf('Folder #%d\nTargetX = %.3f\nkOffGtp = %.4f kOffGdp\n kCatGAP = %.4f kCatGEF',  folderIndex,TARGET_X_TO_MATCH_BY_ITERATION, kOffGtp2KOffGdp, Log2KCatGap);
        title(strTitle);
        for k = 1:length(h)
            zdata = h(k).ZData;
            h(k).CData = zdata;
            h(k).FaceColor = 'interp';
        end
        colorbar;
        xlabel('offRateLevel'), ylabel('OnRateRatioExponent'), zlabel('pValues');
        if sum(isnan(OnRateRatioExponentAxis))
            OnRateRatioExponentAxis = listLog2_kOnGtp2KOnGdp;
        end
        if sum(isnan(offRateLevelAxis))
            offRateLevelAxis = listLog2_kOffGdp;
        end
        set(gca, 'XTickLabel',offRateLevelAxis, 'YTickLabel', OnRateRatioExponentAxis);
        view(-40,85)

        for iSeries = 1:numel(h)
            zData = get(h(iSeries), 'ZData');  % Get the z data
            index = logical(kron(isnan(zData(2:6:end, 2)), ones(6, 1)));  % Find empty bars
            zData(index, :) = nan;                 % Set the z data for empty bars to nan
            set(h(iSeries), 'ZData', zData);   % Update the graphics objects
        end
        saveas(gcf,[sprintf('.\\output%d\\pValues-kCatGap2Gef=%.4f_targetX=%.3f_', folderIndex, Log2KCatGap, TARGET_X_TO_MATCH_BY_ITERATION),'.fig'])
        saveas(gcf,[sprintf('.\\output%d\\pValues-kCatGap2Gef=%.4f_targetX=%.3f_', folderIndex, Log2KCatGap, TARGET_X_TO_MATCH_BY_ITERATION),'.png'])
        %     saveas(gcf,[sprintf('.\\output%d\\pValues-kOnRatio=%.2f_kOffGtp2kCat=%.2f_kOffRatio=%.2f_targetX=%.3f_kCatGap2Gef=%.4f', folderIndex, kOnRatioGtp2Gdp, kOffGdp2KCatRatio, kOffRatioGtp2Gdp, TARGET_X, kCatGap2Gef),'.pdf'])
        fclose(fid);
    end


    if SHOW_DISTRIB_TYPE_MATRIX == true

        fid2 = fopen(distribTypeFileName,'rt');
        C2 = textscan(fid2, '%f %f %f', 'Delimiter', ',');
        OnRateRatioExponentAxis = reshape(C2{1},nKOffs,[]); OnRateRatioExponentAxis = OnRateRatioExponentAxis(1,:);
        offRateLevelAxis = reshape(C2{2},nKOffs,[]); offRateLevelAxis = offRateLevelAxis(:,1);
        distribTypes = reshape(C2{3},nKOffs,[]);
        distribTypes = distribTypes';
        figure;
        distribTypes(distribTypes==0) = 0.001;
        h = bar3(distribTypes, 1);
        title(strTitle);
        lengthOnRateRatioArray = length(OnRateRatioExponentAxis);
        lengthOffRateLevelArray = length(offRateLevelAxis);
        pbaspect([ 1/lengthOnRateRatioArray 1/lengthOffRateLevelArray min(1/lengthOffRateLevelArray, 1/lengthOnRateRatioArray)]);
        customColorMap = [ones(1,3);cool(2000+1)];
        colormap(customColorMap);
        caxis([0 2]);
        view(2);
        for k = 1:length(h)
            zdata = h(k).ZData;
            h(k).CData = zdata;
            h(k).FaceColor = 'interp';
        end
        colorbar;
        xlabel('offRateLevel'), ylabel('OnRateRatioExponent'), zlabel('Type');
        if sum(isnan(OnRateRatioExponentAxis))
            OnRateRatioExponentAxis = listLog2_kOnGtp2KOnGdp;
        end
        if sum(isnan(offRateLevelAxis))
            offRateLevelAxis = listLog2_kOffGdp;
        end
        set(gca, 'XTickLabel',offRateLevelAxis, 'YTickLabel', OnRateRatioExponentAxis);

        for i = 1:numel(h)
            index = logical(kron(isnan(zData(2:6:end, 2)), ones(6, 1)));
            zData = get(h(i), 'ZData');
            zData(index, :) = nan;
            set(h(i), 'ZData', zData);
        end
        saveas(gcf,[sprintf('.\\output%d\\distribType-kCatGap2Gef=%.4f_targetX=%.3f_', folderIndex, Log2KCatGap, TARGET_X_TO_MATCH_BY_ITERATION),'.fig'])
        saveas(gcf,[sprintf('.\\output%d\\distribType-kCatGap2Gef=%.4f_targetX=%.3f_', folderIndex, Log2KCatGap, TARGET_X_TO_MATCH_BY_ITERATION),'.png'])

        %     saveas(gcf,[sprintf('.\\output%d\\distribType-kOnRatio=%.2f_kOffGtp2kCat=%.2f_kOffRatio=%.2f_targetX=%.3f_kCatGap2Gef=%.4f', folderIndex, kOnRatioGtp2Gdp, kOffGdp2KCatRatio, kOffRatioGtp2Gdp, TARGET_X, kCatGap2Gef),'.pdf'])
        fclose(fid2);
    end

end



function [SosRasGdp,SosRasGtp,RasGdp,kOnGdp] = getKOnGtp(kCatGap, kCatGef, kOnGtp2KOnGdp, kOffGdp, kOffGtp, nTargetRasGtp, nRasTotal)

    syms kOnGdp SosRasGdp SosRasGtp RasGdp real
    %eqn1=d(SosRasGdp)/dt
    %eqn2=d(SosRasGdp)/dt
    %eqn3=d(RasGtp)/dt
    %eqn4=d(RasGdp)/dt
    %eqn5=total Ras conserve
    eqn1 = kOnGdp * RasGdp - kOffGdp * SosRasGdp == 0;
    eqn2 = kOnGtp2KOnGdp * kOnGdp * nTargetRasGtp - kOffGtp * SosRasGtp == 0;
    eqn3 = kCatGef * (SosRasGdp + SosRasGtp) * RasGdp - kCatGap * nTargetRasGtp ...
        - kOnGtp2KOnGdp * kOnGdp * nTargetRasGtp + kOffGtp * SosRasGtp == 0;
    eqn4 = -kCatGef * (SosRasGdp + SosRasGtp) * RasGdp + kCatGap * nTargetRasGtp ...
        - kOnGdp * RasGdp + kOffGdp * SosRasGdp == 0;
    eqn5 =  SosRasGdp + SosRasGtp + RasGdp + nTargetRasGtp == nRasTotal;
    %Solve analytically
    [SosRasGdp,SosRasGtp,RasGdp,kOnGdp] = solve([eqn1,eqn2,eqn3,eqn4,eqn5, ...
        SosRasGdp >= 0, 100 > SosRasGtp, SosRasGtp >= 0,RasGdp >= 0,kOnGdp >= 0], ...
        [SosRasGdp,SosRasGtp,RasGdp,kOnGdp]);
    
    kOnGdp = double(kOnGdp);
end

function outNTargetRasGtp = getNTargetRasGtp(kCatGap, kCatGef, kOnGdp, kOnGtp2KOnGdp, kOffGdp, kOffGtp, nRasTotal)
    syms nTargetRasGtp SosRasGdp SosRasGtp RasGdp real
    %eqn1=d(SosRasGdp)/dt
    %eqn2=d(SosRasGdp)/dt
    %eqn3=d(RasGtp)/dt
    %eqn4=d(RasGdp)/dt
    %eqn5=total Ras conserve
    eqn1 = kOnGdp * RasGdp - kOffGdp * SosRasGdp == 0;
    eqn2 = kOnGtp2KOnGdp * kOnGdp * nTargetRasGtp - kOffGtp * SosRasGtp == 0;
    eqn3 = kCatGef * (SosRasGdp + SosRasGtp) * RasGdp - kCatGap * nTargetRasGtp ...
        - kOnGtp2KOnGdp * kOnGdp * nTargetRasGtp + kOffGtp * SosRasGtp == 0;
    eqn4 = -kCatGef * (SosRasGdp + SosRasGtp) * RasGdp + kCatGap * nTargetRasGtp ...
        - kOnGdp * RasGdp + kOffGdp * SosRasGdp == 0;
    eqn5 =  SosRasGdp + SosRasGtp + RasGdp + nTargetRasGtp == nRasTotal;
    %Solve analytically
    [nTargetRasGtp,~,~,~] = solve([eqn1,eqn2,eqn3,eqn4,eqn5, ...
        SosRasGdp >= 0, 100 > SosRasGtp, SosRasGtp >= 0,RasGdp >= 0,kOnGdp >= 0], ...
        [nTargetRasGtp, SosRasGdp,SosRasGtp,RasGdp]);
    
    outNTargetRasGtp = nTargetRasGtp;
end
% 
% function singleTimeTrace = getATimeTrace(flagForKOnAdj, nRasTotal, x_initial, duration, nMinDuration, kOnGtp, kOnGdp, kOffGtp, kOffGdp, kCatGef, kCatGap, gridlen)
% 
%     t = 0;
%     nRasGtp = nRasTotal * x_initial; % rasGTP
%     nRasGdp = nRasTotal * (1-x_initial);
%     nSosRasGtp = 0; % SosRasGtp
%     nSosRasGdp = 0; % SosRasGdp
% 
%     if flagForKOnAdj == true
%         nTrackedVars = 2; % (t, nRasGtp)
%     else
%         nTrackedVars = 5; % (t,nRasGtp,nRasGdp,nSosRasGtp,nSosRasGdp)
%     end
% 
%     % Memory allocation
%     singleTimeTrace = zeros(10 * round(duration), nTrackedVars); 
% 
%     % Writing the initial state
%     if flagForKOnAdj == true
%         singleTimeTrace(1, 1:end) = [t nRasGtp];
%     else
%         singleTimeTrace(1, 1:end) = [t,nRasGtp,nRasGdp,nSosRasGtp,nSosRasGdp]; 
%     end
%     recIdx = 0;
% 
%     nReactionNumber = 0;
%     while (t < duration) || (t < nMinDuration)
%         nReactionNumber = nReactionNumber + 1;
%         r1 = kOnGtp * nRasGtp;
%         r2 = kOnGdp * nRasGdp;
%         r3 = kOffGtp * nSosRasGtp;
%         r4 = kOffGdp * nSosRasGdp;
%         r5 = kCatGef * (nSosRasGtp + nSosRasGdp) / gridlen^2 * nRasGdp / gridlen^2;
%         r6 = kCatGap / gridlen^2 * nRasGtp / gridlen^2;
%         rT=r1+r2+r3+r4+r5+r6;
%         drawT = rand;
%         tau = 1/rT * log(1/drawT);
%         t = t + tau;
%         p1=r1/rT;
%         p2=r2/rT;
%         p3=r3/rT;
%         p4=r4/rT;
%         p5=r5/rT;
%         drawReaction = rand;
%         
%         if drawReaction < p1
%             nRasGtp = nRasGtp - 1;
%             nSosRasGtp = nSosRasGtp + 1;
%         elseif drawReaction < (p1+p2)
%             nSosRasGdp = nSosRasGdp + 1;
%             nRasGdp = nRasGdp - 1;
%         elseif drawReaction < (p1+p2+p3)
%             nSosRasGtp = nSosRasGtp - 1;
%             nRasGtp = nRasGtp + 1;
%         elseif drawReaction < (p1+p2+p3+p4)
%             nSosRasGdp = nSosRasGdp - 1;
%             nRasGdp = nRasGdp + 1;
%         elseif drawReaction < (p1+p2+p3+p4+p5)
%             nRasGtp=nRasGtp+1;
%             nRasGdp=nRasGdp-1;
%         else
%             nRasGtp=nRasGtp-1;
%             nRasGdp=nRasGdp+1;
%         end
%     
%         if floor(t*10) > recIdx % recording every 0.1 second
%             if flagForKOnAdj == true
%                 singleTimeTrace(recIdx+2, 1:nTrackedVars) = [t nRasGtp];
%             else
%                 singleTimeTrace(recIdx+2, 1:nTrackedVars) = [t,nRasGtp,nRasGdp,nSosRasGtp,nSosRasGdp]; 
%             end
%             recIdx = recIdx + 1;
%         end
%     end
%     disp(nReactionNumber);
% end

function traceEverySecond = reorganizeForEverySecond(recordedTimeArray, inputTimeTrace, dataLength, nCorrals)
    traceEverySecond = zeros(dataLength,nCorrals);
    for j = 0:1:dataLength % record for every (integer) second
        for i = 1:nCorrals
            jSecondsMinusRecordedTimes = j - recordedTimeArray(:,i);
            jSecondsMinusRecordedTimes(jSecondsMinusRecordedTimes<=0) = nan;
            [~, idxes] = min(jSecondsMinusRecordedTimes);
            traceEverySecond(j+1,i) = inputTimeTrace(idxes,i);
        end
    end
end
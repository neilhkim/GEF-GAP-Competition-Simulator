clear all
% Define the maximum value for the concentration titration.
max = 11;

% Loop through a range of values
for i = 1:max
    % Calculate k3 and k1 (binding rates, which correspond to concentrations) for each iteration, incrementally increasing them.
    k3 = 0.01 * (i-1) * 0.00003 + 0.0000001;
    k1 = 0.01 * (i-1) * 0.00003 + 0.0000001;

    % Set fixed values for k4 and k2.
    k4 = 0.01 * 0.1;
    k2 = 0.01 * 0.5;

    % Define catalytic rate constants
    kcat1 = 0.01;
    kcat2 = 0.01*0.5;

    % Number of total Ras molecules
    Rtotal = 1000;

    % Declare symbolic variables for species concentrations.
    syms SRd SRt Rd Rt

    % Equations representing the kinetics at equilibrium.
    % eqn1 = dSRd / dt = 0
    eqn1 = k1*Rd - k2*SRd == 0;
    % eqn1 = dSRt / dt = 0
    eqn2 = k3*Rt - k4*SRt == 0;
    % eqn3 = dRt / dt = 0
    eqn3 = kcat1*(SRd+SRt)*Rd - kcat2*Rt - k3*Rt + k4*SRt == 0;
    % eqn3 = dRd / dt = 0
    eqn4 = -kcat1*(SRd+SRt)*Rd + kcat2*Rt - k1*Rd + k2*SRd == 0;
    % eqn5 = total Ras = Rtotal
    eqn5 =  SRd + SRt + Rd + Rt == Rtotal;

    % Solve analytically
    [SRd,SRt,Rd,Rt] = solve([eqn1,eqn2,eqn3,eqn4,eqn5,SRd>=0,SRt>=0,Rd>=0,Rt>=0],[SRd,SRt,Rd,Rt]);

    % Store solutions in arrays, each column corresponding to an iteration.
    SRd_value(1:length(SRd),i) = double(SRd);
    SRt_value(1:length(SRt),i) = double(SRt);
    Rd_value(1:length(Rd),i) = double(Rd);
    Rt_value(1:length(Rt),i) = double(Rt);
end

% Save the results to .mat files

% Create output directory if it doesn't exist
folderIndex = 1; % Start with index 1 for naming output folders
while exist(sprintf('.\\output%d', folderIndex), 'dir')
    folderIndex = folderIndex + 1; % Increment index if folder already exists
end
mkdir(sprintf('.\\output%d', folderIndex)); % Create the new output directory

foldername = sprintf('.\\output%d\\', folderIndex);
save([foldername, 'SRd_value.mat'], 'SRd_value');
save([foldername, 'SRt_value.mat'], 'SRt_value');
save([foldername, 'Rd_value.mat'], 'Rd_value');
save([foldername, 'Rt_value.mat'], 'Rt_value');

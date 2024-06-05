% Script to create a plot comparing results of Moco simulations with different mesh intervals and tolerances
% Kristen Steudel, May 29, 2024

% Close all open figures, clear workspace, and clear command window
close all; clear all; clc;

import org.opensim.modeling.*;

% Define paths to the result files
subjectID = 'Sub003';
rep = '01';
baseDir = 'G:\My Drive\Aa 203\ProjectResults\';
addpath('G:\Shared drives\HPL_Drive\Nord Sprinting Study\Data Processing Scripts\mocap-processing\common')

% Specify the different configurations you tested
meshIntervals = [0.2, 0.15, 0.1, 0.05, 0.02, 0.01]; % example mesh intervals
tolerances = [1e-3, 1e-3]; % example tolerances

files = {'Sub003Nord03_MocoInverse_solution_point2.sto', 'Sub003Nord03_MocoInverse_solution_point15.sto',...
    'Sub003Nord03_MocoInverse_solution_point1.sto','Sub003Nord03_MocoInverse_solution_point05.sto', ...
    'Sub003Nord03_MocoInverse_solution_point02.sto', 'Sub003Nord03_MocoInverse_solution_point01.sto'};
solver_durations = [];
objectives = [];

for j = 1:length(files)
    moco = strcat(baseDir, files{j});
    disp(moco);
    % Read the header to extract solver_duration and objective
    fileID = fopen(moco, 'r');
    headerLine = fgetl(fileID);
    while ischar(headerLine) && ~startsWith(headerLine, '%')
        if contains(headerLine, 'solver_duration')
            solver_duration = sscanf(headerLine, 'solver_duration = %f');
            solver_durations = [solver_durations; solver_duration];
        end
        if contains(headerLine, 'objective')
            objective = sscanf(headerLine, 'objective = %f');
            objectives = [objectives; objective];
        end
        headerLine = fgetl(fileID);
    end
    fclose(fileID);
    
    % Load the .sto file
    stoFile = TimeSeriesTable(moco);
    % Get column labels
    columnLabels = stoFile.getColumnLabels();
    % Get the numeric data as a MATLAB matrix
    data = stoFile.getMatrix();
    % make the data a mat matrix
    data_matrix = data.getAsMat;
    
    % Loop through data matrix to figure out index of each muscle
    for i = 0:size(data_matrix,2)-1
        if strcmp(columnLabels.get(i),'/forceset/addbrev_r/activation')
        display(i, 'addbrev_r') 
        end
        if strcmp(columnLabels.get(i),'/forceset/recfem_r/activation')
        display(i, 'recfem_r')
        end
        if strcmp(columnLabels.get(i),'/forceset/vaslat_r/activation')
        display(i, 'vaslat_r')
        end
        if strcmp(columnLabels.get(i),'/forceset/vasmed_r/activation')
        display(i, 'vasmed_r')
        end
        if strcmp(columnLabels.get(i),'/forceset/tibant_r/activation')
        display(i, 'tibant_r')
        end
        if strcmp(columnLabels.get(i),'/forceset/bflh_r/activation')
        display(i, 'bflh_r')
        end
        if strcmp(columnLabels.get(i),'/forceset/semiten_r/activation')
        display(i, 'semiten_r')
        end
        if strcmp(columnLabels.get(i),'/forceset/gaslat_r/activation')
        display(i, 'gaslat_r')
        end
        if strcmp(columnLabels.get(i),'/forceset/gasmed_r/activation')
        display(i, 'gasmed_r')
        end
        if strcmp(columnLabels.get(i),'/forceset/soleus_r/activation')
        display(i, 'soleus_r')
        end
        if strcmp(columnLabels.get(i),'/forceset/vaslat_l/activation')
        display(i, 'vaslat_l')
        end
        if strcmp(columnLabels.get(i),'/forceset/bflh_l/activation')
        display(i, 'bflh_l')
        end
        if strcmp(columnLabels.get(i),'/forceset/gasmed_l/activation')
        display(i, 'gasmed_l')
        end
        if strcmp(columnLabels.get(i),'/forceset/soleus_l/activation')
        display(i, 'soleus_l')
        end
        % /forceset/semimem_l/activation\
        if strcmp(columnLabels.get(i),'/forceset/semimem_l/activation')
        display(i, 'semimembranosis_l')
        end
    
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Now access each element in the column, and make it into a vector
    
    % Currently manually making the moco time vector, find a way to extract
    % this from the .sto file to make life easier
    % t_moco = 5:(3/length(data_matrix)):8;
    % t_moco = t_moco(1:end-1);
    % t_moco = t_moco';
    % 
    SM_moco = data_matrix(:, 71);

    RF_moco = data_matrix(:,29);
    
    
    VL_moco = data_matrix(:,38);
    
    
    VM_moco = data_matrix(:,39);
    
      
    TA_moco = data_matrix(:,35);
    
    
    BF_moco = data_matrix(:,6);
    
    
    ST_moco = data_matrix(:,32);
    
    
    LG_moco = data_matrix(:,12);
    
    
    MG_moco = data_matrix(:,13);
    
    
    SOL_moco = data_matrix(:,33);
    
    
    VLcont_moco = data_matrix(:,78);
    
    
    BFcont_moco = data_matrix(:,46);
    
    
    MGcont_moco = data_matrix(:,53);
    
    
    SOLcont_moco = data_matrix(:,73);


    commonLength = 100; 

    % Interpolate BF_moco to common length
    BF_moco_interp = interp1(linspace(0, 1, length(BF_moco)), BF_moco, linspace(0, 1, commonLength));

    % Interpolate ST_moco to common length
    ST_moco_interp = interp1(linspace(0, 1, length(ST_moco)), ST_moco, linspace(0, 1, commonLength));

    % Interpolate ST_moco to common length
    SM_moco_interp = interp1(linspace(0, 1, length(SM_moco)), SM_moco, linspace(0, 1, commonLength));

    % Plot the biceps femoris activation for each file
    figure(3);
    hold on;
    colors = lines(length(meshIntervals) * length(tolerances));
    legendEntries = {'0.2', '0.15', '0.1', '0.05', '0.02', '0.01'};
    plot(linspace(0,1, commonLength), BF_moco_interp);
    title('Biceps Femoris Activation Over Trajectory');
    xlabel('Percent of Trajectory');
    ylabel('Biceps Femoris Activation');
    legend(legendEntries, 'Location', 'Best');
    hold off;

    % Plot the semitendinosis activation for each file
    figure(4);
    hold on;
    colors = lines(length(meshIntervals) * length(tolerances));
    legendEntries = {'0.2', '0.15', '0.1', '0.05', '0.02', '0.01'};
    plot(linspace(0,1, commonLength), ST_moco_interp);
    title('Semitendinosis Activation Over Trajectory');
    xlabel('Percent of Trajectory');
    ylabel('Semitendinosis Activation');
    legend(legendEntries, 'Location', 'Best');
    hold off;

    % Plot the semimembranosis activation for each file
    figure(5);
    hold on;
    colors = lines(length(meshIntervals) * length(tolerances));
    legendEntries = {'0.2', '0.15', '0.1', '0.05', '0.02', '0.01'};
    plot(linspace(0,1, commonLength), SM_moco_interp);
    title('Semimembranosis Activation Over Trajectory');
    xlabel('Percent of Trajectory');
    ylabel('Semimembranosis Activation');
    legend(legendEntries, 'Location', 'Best');
    hold off;

end 

% Save the plot
saveas(gcf, strcat(baseDir, 'Nord_BicepsFemorisActivations.png'));

saveas(gcf, strcat(baseDir, 'Nord_SemitendinosisActivations.png'));

figure(1)
sz = 120;
scatter(meshIntervals, objectives, sz,'MarkerEdgeColor',[0 .5 .5],...
              'MarkerFaceColor',[0 .7 .7],...
              'LineWidth',1.5);
xlabel('Mesh Interval');
ylabel('Objective Value');
grid on
saveas(gcf, strcat(baseDir, 'Nord_Moco_Objectives.png'));


figure(2)
sz = 120;
scatter(meshIntervals, solver_durations, sz,'MarkerEdgeColor',[0 .5 .5],...
              'MarkerFaceColor',[0 .7 .7],...
              'LineWidth',1.5);
xlabel('Mesh Interval');
ylabel('Solver Duration (seconds)')
grid on
saveas(gcf, strcat(baseDir, 'Nord_Moco_Durations.png'));





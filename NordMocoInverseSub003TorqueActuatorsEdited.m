% -------------------------------------------------------------------------- %
% OpenSim Moco: NordMocoInverseTorqueActuators.m      
% Kristen Steudel May 29, 2024
% -------------------------------------------------------------------------- %
close all; clear all; clc;
import org.opensim.modeling.*;
subjectID = 'Sub003';
rep = '03';
addpath (strcat('G:\Shared drives\HPL_Drive\Nord Sprinting Study\',subjectID,'\AddBiomechanics\'));

% Define the file paths
kinematicsFile = strcat('G:\Shared drives\HPL_Drive\Nord Sprinting Study\',subjectID,'\AddBiomechanics5_15\IK\ProcessedIK\Nord_03.mot');
modelFile = strcat('G:\Shared drives\HPL_Drive\Nord Sprinting Study\', subjectID, '\AddBiomechanics5_15\Models\match_markers_but_ignore_physics.osim');
externalLoadsFile = strcat('G:\Shared drives\HPL_Drive\Nord Sprinting Study\', subjectID, '\CleanedData\ID\ID_ExternalLoads_Nordics_Rep_3.xml');


% Check if model file exists
if exist(modelFile, 'file') == 2
    tempModel = ModelProcessor(modelFile);
else
    error('Model file does not exist: %s', modelFile);
end

% Ensure external loads file exists
if exist(externalLoadsFile, 'file') == 2
    tempModel.append(ModOpAddExternalLoads(externalLoadsFile));
else
    error('External loads file does not exist: %s', externalLoadsFile);
end

% Remove muscles
tempModel.append(ModOpRemoveMuscles());

% Optional: Add reserve actuators
%tempModel.append(ModOpAddReserves(1.0, 200.0));

myModel = tempModel.process();

% Add reserves with specified optimal force (OF) for each DOF
coordNames = {'pelvis_tilt','pelvis_list','pelvis_rotation',...
    'pelvis_tx','pelvis_ty','pelvis_tz','hip_flexion_r','hip_adduction_r',...
    'hip_rotation_r','knee_angle_r','ankle_angle_r','subtalar_angle_r','mtp_angle_r',...
    'hip_flexion_l','hip_adduction_l','hip_rotation_l','knee_angle_l','ankle_angle_l',...
    'subtalar_angle_l','mtp_angle_l'};

coordControl = [1 1 1 1 1 1,...
    1 1 1 1 1 1 1 1 1 1 1 1 1 1];
coordForce = [2000 2000 2000 2000 2000 2000,...
    2000 2000 2000 2000 2000 2000 2000 2000 2000 2000 2000 2000 2000 2000];

for i = 1:length(coordNames)
    addCoordinateActuator(myModel, coordNames(i), coordForce(i), coordControl(i));
end

myModel.finalizeConnections();

modelProcessor = ModelProcessor(myModel);

inverse = MocoInverse();

% Checks for the kinematics file existence
if exist(kinematicsFile, 'file') == 2
    % File exists, proceed with processing
    inverse.setKinematics(TableProcessor(kinematicsFile));
else
    error('Kinematics file does not exist: %s', kinematicsFile);
end

inverse.setModel(modelProcessor);   %Takes modelProcessor as input type


% Set the kinematics data for the inverse problem
inverse.setKinematics(TableProcessor(kinematicsFile));

% Initial time, final time, and mesh interval
inverse.set_initial_time(2.995);
inverse.set_final_time(8.685);
inverse.set_mesh_interval(0.05);
inverse.set_convergence_tolerance(1e-5);
inverse.set_constraint_tolerance(1e-5);
inverse.set_max_iterations(4000);

% By default, Moco gives an error if the kinematics contains extra columns.
inverse.set_kinematics_allow_extra_columns(true);

% Manually define initial guess
% Create a MocoTrajectory object with the appropriate structure
study = inverse.initialize();
problem = study.updProblem();

excitegoal = MocoControlGoal.safeDownCast(problem.updGoal('excitation_effort'));

% % Option to add in a solution for an initial guess
% guess = MocoTrajectory(study);
% 
% % Define your initial guess values for states and controls
% % For example, setting all states to zero and all controls to a small value
% initial_state_value = 0.0;
% initial_control_value = 0.01;
% 
% % Get the names of the states and controls
% state_names = guess.getStateNames();
% control_names = guess.getControlNames();
% 
% % Set the initial guess for the states
% for i = 0:state_names.size() - 1
%     guess.setState(state_names.get(i), initial_state_value);
% end
% 
% % Set the initial guess for the controls
% for i = 0:control_names.size() - 1
%     guess.setControl(control_names.get(i), initial_control_value);
% end
% 
% % Set the manually defined initial guess for the MocoInverse problem
% inverse.setInitialGuess(guess);


forceSet = myModel.getForceSet();
for i = 0:forceSet.getSize()-1
    forcePath = forceSet.get(i).getAbsolutePathString();
    if contains(string(forcePath), 'pelvis')
        excitegoal.setWeightForControl(forcePath,0); % set the pelvis residuals controls to 0
    end
end

% Nick's current work around to the Inverse Problem's positions/speeds
% issue given that we are 'customizing' the Inverse Problem
problem = study.getProblem();
model = problem.getPhase(0).getModel();
model.initSystem();
position_motion = PositionMotion.safeDownCast(model.getComponent('position_motion'));
coordinatesTable = TimeSeriesTable(kinematicsFile);
kinematics = position_motion.exportToTable(coordinatesTable.getIndependentColumn());

% kinematicsTable = readMOT(kinematicsFile);
% disp(kinematicsTable.labels)

%to resolve insertStatesTrajectory

% Solve the problem and write the solution to a Storage file.
solution = study.solve();
solution.unseal();
solution.insertStatesTrajectory(kinematics); %Got an error saying unable

solutionFilePath = fullfile('G:', 'Shared drives', 'HPL_Drive', 'Nord Sprinting Study', subjectID, [subjectID, 'Nord', rep, '_MocoInverse_solution.sto']);
solution.write(solutionFilePath);

% Check if the solution file exists before generating the report
if exist(solutionFilePath, 'file') ~= 2
    error('Solution file does not exist: %s', solutionFilePath);
end

% Generate a report with plots for the solution trajectory.
model = modelProcessor.process();
% report = osimMocoTrajectoryReport(model, ...
%         strcat('G:\Shared drives\HPL_Drive\Nord Sprinting Study\',subjectID,'\Nord',rep,'_MocoInverse_solution_TighterTolerance_MeshInterval0.05.sto'), 'bilateral', true); % Sends the solution report to the individual's folder
report = osimMocoTrajectoryReport(model, solutionFilePath, 'bilateral', true);

% The report is saved to the working directory.
reportFilepath = report.generate();
open(reportFilepath);



function addCoordinateActuator(modelName, coordName, optForce, controlLim)

import org.opensim.modeling.*;
coordSet = modelName.updCoordinateSet();
actu = CoordinateActuator();
actu.setName(strcat('reserve_jointset_',coordName));
actu.setCoordinate(coordSet.get(coordName));
actu.setOptimalForce(optForce);
actu.setMaxControl(controlLim);
actu.setMinControl(-controlLim);
% Add to ForceSet
modelName.addForce(actu);

end

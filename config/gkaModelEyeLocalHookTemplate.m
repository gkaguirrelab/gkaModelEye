function gkaModelEyeLocalHook
% 
% For use with the ToolboxToolbox.  If you copy this into your
% ToolboxToolbox localToolboxHooks directory (by default,
% ~/localToolboxHooks) and delete "LocalHooksTemplate" from the filename,
% this will get run when you execute tbUse({'gkaModelEye'}) to set up for
% this project.  You then edit your local copy to match your local machine.
%
% The main thing that this does is define Matlab preferences that specify input and output
% directories.
%
% You will need to edit the project location and i/o directory locations
% to match what is true on your computer.

 
%% Define project
toolboxName = 'gkaModelEye';
 
%% Clear out old preferences
if (ispref(toolboxName))
    rmpref(toolboxName);
end

%% Check for required Matlab toolboxes
% The set of Matlab add-on toolboxes being used can be determined by
% running the ExampleTest code, followed by the license function.
%{
    RunExamples(fullfile(userpath(),'toolboxes','gkaModelEye'));
    license('inuse')
%}
% This provides a list of toolbox license names. In the following
% assignment, the license name is given in the comment string after the
% matching version name for each toolbox.
requiredAddOns = {...
    'Mapping Toolbox',...                       % map_toolbox
    'Curve Fitting Toolbox',...                 % curve_fitting_toolbox
    'Optimization Toolbox',...                  % optimization_toolbox
    'Robotics System Toolbox',...               % robotics_system_toolbox
    'Statistics and Machine Learning Toolbox',...   % statistics_toolbox
    'Symbolic Math Toolbox'...                  % symbolic_toolbox
    };
% Given this hard-coded list of add-on toolboxes, we then check for the
% presence of each and issue a warning if absent.
V = ver;
VName = {V.Name};
warnState = warning();
warning off backtrace
for ii=1:length(requiredAddOns)
    if ~any(strcmp(VName, requiredAddOns{ii}))
        warnString = ['The Matlab ' requiredAddOns{ii} ' is missing. ' toolboxName ' may not function properly.'];
        warning('localHook:requiredMatlabToolboxCheck',warnString);
    end
end
warning(warnState);

end

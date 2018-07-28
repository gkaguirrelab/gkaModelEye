function eye = modelEyeParameters( varargin )
% Return the parameters of a model eye
%
% Syntax:
%  eye = modelEyeParameters()
%
% Description:
%   This routine returns the parameters of a model eye used in the
%   sceneGeometry routines.
%
%   The parameters returned by this routine correspond to the eyeWorld
%   coordinate space used in pupilProjection_fwd, which is relative to the
%   optical / pupillary axis, with the apex of the cornea set as zero in
%   depth. The space has the dimensions [depth, horizontal, vertical];
%   negative values of depth are towards the back of the eye. The model
%   assumes the optical and pupil axis of the eye are algined.
%
% Inputs:
%   none
%
% Optional key/value pairs:
%  'sphericalAmetropia'   - Scalar, in units of diopters. Several
%                           parameters of the model eye are adjusted by the
%                           spherical refractive error of the eye. A
%                           negative number is the correction that would be
%                           used for a myopic person.
%  'axialLength'          - Scalar. This is the axial length along the 
%                           optical axis. When set, this fixes the axial
%                           length of the eye to the passed value in
%                           millimeters. As the modeled anterior chamber
%                           depth is not variable, this change is enforced
%                           on the vitreous chamber. The remaining
%                           dimensions of the vitreous chamber are scaled
%                           to fit the proportions predicted by the
%                           Atchison model for the specified degree of
%                           ametropia.
%  'eyeLaterality'        - A text string that specifies which eye (left,
%                           right) to model. Allowed values (in any case)
%                           are {'left','right','L','R','OS','OD'}
%  'species'              - A text string that specifies the species to be
%                           modeled. Supported values (in any case) are
%                           {'human','dog'}
%  'spectralDomain'       - String, options include {'vis','nir'}.
%                           This is the wavelength domain within which
%                           imaging is being performed. The refractive
%                           indices vary based upon this choice.
%  'visualAxisDegRetina'  - 1x3 vector. This is the position of the fovea 
%                           w.r.t. to optical axis in degrees of retina.
%                           The values are [azimuth, elevation, torsion].
%                           Used in model development.
%  'opticDiscAxisDegRegina'  - 1x3 vector. This is the position of the  
%                           optic disc w.r.t. to optical axis in degrees of
%                           retina. The values are [azimuth, elevation,
%                           torsion]. Used in model development.
%
% Outputs:
%   eye                   - A structure with fields that contain the values
%                           for the model eye.
%
% Examples:
%{
    % Default parameters, corresponding to an emmetropic, right, human eye
    eye = modelEyeParameters();
%}
%{
    % Parameters for an myopic (-3), left, human eye
    eye = modelEyeParameters('sphericalAmetropia',-3,'eyeLaterality','left');
%}


%% input parser
p = inputParser; p.KeepUnmatched = true;

% Optional
p.addParameter('sphericalAmetropia',0,@isscalar);
p.addParameter('axialLength',[],@(x)(isempty(x) || isscalar(x)));
p.addParameter('cornealAxis',[],@(x)(isempty(x) || isnumeric(x)));
p.addParameter('eyeLaterality','Right',@ischar);
p.addParameter('species','Human',@ischar);
p.addParameter('spectralDomain','nir',@ischar);
p.addParameter('skipEyeAxes',false,@islogical);

% parse
p.parse(varargin{:})

% Interpret the passed laterality
switch p.Results.eyeLaterality
    case {'right','RIGHT','Right','R','r','od','OD'}
        eyeLaterality = 'Right';
    case {'left','LEFT','Left','L','l','os','OS'}
        eyeLaterality = 'Left';
    otherwise
        error('Please specify a valid eye laterality for the model eye');
end

% Create an empty eye struct
eye = struct();

% Meta data regarding the units of the model
eye.meta.p = p.Results;
eye.meta.units = 'mm';
eye.meta.coordinates = 'eyeWorld';
eye.meta.dimensions = {'depth (axial)' 'horizontal' 'vertical'};
eye.meta.eyeLaterality = eyeLaterality;
eye.meta.sphericalAmetropia = p.Results.sphericalAmetropia;
eye.meta.species = p.Results.species;
eye.meta.spectralDomain = p.Results.spectralDomain;

% Switch parameters at the top level by species
switch eye.meta.species

    %% Human
    case {'human','Human','HUMAN'}
                
        % Eye anatomy
        eye.cornea = human.cornea(eye, p.Results.cornealAxis);
        eye.iris = human.iris(eye);
        eye.pupil = human.pupil(eye);
        eye.retina = human.retina(eye);
        eye.lens = human.lens(eye);
        
        % Rotation centers
        eye.rotationCenters = human.rotationCenters(eye);
        
        % Computation of the eye axes is optional
        if ~p.Results.skipEyeAxes
           eye.axes = human.axes(eye);
        end
        
        % Refractive indices
        eye.index.vitreous = returnRefractiveIndex( 'vitreous', p.Results.spectralDomain );
        eye.index.aqueous = returnRefractiveIndex( 'aqueous', p.Results.spectralDomain );

        
    %% Canine
    case {'dog','Dog','canine','Canine'}


        
    otherwise
        error('Please specify a valid species for the eye model');
end



end % function



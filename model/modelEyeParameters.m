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
%   assumes the optical and pupil axes of the eye are algined.
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
%                           optical axis. This value is converted into an
%                           equivalent spherical error and then used to set
%                           the eye biometry. If both sphericalAmetropia
%                           and axialLength are passed, then inly spherical
%                           error will be used.
%  'eyeLaterality'        - A text string that specifies which eye (left,
%                           right) to model. Allowed values (in any case)
%                           are {'left','right','L','R','OS','OD'}
%  'species'              - A text string that specifies the species to be
%                           modeled. Supported values (in any case) are
%                           {'human','dog'}
%  'ageYears'             - Scalar that supplies the age in years of the
%                           eye to be modeled. Modifies the lens.
%  'accommodationDiopeters' - Scalar that supplies the accommodation state
%                           of the eye. Valid values range from zero
%                           (unaccommodated) to +10.
%  'measuredCornealCurvature' - 1x2 or 1x3 vector. Provides the horizontal 
%                           and vertical curvature of the cornea (diopters;
%                           K1 and K2). The third value is the rotation of
%                           the "horizontal" axis away from horizontal in
%                           degrees to allow specification of oblique
%                           asigmatism. Only the modeling of regular
%                           corneal astigmatism is supported. In subsequent
%                           routines the curvature in diopters is converted
%                           to a radius of curvature in mm using:
%                               r = 1000*(1.3375-1)/K
%                           If left undefined, the "canonical" Navarro 2006
%                           corneal parameters will be used.
%  'spectralDomain'       - String, options include {'vis','nir'}.
%                           This is the wavelength domain within which
%                           imaging is being performed. The refractive
%                           indices vary based upon this choice.
%  'skipEyeAxes'          - Logical. If set to true, the computation of
%                           retinal landmarks and eye axes is skipped. This
%                           is used internally by routines during code
%                           complilation.
%  'skipNodalPoint'          - Logical. If set to true, the computation of
%                           the effective nodal point is skipped. This is
%                           used internally by routines during code
%                           complilation.
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
p.addParameter('sphericalAmetropia',[],@(x)(isempty(x) || isscalar(x)));
p.addParameter('axialLength',[],@(x)(isempty(x) || isscalar(x)));
p.addParameter('eyeLaterality','Right',@ischar);
p.addParameter('species','Human',@ischar);
p.addParameter('ageYears',20,@isscalar);
p.addParameter('accommodationDiopeters',0,@isscalar);
p.addParameter('measuredCornealCurvature',[],@(x)(isempty(x) || isnumeric(x)));
p.addParameter('spectralDomain','nir',@ischar);
p.addParameter('skipEyeAxes',false,@islogical);
p.addParameter('skipNodalPoint',false,@islogical);

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

% Spherical refractive error (ametropia) is strongly correlated with axial
% length. The model uses spherical ametropia to determine several
% biometric properties of the eye. If axial length is provided instead of
% spherical error, then Eq 9 of Atchison (2006) Vision Research is used to
% assign a spherical ametropia. If both ametropia and axial length are
% provided, the curvature of the cornea (if not specified) and the
% curvature of the retina are specified by spherical ametropia. The
% absolute size of the vitreous chamber, however, is scaled up or down to
% force the total axial length to be equal to the measured value.
if isempty(p.Results.sphericalAmetropia) && isempty(p.Results.axialLength)
    sphericalAmetropia = 0;
    notes = 'Spherical refractive error not set; assuming emmetropia';
end
if isempty(p.Results.sphericalAmetropia) && ~isempty(p.Results.axialLength)    
    ametropiaFromLength = @(x) (23.58 - x)./0.299; % Atchison 2006 Eq 9.
    sphericalAmetropia = ametropiaFromLength(p.Results.axialLength);
    notes = 'Spherical refractive error derived from axial length';
end
if ~isempty(p.Results.sphericalAmetropia) && isempty(p.Results.axialLength)
    sphericalAmetropia = p.Results.sphericalAmetropia;
    notes = 'Spherical refractive error provided; axial length derived from SR';
end
if ~isempty(p.Results.sphericalAmetropia) && ~isempty(p.Results.axialLength)
    sphericalAmetropia = p.Results.sphericalAmetropia;
    notes = 'Spherical refractive error and axial length provided. Retinal curvature set by SR, absolute size by axial length';
end

% Meta data regarding properties of the model
eye.meta.p = p.Results;
eye.meta.units = 'mm';
eye.meta.coordinates = 'eyeWorld';
eye.meta.dimensions = {'depth (axial)' 'horizontal' 'vertical'};
eye.meta.eyeLaterality = eyeLaterality;
eye.meta.sphericalAmetropia = sphericalAmetropia;
eye.meta.axialLength = p.Results.axialLength;
eye.meta.species = p.Results.species;
eye.meta.ageYears = p.Results.ageYears;
eye.meta.accommodationDiopeters = p.Results.accommodationDiopeters;
eye.meta.measuredCornealCurvature = p.Results.measuredCornealCurvature;
eye.meta.spectralDomain = p.Results.spectralDomain;
eye.meta.notes = notes;

% Switch parameters at the top level by species
switch eye.meta.species

    %% Human
    case {'human','Human','HUMAN'}
                
        % Eye anatomy
        eye.cornea = human.cornea(eye);
        eye.iris = human.iris(eye);
        eye.pupil = human.pupil(eye);
        eye.retina = human.retina(eye);
        eye.lens = human.lens(eye);
        
        % If the axial length was not passed, calculate and store the value
        if isempty(eye.meta.axialLength)
            retinaRadii = quadric.radii(eye.retina.S);
            retinaCenter = quadric.center(eye.retina.S);
            eye.meta.axialLength = retinaRadii(1)-retinaCenter(1);
        end

        % Rotation centers
        eye.rotationCenters = human.rotationCenters(eye);
        
        % Identification of the eye axes is optional
        if ~p.Results.skipEyeAxes
           eye.axes = human.axes(eye);
        end
        
        % Refractive indices
        eye.index.vitreous = returnRefractiveIndex( 'vitreous', p.Results.spectralDomain );
        eye.index.aqueous = returnRefractiveIndex( 'aqueous', p.Results.spectralDomain );
         
        % Identification of the effective nodal point is optional
        if ~p.Results.skipNodalPoint
           eye.meta.nodalPoint = calcEffectiveNodalPoint(eye);
        end
        
    %% Canine
    case {'dog','Dog','canine','Canine'}
        error('Geoff needs to implement the canine model here');
        
        
    otherwise
        error('Please specify a valid species for the eye model');
end



end % function



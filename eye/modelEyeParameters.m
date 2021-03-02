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
%   coordinate space used in projectModelEye, which is relative to the
%   optical axis, with the apex of the cornea set as zero in depth. The
%   space has the dimensions [depth, horizontal, vertical]; negative values
%   of depth are towards the back of the eye.
%
%   Note that not all parameters are applicable to all "species" of model
%   eye. For example, the Drasdo & Fowler 1974 model has fixed size and
%   accommodation, so most of these optional key-value parameters are
%   ignored.
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
%  'accommodation'        - Scalar. The accommodative state of the eye, in
%                           diopters. If set, the model will search for the
%                           navarroD parameter that produces the requested
%                           accommodation for points on the longitudinal
%                           axis. Not all model eyes are capable of
%                           accommodating at all distances. For example, a
%                           myopic eye cannot be brought into focus at far
%                           distances.
%  'navarroD'             - Scalar. A parameter of the model that
%                           influences the shape of the crystalline lens.
%                           This parameter is passed by internal functions,
%                           and is typically not set by the user.
%  'axialLength'          - Scalar. This is the axial length (in mm) along
%                           the optical axis. This value is converted into
%                           an equivalent spherical error and then used to
%                           set the eye biometry. If both
%                           sphericalAmetropia and axialLength are passed,
%                           then only spherical error will be used.
%  'eyeLaterality'        - A text string that specifies which eye (left,
%                           right) to model. Allowed values (in any case)
%                           are {'left','right','L','R','OS','OD'}
%  'species'              - A text string that specifies the species to be
%                           modeled. Supported values (in any case) are
%                           {'human','dog'}
%  'ageYears'             - Scalar that supplies the age in years of the
%                           eye to be modeled. Influences the refractive
%                           index values of the lens.
%  'derivedParams'        - Struct that contains fields with parameters
%                           used by various eye model components. If left
%                           empty, the parameters will be obtained from a
%                           stored set of values.
%  'kvals' -                1x2 to 1x5 vector. Provides the horizontal
%                           and vertical curvature of the cornea (diopters;
%                           K1 and K2). The first value is always the
%                           smaller, and thus describes the "flatter"
%                           surface of the retina. The second value is
%                           always larger, and describes the curvature of
%                           the surface of the retina oriented 90 degrees
%                           away from the flatest surface. In subsequent
%                           routines the curvature in diopters is converted
%                           to a radius of curvature in mm using:
%                               r = 1000*(1.3375-1)/K
%                           The additional parameter values give the
%                           rotation of the corneal ellipsoid. In the order
%                           of:
%                               [torsion, tilt, tip]
%                           which are (respectively), the rotation of the
%                           "horizontal" axis away from horizontal in
%                           degrees (a.k.a., oblique asigmatism; only the
%                           modeling of regular corneal astigmatism is
%                           supported), then the rotation of the cornea
%                           around the vertical axis (tilt) and around the
%                           horizontal axis (tip). If left undefined, the
%                           "canonical" Navarro 2006 corneal parameters
%                           will be used.
%  'corneaAxialRadius'    - Scalar. The length (in mm) of the axial
%                           semi-radius of the corneal ellipsoid. If not
%                           left empty, this value will override the
%                           default value implemented within the cornea
%                           function.
%  'rotationCenterScalers' - 1x2 numeric vector. These values apply joint
%                           and differential multiplicative scaling on
%                           the positions of the azimuthal and elevational
%                           rotation centers of the eye.
%  'primaryPosition'      - 1x2 numeric vector. This is the [azi ele] pose
%                           of the eye for which torsion is zero, and for
%                           which every eye movement consistent with
%                           Listing's Law will result in an eye pose for
%                           which torsion is also zero. This parameter has
%                           the odd property that it is not intrinsic to
%                           the eye, but is defined within the coordinate
%                           frame in which the alignment of the camera and
%                           eye optical systems has eyePose values of
%                           [0, 0].
%  'spectralDomain'       - String or numerica scalar. This is the
%                           wavelength domain within which imaging is being
%                           performed. The refractive indices vary based
%                           upon this choice. Either a wavelength (in nm)
%                           may be provided, or one of the char vectors:
%                           {'vis','nir'}.
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
    % Parameters for a myopic (-3), left, human eye
    eye = modelEyeParameters('sphericalAmetropia',-3,'eyeLaterality','left');
%}


%% input parser
p = inputParser; p.KeepUnmatched = false; p.PartialMatching = false;

% Optional
p.addParameter('sphericalAmetropia',[],@(x)(isempty(x) || isscalar(x)));
p.addParameter('accommodation',[],@(x)(isempty(x) || isscalar(x)));
p.addParameter('navarroD',[],@(x)(isempty(x) || isscalar(x)));
p.addParameter('axialLength',[],@(x)(isempty(x) || isscalar(x)));
p.addParameter('eyeLaterality','Right',@ischar);
p.addParameter('species','Human',@ischar);
p.addParameter('ageYears',18,@isscalar);
p.addParameter('derivedParams',[],@(x)(isstruct(x) || isempty(x)));
p.addParameter('kvals',[],@(x)(isempty(x) || isnumeric(x)));
p.addParameter('corneaAxialRadius',[],@(x)(isempty(x) || isnumeric(x)));
p.addParameter('rotationCenterScalers',[1 1],@isnumeric);
p.addParameter('primaryPosition',[0 0],@isnumeric);
p.addParameter('spectralDomain','nir',@(x)(ischar(x) || isnumeric(x)));

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
    notes = 'Spherical refractive error and axial length provided. Relative retinal curvature set by SR, absolute size by axial length';
end

% Meta data regarding properties of the model
eye.meta.p = p.Results;
eye.meta.units = 'mm';
eye.meta.coordinates = 'eyeWorld';
eye.meta.dimensions = {'depth (axial)' 'horizontal' 'vertical'};
eye.meta.eyeLaterality = eyeLaterality;
eye.meta.sphericalAmetropia = sphericalAmetropia;
eye.meta.accommodation = p.Results.accommodation;
eye.meta.navarroD = p.Results.navarroD;
eye.meta.axialLength = p.Results.axialLength;
eye.meta.species = p.Results.species;
eye.meta.ageYears = p.Results.ageYears;
eye.meta.kvals = p.Results.kvals;
eye.meta.corneaAxialRadius = p.Results.corneaAxialRadius;
eye.meta.rotationCenterScalers = p.Results.rotationCenterScalers;
eye.meta.primaryPosition = p.Results.primaryPosition;
eye.meta.spectralDomain = p.Results.spectralDomain;
eye.meta.notes = notes;
eye.meta.varargin = varargin;

% Switch parameters at the top level by species
switch eye.meta.species
    
    %% Human
    case {'human','Human','HUMAN'}
        
        % Obtain the derived params
        if isempty(p.Results.derivedParams)
            filename = fullfile(replace(mfilename('fullpath'),mfilename(),''),'species','+human','derivedParams.mat');
            load(filename,'derivedParams');
            eye.derivedParams = derivedParams;
        else
            eye.derivedParams = p.Results.derivedParams;
        end

        % Refractive indices
        eye.index.vitreous = returnRefractiveIndex( 'vitreous', p.Results.spectralDomain );
        eye.index.aqueous = returnRefractiveIndex( 'aqueous', p.Results.spectralDomain );

        % Eye anatomy
        eye.cornea = human.cornea(eye);
        eye.iris = human.iris(eye);
        eye.stop = human.stop(eye);
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
                
        % Anatomical and optical landmarks
        eye.landmarks.medialCanthus = human.landmarks.medialCanthus(eye);
        eye.landmarks.lateralCanthus = human.landmarks.lateralCanthus(eye);
        eye.landmarks.vertex = human.landmarks.vertex(eye);
        [eye.landmarks.incidentNode,eye.landmarks.emergentNode] = human.landmarks.nodes(eye);
        eye.landmarks.fovea = human.landmarks.fovea(eye);
        eye.landmarks.opticDisc = human.landmarks.opticDisc(eye);
        
        
    %% Drasdo & Fowler 1974 model eye
    case {'Drasdo','drasdo'}
        
        % Eye anatomy
        eye.cornea = drasdo.cornea(eye);
        eye.stop = drasdo.stop(eye);
        eye.lens = drasdo.lens(eye);
        eye.retina = drasdo.retina(eye);
        
        % Index of refraction for the vitreous and aqueous
        eye.index.vitreous = 1.336;
        eye.index.aqueous = 1.336;
        
        % Landmarks
        eye.landmarks.vertex = drasdo.landmarks.vertex(eye);
        [eye.landmarks.incidentNode,eye.landmarks.emergentNode] = drasdo.landmarks.nodes(eye);
        eye.landmarks.fovea = drasdo.landmarks.fovea(eye);
        
    %% Canine
    case {'dog','Dog','canine','Canine'}
        error('Geoff needs to implement the canine model here');
        
    otherwise
        error('Please specify a valid species for the eye model');
end



end % function



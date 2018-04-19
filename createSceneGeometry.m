function sceneGeometry = createSceneGeometry(varargin)
% Create and return a sceneGeometry structure
%
% Syntax:
%  sceneGeometry = createSceneGeometry()
%
% Description:
%   Using default values and passed key/value pairs, this routine creates a
%   sceneGeometry structure, with fields that the describe a camera, an
%   eye, corrective lenses, and the geometric relationship between them.
%   The fields are:
%
%  'cameraIntrinsic' - A structure that defines the properties of a pinhole
%       camera model. Sub-fields:
%
%      'matrix' - A 3x3 matrix of the form:
%
%               [fx  s x0; 0  fy y0; 0   0  1]
%
%           where fx and fy are the focal lengths of the camera in the x
%           and y image dimensions, s is the axis skew, and x0, y0 define
%           the principle offset point. For a camera sensor with square
%           pixels, fx = fy. Ideally, skew should be zero. The principle
%           offset point should be in the center of the sensor. Units are
%           traditionally in pixels. These values can be empirically
%           measured for a camera using a calibration approach
%           (https://www.mathworks.com/help/vision/ref/cameramatrix.html).
%           Note that the values in the matrix returned by the Matlab
%           camera calibration routine must be re-arranged to correspond to
%           the X Y image dimension specifications we use here.
%       
%      'radialDistortion' - A 1x2 vector that models the radial distortion
%           introduced the lens. This is an empirically measured
%           property of the camera system.
%
%      'sensorResolution - A 1x2 vector that provides the dimension of the
%           camera image in pixels, along the X and Y dimensions,
%           respectively.
%
%  'cameraExtrinsic' - A structure that defines the spatial position of the
%       camera w.r.t. to eye. Sub-fields:
%
%      'translation' - A 3x1 vector of the form [x; y; z], with the values
%           specifying the location (horizontal, vertical, and depth,
%           respectively) of the principle offset point of the camera in mm
%           relative to the scene coordinate system. We define the origin
%           of the scene coordinate system to be x=0, y=0 along the optical
%           axis of the eye, and z=0 to be the apex of the corneal surface.
%
%      'rotationZ' - Scalar in units of degrees. This is used to create the
%           camera rotation matrix. The rotation matrix specifies the
%           rotation of the camera relative to the axes of the world
%           coordinate system. Because eye rotations are set have a value
%           of zero when the camera axis is aligned with the pupil axis of
%           the eye, the camera rotation around the X and Y axes of the
%           coordinate system are fixed at zero. Rotation about the Z axis
%           is meaningful, as the model eye has different rotation
%           properties in the azimuthal and elevational directions, and has
%           a non circular exit pupil. We specify the camera rotation
%           around the Z axis in degrees. In pupilProjection_fwd, this is
%           used to construct the rotation matrix.
%
%      'primaryPosition' - A 1x2 vector of [eyeAzimuth, eyeElevation] at
%           which the eye is in primary position (as defined by Listing's
%           Law) and thus has zero torsion.
%
%  'eye' - A structure that is returned by the function modelEyeParameters.
%       The parameters define the anatomical properties of the eye. These
%       parameters are adjusted for the measured spherical refractive error
%       of the subject and (optionally) measured axial length. Unmatched
%       key-value pairs passed to createSceneGeometry are passed to
%       modelEyeParameters.
%
%  'refraction' - A structure that identifies the function to be used
%       to compute the virtual image location of eyeWorldPoints subject to
%       refraction by the optical system (cornea and corrective lenses, if
%       any). The MATLAB virtualImageFunc is specified if the compiled MEX
%       version is not available. Sub-fields:
%           
%          'handle'       - Handle for the function.
%          'path'         - Full path to the function.
%          'opticalSystem' - A structure with the sub-fields 'p1p2' and
%                           'p1p3', each of which is an 10x4 matrix. The
%                           first m rows contain information regarding the
%                           m surfaces of the cornea, and any corrective
%                           lenses, into a format needed for ray tracing.
%                           The trailing (10-m) rows contain nans. The
%                           optical system must be a fixed size matrix so
%                           that the compiled virtualImageFuncMex is able
%                           to pre-allocate variables. The nan rows are
%                           later stripped by the rayTraceCenteredSurfaces
%                           function. The p1p2 and p1p3 sub-fields
%                           correspond to the two planes in eyeWorld
%                           coordinates. This allows the model to account
%                           for optical surfaces having a different
%                           elliptical cross section in the axial and
%                           sagittal dimension.
%
%  'lenses' - An optional structure that describes the properties of 
%       refractive lenses that are present between the eye and the camera.
%       Possible sub-fields are 'contact' and 'spectacle', each of which
%       holds the parameters that were used to add a refractive lens to the
%       optical path.
%
%   constraintTolerance - A scalar. This value is used by the function 
%       pupilProjection_inv. The inverse projection from an ellipse on the
%       image plane to eye params (azimuth, elevation) imposes a constraint
%       on how well the solution must match the shape and area of the
%       ellipse. This constraint is expressed as a proportion of error,
%       relative to either the largest possible area in ellipse shape or an
%       error in area equal to the area of the ellipse itself (i.e.,
%       unity). If the constraint is made too stringent, then in an effort
%       to perfectly match the shape of the ellipse, error will increase in
%       matching the position of the ellipse center. It should be noted
%       that even in a noise-free simulation, it is not possible to
%       perfectly match ellipse shape and area while matching ellipse
%       center position, as the shape of the projection of the pupil upon
%       the image plane deviates from perfectly elliptical due to
%       perspective effects. We find that a value in the range 0.01 - 0.03
%       provides an acceptable compromise in empirical data.
%
%   meta - A structure that contains information regarding the creation and
%       modification of the sceneGeometry.
%
%
% Inputs:
%   none
%
% Optional key/value pairs
%  'sceneGeometryFileName' - Full path to file 
%  'intrinsicCameraMatrix' - 3x3 matrix
%  'radialDistortionVector' - 1x2 vector of radial distortion parameters
%  'cameraRotationZ'      - Scalar.
%  'extrinsicTranslationVector' - 3x1 vector
%  'primaryPosition'      - 1x3 vector
%  'constraintTolerance'  - Scalar. Range 0-1. Typical value 0.01 - 0.03
%  'contactLens'          - Scalar or 1x2 vector, with values for the lens
%                           refraction in diopters, and (optionally) the
%                           index of refraction of the lens material. If
%                           left empty, no contact lens is added to the
%                           model.
%  'spectacleLens'        - Scalar, 1x2, or 1x3 vector, with values for the
%                           lens refraction in diopters, (optionally) the
%                           index of refraction of the lens material, and
%                           (optinally) the vertex distance in mm. If left
%                           empty, no spectacle is added to the model.
%  'medium'               - String, options include:
%                           {'air','water','vacuum'}. This sets the index
%                           of refraction of the medium between the eye an
%                           the camera.
%  'spectralDomain'       - String, options include {'vis','nir'}.
%                           This is the light domain within which imaging
%                           is being performed. The refractive indices vary
%                           based upon this choice.
%  'forceMATLABVirtualImageFunc' - Logical, default false. If set to
%                           true, the native MATLAB code for the
%                           virtualImageFunc is used for refraction,
%                           instead of a compiled MEX file. This is used
%                           for debugging and demonstration purposes.
%
% Outputs
%	sceneGeometry         - A structure.
%
% Examples:
%{
    % Create a scene geometry file for a myopic eye wearing a contact lens.
    % The key-value sphericalAmetropia is passed to modelEyeParameters
    % within the routine
    sceneGeometry = createSceneGeometry('sphericalAmetropia',-2,'contactLens',-2);
%}
%{
    % Create a scene geometry file for a hyperopic eye wearing spectacles
    % that provide appropriate correction when underwater.
    sceneGeometry = createSceneGeometry('sphericalAmetropia',-2,'spectacleLens',2,'medium','water');
    % Plot a figure that traces a ray arising from the optical axis at the
    % pupil plane, departing at 15 degrees.    
    figureFlag.zLim = [-15 20]; figureFlag.hLim = [-10 10];
    rayTraceCenteredSurfaces([-3.7 2], deg2rad(15), sceneGeometry.refraction.opticalSystem.p1p2,figureFlag);
%}


%% input parser
p = inputParser; p.KeepUnmatched = true;

% Optional analysis params
p.addParameter('sceneGeometryFileName','', @(x)(isempty(x) | ischar(x)));
p.addParameter('intrinsicCameraMatrix',[2600 0 320; 0 2600 240; 0 0 1],@isnumeric);
p.addParameter('sensorResolution',[640 480],@isnumeric);
p.addParameter('radialDistortionVector',[0 0],@isnumeric);
p.addParameter('extrinsicTranslationVector',[0; 0; 120],@isnumeric);
p.addParameter('cameraRotationZ',0,@isnumeric);
p.addParameter('constraintTolerance',0.02,@isscalar);
p.addParameter('contactLens',[], @(x)(isempty(x) | isnumeric(x)));
p.addParameter('spectacleLens',[], @(x)(isempty(x) | isnumeric(x)));
p.addParameter('medium','air',@ischar);
p.addParameter('spectralDomain','nir',@ischar);
p.addParameter('forceMATLABVirtualImageFunc',false,@islogical);

% parse
p.parse(varargin{:})


%% cameraIntrinsic
sceneGeometry.cameraIntrinsic.matrix = p.Results.intrinsicCameraMatrix;
sceneGeometry.cameraIntrinsic.radialDistortion = p.Results.radialDistortionVector;
sceneGeometry.cameraIntrinsic.sensorResolution = p.Results.sensorResolution;

%% cameraExtrinsic
sceneGeometry.cameraExtrinsic.translation = p.Results.extrinsicTranslationVector;
sceneGeometry.cameraExtrinsic.rotationZ = p.Results.cameraRotationZ;
sceneGeometry.cameraExtrinsic.primaryPosition = [0,0];

%% eye
sceneGeometry.eye = modelEyeParameters('spectralDomain',p.Results.spectralDomain,varargin{:});

%% refraction - handle and path
% Handle to the function; use the MEX version if available
if exist('virtualImageFuncMex')==3 && ~p.Results.forceMATLABVirtualImageFunc
    sceneGeometry.refraction.handle = @virtualImageFuncMex;
    sceneGeometry.refraction.path = which('virtualImageFuncMex');
else
    sceneGeometry.refraction.handle = @virtualImageFunc;
    sceneGeometry.refraction.path = which('virtualImageFunc');
end

%% refraction - optical system

% Assemble the opticalSystem. First get the refractive index of the medium
% between the eye and the camera
mediumRefractiveIndex = returnRefractiveIndex( p.Results.medium, p.Results.spectralDomain );

% The center of the cornea front surface is at a position equal to its
% radius of curvature, thus placing the apex of the front corneal surface
% at a z position of zero. The back surface is shifted back to produce
% the appropriate corneal thickness.
cornealThickness = -sceneGeometry.eye.cornea.back.center(1)-sceneGeometry.eye.cornea.back.radii(1);

% The axis of the cornea is rotated w.r.t. the optical axis of the eye.
% Here, we derive the radii for the ellipse that is the intersection of the
% p1p2 and p1p3 planes with the ellipsoids for the back and front corneal
% surfaces
corneaBackRotRadii=ellipsesFromEllipsoid(sceneGeometry.eye.cornea.back.radii,sceneGeometry.eye.cornea.axis);
corneaFrontRotRadii=ellipsesFromEllipsoid(sceneGeometry.eye.cornea.front.radii,sceneGeometry.eye.cornea.axis);

% Build the optical system matrix for the p1p2 and p1p3 planes. We require
% both as the cornea is not radially symmetric. The p1p2 system is for the
% horizontal (axial) plane of the eye, and the p1p3 system for the vertical
% (sagittal) plane of the eye.
sceneGeometry.refraction.opticalSystem.p1p2 = [nan, nan, nan, sceneGeometry.eye.index.aqueous; ...
    -sceneGeometry.eye.cornea.back.radii(1)-cornealThickness, -corneaBackRotRadii(1), -corneaBackRotRadii(2),  sceneGeometry.eye.index.cornea; ...
    -sceneGeometry.eye.cornea.front.radii(1), -corneaFrontRotRadii(1), -corneaFrontRotRadii(2), mediumRefractiveIndex];
sceneGeometry.refraction.opticalSystem.p1p3 = [nan, nan, nan, sceneGeometry.eye.index.aqueous; ...
    -sceneGeometry.eye.cornea.back.radii(1)-cornealThickness, -corneaBackRotRadii(1), -corneaBackRotRadii(3),  sceneGeometry.eye.index.cornea; ...
    -sceneGeometry.eye.cornea.front.radii(1), -corneaFrontRotRadii(1), -corneaFrontRotRadii(3), mediumRefractiveIndex];

%% Lenses
% Add a contact lens if requested
if ~isempty(p.Results.contactLens)
    switch length(p.Results.contactLens)
        case 1
            lensRefractiveIndex=returnRefractiveIndex( 'hydrogel', p.Results.spectralDomain );
            [sceneGeometry.refraction.opticalSystem.p1p2, ~] = addContactLens(sceneGeometry.refraction.opticalSystem.p1p2, p.Results.contactLens, 'lensRefractiveIndex', lensRefractiveIndex);
            [sceneGeometry.refraction.opticalSystem.p1p3, pOutFun] = addContactLens(sceneGeometry.refraction.opticalSystem.p1p3, p.Results.contactLens, 'lensRefractiveIndex', lensRefractiveIndex);
        case 2
            [sceneGeometry.refraction.opticalSystem.p1p2, ~] = addContactLens(sceneGeometry.refraction.opticalSystem.p1p2, p.Results.contactLens(1), 'lensRefractiveIndex', p.Results.contactLens(2));
            [sceneGeometry.refraction.opticalSystem.p1p3, pOutFun] = addContactLens(sceneGeometry.refraction.opticalSystem.p1p3, p.Results.contactLens(1), 'lensRefractiveIndex', p.Results.contactLens(2));
        otherwise
            error('The key-value pair contactLens is limited to two elements: [refractionDiopters, refractionIndex]');
    end
    sceneGeometry.lenses.contact = pOutFun.Results;
end

% Add a spectacle lens if requested
if ~isempty(p.Results.spectacleLens)
    switch length(p.Results.spectacleLens)
        case 1
            lensRefractiveIndex=returnRefractiveIndex( 'polycarbonate', p.Results.spectralDomain );
            [sceneGeometry.refraction.opticalSystem.p1p2, ~] = addSpectacleLens(sceneGeometry.refraction.opticalSystem.p1p2, p.Results.spectacleLens, 'lensRefractiveIndex', lensRefractiveIndex);
            [sceneGeometry.refraction.opticalSystem.p1p3, pOutFun] = addSpectacleLens(sceneGeometry.refraction.opticalSystem.p1p3, p.Results.spectacleLens, 'lensRefractiveIndex', lensRefractiveIndex);
        case 2
            [sceneGeometry.refraction.opticalSystem.p1p2, ~] = addSpectacleLens(sceneGeometry.refraction.opticalSystem.p1p2, p.Results.spectacleLens, 'lensRefractiveIndex', p.Results.spectacleLens(2));
            [sceneGeometry.refraction.opticalSystem.p1p3, pOutFun] = addSpectacleLens(sceneGeometry.refraction.opticalSystem.p1p3, p.Results.spectacleLens, 'lensRefractiveIndex', p.Results.spectacleLens(2));
        case 3
            [sceneGeometry.refraction.opticalSystem.p1p2, ~] = addSpectacleLens(sceneGeometry.refraction.opticalSystem.p1p2, p.Results.spectacleLens, 'lensRefractiveIndex', p.Results.spectacleLens(2),'lensVertexDistance', p.Results.spectacleLens(3));
            [sceneGeometry.refraction.opticalSystem.p1p3, pOutFun] = addSpectacleLens(sceneGeometry.refraction.opticalSystem.p1p3, p.Results.spectacleLens, 'lensRefractiveIndex', p.Results.spectacleLens(2),'lensVertexDistance', p.Results.spectacleLens(3));
        otherwise
            error('The key-value pair spectacleLens is limited to three elements: [refractionDiopters, refractionIndex, vertexDistance]');
    end
    sceneGeometry.lenses.spectacle = pOutFun.Results;
end

% Pad the optical system with nan rows to reach a fixed 10x4 size
sceneGeometry.refraction.opticalSystem.p1p2 = [sceneGeometry.refraction.opticalSystem.p1p2; ...
    nan(10-size(sceneGeometry.refraction.opticalSystem.p1p2,1),4)];
sceneGeometry.refraction.opticalSystem.p1p3 = [sceneGeometry.refraction.opticalSystem.p1p3; ...
    nan(10-size(sceneGeometry.refraction.opticalSystem.p1p3,1),4)];


%% constraintTolerance
sceneGeometry.constraintTolerance = p.Results.constraintTolerance;

%% meta
sceneGeometry.meta.createSceneGeometry = p.Results;


%% Save the sceneGeometry file
if ~isempty(p.Results.sceneGeometryFileName)
    save(p.Results.sceneGeometryFileName,'sceneGeometry');
end

end % createSceneGeometry


%% LOCAL FUNCTIONS

function rotRadii = ellipsesFromEllipsoid(radii,angles)
% Returns ellipse radii that are derived from a rotated ellipsoid
%
% Syntax:
%  rotRadii = ellipsesFromEllipsoid(radii,angles)
%
% Description:
%   The ellipsoids that describe the back and front surface of the cornea
%   are rotated with respect to the optical axis of the eye. Here, we
%   calculate the ellipse radii that correspond to the p1p2 and p1p3 axes
%   as they intersect with the rotated ellipsoid.
%


% The angles specify the rotation of the corneal ellipsoid w.r.t. the
% optical axis of the eye. As we are rotating the axes here, we need to
% take the negative of the angles.
angles = -angles;

R.azi = [cosd(angles(1)) -sind(angles(1)) 0; sind(angles(1)) cosd(angles(1)) 0; 0 0 1];
R.ele = [cosd(angles(2)) 0 sind(angles(2)); 0 1 0; -sind(angles(2)) 0 cosd(angles(2))];
R.tor = [1 0 0; 0 cosd(angles(3)) -sind(angles(3)); 0 sind(angles(3)) cosd(angles(3))];

rotMat = R.tor * R.ele * R.azi;

% Obtain the semi-axes of the ellipses in each of the planes p1p2 and p1p3

% p1p2
rotPlane = rotMat * [0; 1; 0];
[Aye,Bye]=EllipsoidPlaneIntersection(rotPlane(1),rotPlane(2),rotPlane(3),0,radii(1),radii(2),radii(3));
rotRadii(1:2) = [Aye,Bye];

% p1p3
rotPlane = rotMat * [0; 0; 1];
[~,Bye]=EllipsoidPlaneIntersection(rotPlane(1),rotPlane(2),rotPlane(3),0,radii(1),radii(2),radii(3));

rotRadii(3) = Bye;

end

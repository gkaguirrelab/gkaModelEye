function [M, stdM] = calcAngularMagnification(eye, varargin)
% Calcuates the percept angular magnification produced by artificial lenses
%
% Syntax:
%  magnification = calcAngularMagnification(eye)
%
% Description
%   Calculates the angular magnification (as seen by the retina) produced
%   by the addition of artificial lenses (contacts or spectacles). Values
%   greater than 1 indicate magnification, values less than 1 reflect
%   minification.
%
% Inputs:
%   eye                   - Structure. SEE: modelEyeParameters
%
% Optional key/value pairs:
%  'cameraMedium'         - String, options include:
%                           {'air','water','vacuum'}. This sets the index
%                           of refraction of the medium between the eye and
%                           the camera.
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
%  'targetDistance'       - Scalar. Distance in mm of the point that is the
%                           source of the magnified image.
%
% Outputs:
%   M                     - Scalar. The central tendency of the
%                           magnification produced by the system.
%   stdM                  - Scalar. The standard deviation of the
%                           magnification as measured from a range of
%                           positions (±20 degrees) in the visual field.
%
% Examples:
%{
    eye = modelEyeParameters;
    diopters = -4;
    [M, stdM] = calcAngularMagnification(eye,'spectacleLens',diopters);
    outline = sprintf('The angular magnfication on the retina produced by a %d diopter spectacle lens is x%2.2f \n',diopters,M);
    fprintf(outline)
    % Compare the value to the value in the table reported here:
    %   https://www.opticianonline.net/cet-archive/41
    assert(abs(M-0.942)<0.01)
%}
%{
    % Replicate Figure 1 of: WESTHEIMER, GERALD. "The visual world of the
    % new contact lens wearer." The Australian Journal of Optometry 46.5
    % (1963): 124-127.
    % Using the same kvals as Westheimer
    kvals= [43.5, 43.5, 0, 0, 0];
    eye = modelEyeParameters('kvals',kvals);
    d = -11:1:7;
    for ii = 1:length(d)
        Mspectacle(ii) = calcAngularMagnification(eye,'spectacleLens',d(ii));
        Mcontact(ii) = calcAngularMagnification(eye,'contactLens',d(ii));
    end
    figure
    pHandle(1) = plot(d,100*(Mspectacle-1),'-r');
    hold on
    pHandle(2) = plot(d,100*(Mcontact-1),'-b');
    xlim([-12.5 8]);
    ylim([-17.5 20]);
    plot([0 0],[-17.5 20],':k')
    plot([-12.5 8],[0 0],':k')
    ylabel('Magnification [%]');
    xlabel('Spherical correction [diopters]');
    legend( pHandle,{'spectacle','contact'},'Location','northwest');
    title('Replicate Westheimer 1963');
%}


%% input parser
p = inputParser; p.KeepUnmatched = true;

p.addRequired('eye',@isstruct);

% Optional analysis params
p.addParameter('cameraMedium','air',@ischar);
p.addParameter('contactLens',[], @(x)(isempty(x) | isnumeric(x)));
p.addParameter('spectacleLens',[], @(x)(isempty(x) | isnumeric(x)));
p.addParameter('targetDistance',1500,@isnumeric);

% parse
p.parse(eye, varargin{:})


%% Create the initial optical systems
opticalSystemRotInitial = assembleOpticalSystem(eye, 'surfaceSetName', 'stopToMedium', 'cameraMedium', p.Results.cameraMedium, 'opticalSystemNumRows', []);
opticalSystemFixInitial = initializeOpticalSystem(returnRefractiveIndex( p.Results.cameraMedium, eye.meta.spectralDomain));
opticalSystemRot = opticalSystemRotInitial;
opticalSystemFix = opticalSystemFixInitial;


%% Add the contact lens
if ~isempty(p.Results.contactLens)
    
    switch length(p.Results.contactLens)
        case 1
            lensRefractiveIndex=returnRefractiveIndex( 'hydrogel', eye.meta.spectralDomain );
            opticalSystemRot = addContactLens(opticalSystemRotInitial, p.Results.contactLens, 'lensRefractiveIndex', lensRefractiveIndex,'cornealRotation',eye.cornea.rotation);
        case 2
            opticalSystemRot = addContactLens(opticalSystemRotInitial, p.Results.contactLens(1), 'lensRefractiveIndex', p.Results.contactLens(2),'cornealRotation',eye.cornea.rotation);
        otherwise
            error('The key-value pair contactLens is limited to two elements: [refractionDiopters, refractionIndex]');
    end
    
end


%% Add the spectacle lens
if ~isempty(p.Results.spectacleLens)
    
    switch length(p.Results.spectacleLens)
        case 1
            lensRefractiveIndex=returnRefractiveIndex( 'polycarbonate', eye.meta.spectralDomain );
            opticalSystemFix = addSpectacleLens(opticalSystemFixInitial, p.Results.spectacleLens, 'lensRefractiveIndex', lensRefractiveIndex);
        case 2
            opticalSystemFix = addSpectacleLens(opticalSystemFixInitial, p.Results.spectacleLens(1), 'lensRefractiveIndex', p.Results.spectacleLens(2));
        case 3
            opticalSystemFix = addSpectacleLens(opticalSystemFixInitial, p.Results.spectacleLens(1), 'lensRefractiveIndex', p.Results.spectacleLens(2),'lensVertexDistance', p.Results.spectacleLens(3));
        case 4
            opticalSystemFix = addSpectacleLens(opticalSystemFixInitial, p.Results.spectacleLens(1), 'lensRefractiveIndex', p.Results.spectacleLens(2),'lensVertexDistance', p.Results.spectacleLens(3), 'baseCurve', p.Results.spectacleLens(4));
        otherwise
            error('The key-value pair spectacleLens is limited to four elements: [refractionDiopters, refractionIndex, vertexDistance, baseCurve]');
    end
    
end


%% Define some arguments for findPupilRay
eyeCoordOrigin = eye.stop.center;
eyePose = [0 0 0 0];
rotationCenters = eye.rotationCenters;


%% Loop through world targets
% The calculated magnification varies a bit depending upon the point from
% which it is calculated, so perform the calculation from a few points

horiz = [-p.Results.targetDistance/5 p.Results.targetDistance/5];
vert = [-p.Results.targetDistance/5 p.Results.targetDistance/5];
M = [];

for hh=1:length(horiz)
    for vv = 1:length(vert)
        
        %% Assemble the world target
        worldTarget = [horiz(hh);vert(vv);p.Results.targetDistance];
        
        
        %% The angle with which a ray departs the entrance pupil
        [~, initialRay]  = findPupilRay( eyeCoordOrigin, eyePose, worldTarget, rotationCenters, opticalSystemRotInitial, opticalSystemFixInitial );
        [p1p2A,p1p3A] = quadric.rayToAngles(initialRay);
        
        
        %% Repeat the ray trace with the lens correction
        [~, initialRay]  = findPupilRay( eyeCoordOrigin, eyePose, worldTarget, rotationCenters, opticalSystemRot, opticalSystemFix );
        [p1p2B,p1p3B] = quadric.rayToAngles(initialRay);
        
        
        %% Calculate magnification
        % The magnification is given by the ratio of the corrected and uncorrected
        % angles
        M(end+1) = mean([p1p2B/p1p2A, p1p3B/p1p3A]);
        
    end
end

% Report the std and mean magnification across the sampled positions
stdM = std(M);
M = mean(M);

end


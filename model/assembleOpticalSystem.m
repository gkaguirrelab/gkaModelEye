function [opticalSystem, surfaceLabels, surfaceColors] = assembleOpticalSystem( eye, varargin )


% Assemble the system from retina to camera

%% input parser
p = inputParser; p.KeepUnmatched = true;

p.addRequired('eye',@isstruct);

% Optional analysis params
p.addParameter('surfaceSetName','pupilToCamera',@ischar);
p.addParameter('cameraMedium','air',@ischar);
p.addParameter('contactLens',[], @(x)(isempty(x) | isnumeric(x)));
p.addParameter('spectacleLens',[], @(x)(isempty(x) | isnumeric(x)));
p.addParameter('opticalSystemNumRows',100,@isnumeric);

% parse
p.parse(eye, varargin{:})

% Get the refractive index of the medium in which the camera resides
mediumRefractiveIndex = returnRefractiveIndex( p.Results.cameraMedium, eye.meta.spectralDomain );

% Initialize an empty optical system matrix
opticalSystem=[];

switch p.Results.surfaceSetName
    case 'retinaToCamera'
        
        % Start in the retina
        opticalSystem = [opticalSystem; ...
            [eye.posteriorChamber.S eye.posteriorChamber.side eye.posteriorChamber.boundingBox eye.posteriorChamber.mustIntersect eye.index.vitreous]];
        
        % Add the lens
        opticalSystem = [opticalSystem; ...
            [eye.lens.S eye.lens.side eye.lens.boundingBox eye.lens.mustIntersect [eye.lens.index; eye.index.aqueous]]];
        
        % Add the cornea
        opticalSystem = [opticalSystem; ...
            [eye.cornea.S eye.cornea.side eye.cornea.boundingBox eye.cornea.mustIntersect [eye.cornea.index; mediumRefractiveIndex]]];
        
        % Assemble the labels
        surfaceLabels = [eye.posteriorChamber.label; eye.lens.label; eye.cornea.label];
        
        % Assemble the surface plot colors
        surfaceColors = [eye.posteriorChamber.plot.color; eye.lens.plot.color; eye.cornea.plot.color];
        
        % % Add a contact lens if requested
        % if ~isempty(p.Results.contactLens)
        %     switch length(p.Results.contactLens)
        %         case 1
        %             lensRefractiveIndex=returnRefractiveIndex( 'hydrogel', p.Results.spectralDomain );
        %             [opticalSystem, pOutFun] = addContactLens(opticalSystem, p.Results.contactLens, 'lensRefractiveIndex', lensRefractiveIndex);
        %         case 2
        %             [opticalSystem, pOutFun] = addContactLens(opticalSystem, p.Results.contactLens(1), 'lensRefractiveIndex', p.Results.contactLens(2));
        %         otherwise
        %             error('The key-value pair contactLens is limited to two elements: [refractionDiopters, refractionIndex]');
        %     end
        %     sceneGeometry.lenses.contact = pOutFun.Results;
        % end
        %
        % % Add a spectacle lens if requested
        % if ~isempty(p.Results.spectacleLens)
        %     switch length(p.Results.spectacleLens)
        %         case 1
        %             lensRefractiveIndex=returnRefractiveIndex( 'polycarbonate', p.Results.spectralDomain );
        %             [opticalSystem, pOutFun] = addSpectacleLens(opticalSystem, p.Results.spectacleLens, 'lensRefractiveIndex', lensRefractiveIndex);
        %         case 2
        %             [opticalSystem, pOutFun] = addSpectacleLens(opticalSystem, p.Results.spectacleLens, 'lensRefractiveIndex', p.Results.spectacleLens(2));
        %         case 3
        %             [opticalSystem, pOutFun] = addSpectacleLens(opticalSystem, p.Results.spectacleLens, 'lensRefractiveIndex', p.Results.spectacleLens(2),'lensVertexDistance', p.Results.spectacleLens(3));
        %         otherwise
        %             error('The key-value pair spectacleLens is limited to three elements: [refractionDiopters, refractionIndex, vertexDistance]');
        %     end
        %     sceneGeometry.lenses.spectacle = pOutFun.Results;
        % end
        
    case 'pupilToCamera'
        
        % Start in the aqueous
        opticalSystem = [opticalSystem; ...
            [nan(1,10) nan nan(1,6) nan eye.index.aqueous]];
        
        % Add the cornea
        opticalSystem = [opticalSystem; ...
            [eye.cornea.S eye.cornea.side eye.cornea.boundingBox eye.cornea.mustIntersect [eye.cornea.index; mediumRefractiveIndex]]];
        
        % Assemble the labels
        surfaceLabels = [{'aqueous'}; eye.cornea.label];
        
        % Assemble the surface plot colors
        surfaceColors = [{[nan nan nan]}; eye.cornea.plot.color];
        
    otherwise
        error('Unrecognized surfaceSetName');
        
end % switch statement


%% Pad the optical system matrix
% The number of rows in the optical system matrix is set to a fixed value
% so that the compiled ray-tracing routines can expect a constant size for
% the input variables. The nan rows are stripped out at the time of ray
% tracing.
osRowLength = size(opticalSystem,2);
opticalSystem = [opticalSystem; ...
    nan(p.Results.opticalSystemNumRows-size(opticalSystem,1),osRowLength)];

end % assembleOpticalSystem

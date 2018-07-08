function opticalSystem = assembleOpticalSystem( eye, varargin )


% Assemble the system from retina to camera


%% input parser
p = inputParser; p.KeepUnmatched = true;

p.addRequired('eye',@isstruct);

% Optional analysis params
p.addParameter('contactLens',[], @(x)(isempty(x) | isnumeric(x)));
p.addParameter('spectacleLens',[], @(x)(isempty(x) | isnumeric(x)));
p.addParameter('cameraMedium','air',@ischar);
p.addParameter('spectralDomain','nir',@ischar);

% parse
p.parse(eye, varargin{:})

mediumRefractiveIndex = returnRefractiveIndex( p.Results.cameraMedium, p.Results.spectralDomain );

%% Build the optical system matrix

% Start in the retina
opticalSystem(1,:)     = [nan(1,10) nan nan(1,6) nan eye.index.aqueous];

% Add the lens
opticalSystem = [opticalSystem; ...
    [eye.lens.S eye.lens.side eye.lens.boundingBox eye.lens.mustIntersect [eye.lens.index; eye.index.vitreous]]];

% Add the cornea
opticalSystem = [opticalSystem; ...
    [eye.cornea.S eye.cornea.side eye.cornea.boundingBox eye.cornea.mustIntersect [eye.cornea.index; mediumRefractiveIndex]]];

% Add a contact lens if requested
if ~isempty(p.Results.contactLens)
    switch length(p.Results.contactLens)
        case 1
            lensRefractiveIndex=returnRefractiveIndex( 'hydrogel', p.Results.spectralDomain );
            [opticalSystem, pOutFun] = addContactLens(opticalSystem, p.Results.contactLens, 'lensRefractiveIndex', lensRefractiveIndex);
        case 2
            [opticalSystem, pOutFun] = addContactLens(opticalSystem, p.Results.contactLens(1), 'lensRefractiveIndex', p.Results.contactLens(2));
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
            [opticalSystem, pOutFun] = addSpectacleLens(opticalSystem, p.Results.spectacleLens, 'lensRefractiveIndex', lensRefractiveIndex);
        case 2
            [opticalSystem, pOutFun] = addSpectacleLens(opticalSystem, p.Results.spectacleLens, 'lensRefractiveIndex', p.Results.spectacleLens(2));
        case 3
            [opticalSystem, pOutFun] = addSpectacleLens(opticalSystem, p.Results.spectacleLens, 'lensRefractiveIndex', p.Results.spectacleLens(2),'lensVertexDistance', p.Results.spectacleLens(3));
        otherwise
            error('The key-value pair spectacleLens is limited to three elements: [refractionDiopters, refractionIndex, vertexDistance]');
    end
    sceneGeometry.lenses.spectacle = pOutFun.Results;
end

% Pad the optical system with nan rows to reach a fixed 20x19 size
opticalSystem = [opticalSystem; ...
    nan(20-size(opticalSystem,1),19)];


end


function [azimuth, elevation] = retinaCoordToVisualFieldAngles( eye, eyePoint )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    eye = modelEyeParameters();

        opticalSystem = assembleOpticalSystem( eye, 'surfaceSetName', 'retinaToPupil' );

        args = {eye.pupil.center([2 3 1])', ...
                eye.rotationCenters, ...
                opticalSystem};

        eyePose = [0 0 0 2];

        vitreousChamberApex = quadric.ellipsoidalGeoToCart( [-90 -90 0], opticalSystem(1,1:10) );
        eyePoint = quadric.ellipsoidalGeoToCart( [-90 -65 0], opticalSystem(1,1:10) );

        [distance,startAngle,endAngle,~] = quadric.panouGeodesicDistance(opticalSystem(1,1:10),vitreousChamberApex,eyePoint);
        
        
            [~, ~, angle_p1p2, angle_p1p3] = ...
                virtualImageFunc(...
                eyePoint', eyePose, args{:});


end


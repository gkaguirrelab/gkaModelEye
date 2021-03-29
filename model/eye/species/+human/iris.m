function iris = iris( eye )
% Returns the iris sub-field of an eye model structure
%
% Syntax:
%  iris = human.iris( eye )
%
% Description:
%   Returns parameters that define the properties of the iris.
%
%   The iris can be modeled has having a thickness. This thickness
%   influences the properties of the entrance pupil, as when the eye is
%   rotated w.r.t. the camera either the front or back surface of the iris
%   aperture stop defines the near or far edge of the entrance pupil.
%
%   I position the anterior surface of the iris at a depth of 3.9 mm,
%   which reflects a cycloplegic eye. I model the eye with zero iris angle,
%   thus making the iris a plane. The iris is centered on the optical axis.
%
% Inputs:
%   eye                   - Structure.
%
% Outputs:
%   retina                - Structure.
%


% The model supports specification of a non-zero iris thickness. This was
% not found to improve the estimation of the appearance of the entrance
% pupil, however, so this is set to zero.
iris.thickness = 0.0;

switch eye.meta.eyeLaterality
    case 'Right'
        iris.center = [-3.9 0 0];
    case 'Left'
        iris.center = [-3.9 0 0];
end

% Define the iris radius. One study measured the horizontal visible iris
% diameter (HVID) in 200 people, and found a mean of 11.8 with a range of
% 10.2 - 13.0.
%
%    PJ Caroline & MP Andrew. "The Effect of Corneal Diameter on Soft Lens
%    Fitting, Part 1" Contact Lens Spectrum, Issue: April 2002
%    https://www.clspectrum.com/issues/2002/april-2002/contact-lens-case-reports
%
% Bernd Bruckner of the company Appenzeller Kontaktlinsen AG supplied me
% with a tech report from his company (HVID & SimK study) that measured
% HVID in 461 people. These data yield a mean iris radius of 5.92 mm, 0.28
% SD. The values from the histogram are represented here, along with a
% Gaussian fit to the distribution
%{
    counts = [0 2 2 0 0 4 5 12 19 23 36 44 52 41 39 37 43 30 28 12 15 10 4 1 0 2 0];
    HVIDRadiusmm = (10.5:0.1:13.1)/2;
    hvidGaussFit = fit(HVIDRadiusmm', counts', 'gauss1');
    hvidRadiusMean = hvidGaussFit.b1;
    hvidRadiusSD =  hvidGaussFit.c1;
    figure
    plot(HVIDRadiusmm, hvidGaussFit(HVIDRadiusmm), '-r')
    hold on
    plot(HVIDRadiusmm, counts, '*k')
    xlabel('HVID radius in mm')
    ylabel('counts')
%}
% The HVID is the refracted iris size. We can use the forward model to find
% the size of the true iris.
%{
    % Mean and SD radius value from prior block of code
    hvidRadiusMean = 5.9161;
    hvidRadiusSD = 0.2830;
    sceneGeometry = createSceneGeometry();
    % Inactivate ray-tracing
    sceneGeometry.refraction = [];
    % Get the area in pixels of a "pupil" that is the same radius
    % as the HVID when there is no ray tracing
    hvidP=projectModelEye([0 0 0 hvidRadiusMean],sceneGeometry);
    % Restore ray tracing
    sceneGeometry = createSceneGeometry();
    % Set up the objective function
    myArea = @(p) p(3);
    myObj = @(r) (hvidP(3) - myArea(projectModelEye([0 0 0 r],sceneGeometry)))^2;
    [r,pixelError] = fminsearch(myObj,5.5);
    fprintf('An unrefracted iris radius of %4.2f yields a refracted HVID of %4.2f \n',r,hvidRadiusMean)
%}
% We use this true iris size and then subject the iris perimeter points to
% refraction
iris.radius = 5.55;

end


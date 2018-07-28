function iris = iris( eye )

% The iris has a thickness. This thickness influences the properties of the
% entrance pupil, as when the eye is rotated w.r.t. the camera either the
% front or back surface of the iris aperture defines the near or far edge
% of the entrance pupil. A thickness of 0.15 mm was found to produce an
% entrance pupil ellipse that well matches the Mathur 2013 empirical
% results.
iris.thickness = 0.15;

% We position the anterior surface of the iris at a depth of 3.925 mm,
% which reflects a cycloplegic eye. We model an eye with zero iris angle,
% thus making the iris a plane. We adjust the position of the iris so that
% it is centered within the rotated corneal ellipse. This is consistent
% with reports that the iris is shifted slightly upward with respect to the
% pupil center, although inconsistent with the report that it is shifted
% temporally:
%
%   ...the typical entrance pupil is decentered approximately 0.15 mm
%   nasally and 0.1 mm inferior to the geometric center of the visible iris
%   circumference
%
% Bennett, Edward S., and Barry A. Weissman, eds. Clinical contact lens
% practice. Lippincott Williams & Wilkins, 2005, p119
switch eye.meta.eyeLaterality
    case 'Right'
        iris.center = [-4+iris.thickness/2 0.35 0.35];
    case 'Left'
        iris.center = [-4+iris.thickness/2 -0.35 0.35];
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
    sceneGeometry = createSceneGeometry();
    sceneGeometry.refraction = [];
    % Get the area in pixels of a "pupil" that is the same radius
    % as the HVID when there is no ray tracing
    hvidP=pupilProjection_fwd([0 0 0 hvidRadiusMean],sceneGeometry);
    % Restore ray tracing
    sceneGeometry = createSceneGeometry();
    % Set up the objective function
    myArea = @(p) p(3);
    myObj = @(r) (hvidP(3) - myArea(pupilProjection_fwd([0 0 0 r],sceneGeometry)))^2;
    [r,pixelError] = fminsearch(myObj,5.5);
    fprintf('An unrefracted iris radius of %4.2f yields a refracted HVID of %4.2f \n',r,hvidRadiusMean)
%}
% We use this true iris size and then subject the iris perimeter points to
% refraction
iris.radius = 5.57;

end


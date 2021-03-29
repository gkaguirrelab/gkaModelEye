
function plotConicSection(S,plane,colorCode,boundingBox)
F = quadric.vecToFunc(S);
switch plane
    case 'Horizontal'
        fh = @(x,y) F(x,y,0);
        rangeVec = boundingBox([1 2 3 4]);
    case 'Vertical'
        fh = @(x,y) F(x,0,y);
        rangeVec = boundingBox([1 2 5 6]);
    case 'Coronal'
        fh = @(x,y) F(0,x,y);
        rangeVec = boundingBox([3 4 5 6]);
end
fimplicit(fh,rangeVec,'Color', colorCode,'LineWidth',1);
end

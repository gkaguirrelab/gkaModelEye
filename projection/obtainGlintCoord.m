function glintCoord = obtainGlintCoord(imagePoints,pointLabels)

glintIdx = strncmp(pointLabels,'glint',5);
if any(glintIdx)
    glintCoord = imagePoints(glintIdx,:);
else
    glintCoord = [];
end

end


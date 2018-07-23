% This function finds the vertices of new material that is added. The
% centre where new material is added (insertion point) should be provided

% INPUTS:
% centre; 2x1 row vector specifying the location of the insertion point
% that added the new material
% numVertices; scalar, number of vertices that make up the polygon that describes
% this region
% amountMaterial; scalar, amount of material that should be inserted in one
% timestep delta-t

% OUTPUTS:
% vertices; numVerticesx2 matrix with column 1 being the x-coords and
% column 2 the y-coords of the vertices that make up the polygon that
% describes the region of newly-inserted material

function vertices = findVerticesNewMaterial(centre,numVertices,amountMaterial)

% we calculate the radius of a regular numVertices-sided polygon that would
% have an area of amountMaterial

radius = sqrt((2*amountMaterial)/(numVertices*sin(2*pi/numVertices)));

% produce a polygon and access its vertices

pgon = nsidedpoly(numVertices);

% this produces a polygon with the right number of sides but radius 1, so
% scale it so that it has the correct radius to have the right area
% (calculated above)

vertices = pgon.Vertices*radius;

% the polygon is also centred at (0,0), so shift it to the correct centre

vertices(:,1) = vertices(:,1) + centre(1);
vertices(:,2) = vertices(:,2) + centre(2);
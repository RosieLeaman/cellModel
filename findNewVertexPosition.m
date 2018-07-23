function newVertex = findNewVertexPosition(vertex,flow,dt)

newVertex = vertex + flow*dt;
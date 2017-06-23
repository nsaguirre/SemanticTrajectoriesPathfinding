Semantic Trajectories
Nathan Aguirre and Ben Meline
Work performed under Dr. Goce Trajcevski

Refer to SemanticTrajectoryWriteup for full writeup and description.

share/data contains .osm files retrieved from OpenStreetMap. Method of retrieval was URL download through bounding box: GET /api/0.6/map.
Refer to http://wiki.openstreetmap.org/wiki/API_v0.6 for more information.

share/code contains relevant code.

traj_range_query.py contains all functions used in the semantic trajectory query.
Contents are separated under headings

1. Solution - functions that perform the query
2. Time Calculation - node-to-node time calculation functions
3. File parsing - parsing input files
4. Speeds - dictionary describing speeds to associate to various types of streets
5. Graph generation - Using networkx, generate the graphs per query files
6. Node matching - Function to find the closest node
7. Points-of-Interest (POI) injection - Insertion of POI attributes and nodes into graph
8. Path searching - path searching functions, 
9. Utilities - supporting functions
10. Plotting - plotting function

SemanticTrajectoryNotebook.ipynb runs through examples of queries.

Requires networkx, shapely, xml, regex, and matplotlib
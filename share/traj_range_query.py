'''
Semantic trajectory. Compiled list of functions and classes for semantic trajectories

Nathan Aguirre
Ben Meline
Work performed under Dr. Goce Trajcevski

All relevant functions have been assembled here.

Table of Contents
1. Solution
2. Time Calculation
3. File parsing
4. Speeds
5. Graph generation
6. Node matching
7. Points-of-Interest (POI) injection
8. Path searching
9. Utilities
10. Plotting
'''

import xml.etree.ElementTree as ET
import networkx as nx
import itertools
import heapq
from math import radians, sin, cos, sqrt, asin
from shapely.geometry import Polygon, LineString, MultiLineString, Point
try:
    import regex as re
except:
    import re
import matplotlib.pyplot as plt

'''
1. Solution

'''

def run_solution(traj_file, query_file, poi_file):
    # Run full solution starting with files
    traj_dict = trajectory_parser(traj_file)
    poi_dict = poi_parser(poi_file)
    vertices, min_time, max_time, time_limit = query_parser(query_file)
        
    print('Files parsed')
    G, pos = generate_graph(vertices)
    print('Graph generated')
    G_poi, poi_node_dict, pos = inject_multiple_pois(G, poi_dict, pos, vertices)
    print('Injected POIs')

    G_poly = nodes_in_poly(G_poi)
    valid_traj = run_query(G_poi, G_poly, poi_node_dict, traj_dict, min_time, max_time, time_limit)
    
    # return G_poi, pos, valid_traj, traj_dict
    return valid_traj

def run_query(G, G_poly, poi_dict, traj_dict, min_time, max_time, time_limit):
    '''Run single query on the graphs generated from traj_file, query_file, and poi_file
    Given a trajectory file, does a path exist such that the '''
    if min_time > max_time:
        print('Maximum time is smaller than minimum time. Check values/file')
    if max_time > time_limit:
        print('Time limit is smaller than minimum time. Check values/file')
    
    valid_traj = {}
    for i in traj_dict:

        # check all possible trajectories
        is_path_valid = False
        traj_array = traj_dict[i] # coordinates along the trajectory        

        valid_path = []
        
        starting_time = 0
        steps = []
        flat_list = []
        
        # try every point in trajectory
        for a in traj_array:
            lat1, lon1 = a.lat_start, a.lon_start
            lat2, lon2 = a.lat_end, a.lon_end
            
            if a.type != 'stop':
                # not a POI, then calculate nearest node
                a_node, a_dist = closest_node(lat1,lon1, G)
                b_node, b_dist = closest_node(lat2,lon2, G)
                print(a_node,b_node)
            else:
                # if stopping at POI, then b_node is the same node
                a_node = poi_dict[a.poi]
                b_node = poi_dict[a.poi]

            steps.append((a_node, b_node))
            flat_list.extend([a_node, b_node])
        
        first_node = steps[0][0] # select first node of first pair of points
        last_node = steps[-1][1] # select last node of last pair of points
        
        # sum wait times along trajectory
        min_wait = 0
        for node in list(set(flat_list)):
            if 'wait_time' in G.node[node]:
                min_wait += G.node[node]['wait_time']
        
        # Calculate time from first node to all nodes and last node to all nodes
        global_source_times, global_target_times = make_time_dict(G, first_node, last_node, cutoff = None, weight = 'time_tag')
        # Calculate minimum time to polygon from first node and last node
        global_min_to_poly, global_min_to_target = get_min_times_through_poly_dict(G_poly, first_node, last_node, global_source_times, global_target_times)
        
        # needs to be possible to reach polygon and go to POIs on trajectory. Commented out due to difficulty in 
        # if global_min_to_target + min_wait > time_limit:
        #     print('Cannot reach polygon under time limit')
        #     continue
            
        # Make dictionary of shortest time to source and target from other nodes for each leg
        leg_dict = {}
        
        print('First steps')
        # For each leg of the trajectory, store all the shortest times fron each node to the polygon, from node to target node through polygon, and to target node
        for step in steps:
            s_node, t_node = step[0], step[1] # start and end nodes on leg
            local_source_times, local_target_times = make_time_dict(G, s_node, t_node, cutoff = None, weight = 'time_tag')
            local_min_to_poly, local_min_to_target = get_min_times_through_poly_dict(G_poly, s_node, t_node, local_source_times, local_target_times)
            to_target = local_source_times.get(t_node, float('inf'))
            leg_dict[step] = {'to_poly': local_min_to_poly, 'to_target_through_poly': local_min_to_target, 'to_target': to_target}
        print(leg_dict)

        shortest_pathing = []
        # Compare the possible order of paths. In other words, where in the trajectory can the
        # pass through the polygon occur, given the time constraints? Prints all options for shortest pathing.
        for step in steps:
            s_node, t_node = step[0], step[1]
            through_poly = leg_dict[step]['to_target_through_poly']
            leg_sum = through_poly
            for leg in leg_dict:
                if leg != step:
                    leg_sum += leg_dict[step]['to_target']
            leg_dict[step]['leg_shortest_path'] = leg_sum
            shortest_pathing.append(leg_sum)
        
        print('Shortest Pathing: ')
        print(shortest_pathing)
        if min(shortest_pathing) > time_limit:
            # impossible trajectory - move to next trajectory
            print('Impossible trajectory, time_limit too low to make valid trajectory')
            continue
            
        print('----')
        print('Passing steps')
        
        # Finally iterate over steps to find valid paths. If-else statements prevent needless path-searching
        
        for step in steps:
            
            s_node, t_node = step[0], step[1]
            
            local_source_times, local_target_times = make_time_dict(G, s_node, t_node, cutoff = None, weight = 'time_tag')
            local_min_to_poly, local_min_to_target = get_min_times_through_poly_dict(G_poly, s_node, t_node, local_source_times, local_target_times)
            print(local_min_to_poly, local_min_to_target)
            print(s_node, t_node)
            print('------------')
            leg_shortest_path = leg_dict[step]['leg_shortest_path']
            
            # Determine the minimum time required to reach polygon and minimum time required to reach target
            print((local_min_to_poly + starting_time))
            if ((local_min_to_poly + starting_time) > max_time):
                print('Cannot reach polygon under upper time bound for this leg')
                if t_node in local_source_times:
                    starting_time += local_source_times[t_node]
                else:
                    starting_time +=float('inf')
                continue
                
            print((local_min_to_target + starting_time))
            if ((local_min_to_target + starting_time) > time_limit):
                print('Cannot reach end through polygon under time limit for this leg')
                if t_node in local_source_times:
                    starting_time += local_source_times[t_node]
                else:
                    starting_time +=float('inf')
                continue
            
            if (leg_shortest_path > time_limit):
                print('Shortest path from start to finish with this leg passing through the polygon is over the time limit')
                if t_node in local_source_times:
                    starting_time += local_source_times[t_node]
                else:
                    starting_time +=float('inf')
                continue
            
            valid_path = find_path(G, s_node, t_node, starting_time, min_time, max_time, local_target_times, time_limit)

            if len(valid_path)>0:
                is_path_valid = True
                
                print('Found')
                polygon_step = step
                break
                
        if is_path_valid:
            valid_traj[i] = valid_path
            poly_node = valid_path[-1]
            full_traj = steps
            split_step = steps.index(polygon_step)
            to_poly_node, from_poly_node = steps[split_step][0], steps[split_step][1]
            full_traj[split_step] = (to_poly_node, poly_node)
            full_traj.insert(split_step+1, (poly_node, from_poly_node))
            #for j in range(len(steps)):
            #    if j == split_step:
            #        full_traj.append((steps[j][0], poly_node))
            #        full_traj.append((poly_node, steps[j][1]))
            #    else:
            #        full_traj.append(steps[i])
            print('Full trajectory:')
            print(full_traj)
            for k in range(len(full_traj)):
                if k >= split_step+1:
                    s_node, t_node = full_traj[k][0], full_traj[k][1]
                    print(s_node, t_node)
                    sub_path = nx.shortest_path(G, s_node, t_node, weight = 'time_tag')
                    valid_traj[i].extend(sub_path)

    return valid_traj

'''
2. Time Calculation
'''

def make_time_dict(G, source, target, cutoff = None, weight = 'weight'):
    # Calculates shortest path according to Dijkstra algorithm from given source node and to target node
    source_times, pair_paths_source = single_source_dijkstra_path_2(G, source, cutoff=cutoff,weight=weight)
    target_times, pair_paths_target = single_target_dijkstra_path_2(G, target, cutoff=cutoff, weight=weight)
    return source_times, target_times

def single_source_dijkstra_path_2(G,source, cutoff=None, weight='weight'):
    # Single source shortest pathfinding
    (length,path)=single_source_dijkstra_2(G,source, weight = weight)
    return length, path

def single_target_dijkstra_path_2(G, target, cutoff=None, weight='weight'):
    # Simple modification - switch the direction of all the edges and calculate the SSSP from the target
    G_rev = G.reverse(G)
    (length,path)=single_source_dijkstra_2(G_rev, target, weight = weight)
    for i in path:
        path[i] = path[i].reverse()
    return length, path

def single_source_dijkstra_2(G,source, target=None,cutoff=None,weight='weight', poi_nodes = None):
    # Modification of Networkx's single_source_dijkstra and _dijkstra functions 
    if source==target:
        return ({source:0}, {source:[source]})
    dist = {}  # dictionary of final distances
    paths = {source:[source]}  # dictionary of paths
    seen = {source:0}
    fringe=[] # use heapq with (distance,label) tuples
    heapq.heappush(fringe,(0,source))
    while fringe:
        (d,v)=heapq.heappop(fringe)
        
        is_poi = False
        if poi_nodes != None:
            if v in set(poi_nodes):
                is_poi = True
            
        if v in dist:
            continue # already searched this node.
        dist[v] = d
        if v == target:
            break
        #for ignore,w,edgedata in G.edges_iter(v,data=True):
        #is about 30% slower than the following
        if G.is_multigraph():
            edata=[]
            for w,keydata in G[v].items():
                minweight=min((dd.get(weight,1)
                               for k,dd in keydata.items()))
                edata.append((w,{weight:minweight}))
        else:
            edata=iter(G[v].items())

        for w,edgedata in edata:
            
            wait_time = 0
            if is_poi:
                wait_time = G.node[v]['wait_time']
                
            vw_dist = dist[v] + edgedata.get(weight,1) + wait_time
            if cutoff is not None:
                if vw_dist>cutoff:
                    continue
            if w in dist:
                if vw_dist < dist[w]:
                    raise ValueError('Contradictory paths found:',
                                     'negative weights?')
            elif w not in seen or vw_dist < seen[w]:
                seen[w] = vw_dist
                heapq.heappush(fringe,(vw_dist,w))
                paths[w] = paths[v]+[w]

    return (dist,paths) 
    
def get_min_times_through_poly_dict(G_poly, source, target, source_dict, target_dict):
    # Calculate shortest path from source to target through any node in the polygon region
    time_from_source = float('inf')
    time_to_target = float('inf')
    
    min_time_to_poly = float('inf')
    min_total_time = float('inf')
    
    for node in G_poly:
        if node in source_dict and node in target_dict:
            time_from_source = source_dict[node]
            time_to_target = target_dict[node]
            if time_from_source < min_time_to_poly:
                min_time_to_poly = time_from_source 
            if (time_from_source+time_to_target) < min_total_time:
                min_total_time = time_from_source+time_to_target
    
    return min_time_to_poly, min_total_time,            

'''
3. File parsing
'''

class TrajPoint:
    # TrajPoint class keeps track of results from trajectory parser.
    def __init__(self, array):
        if 'stop' in array[1]:
            self.type = 'stop'
            self.poi = array[1][5:-1]
        else:
            self.type = array[1]
            
        self.id = int(array[0])
        self.lat_start = float(array[3])
        self.lon_start = float(array[4])
        self.lat_end = float(array[6])
        self.lon_end = float(array[7])

def trajectory_parser(f):
    # Parse trajectory file. Object_Id corresponds to full trajectory. 
    fr = open(f, 'r')
    next(fr)
    traj_dict = {}

    for line in fr:
        if 'stop' in line:
            # Remove spaces between parentheses
            m = re.search(r'\([^)]*\)', line)
            m_remove = re.sub('\s','', m[0])
            line = re.sub(r'\([^)]*\)', m_remove, line)
        s = line.split()
        new_traj_point = TrajPoint(s)
        
        if new_traj_point.id not in traj_dict:
            traj_dict[new_traj_point.id] = [new_traj_point]
        else:
            array = traj_dict[new_traj_point.id]
            array.append(new_traj_point)    
    fr.close()

    return traj_dict

def poi_parser(f):
    # Parse the POI text file. Contains list of POIs and their coordinates
    poi_dict = {}
    fr = open(f, 'r')

    while True:
        line1 = fr.readline()
        if not line1:
            break
        if line1.strip() == '':
            continue
        line2 = fr.readline()
        if not line2:
            break
        
        lat, lon, wait_time = float(line2.split()[0]), float(line2.split()[1]), float(line2.split()[2])/180.0
        poi_dict[re.sub('[\s+]', '', line1)] = {'lat': lat, 'lon': lon, 'wait_time': wait_time}
    
    fr.close()
    return poi_dict

def query_parser(f):
    # Parse query file. contains polygon vertices in counterclockwise order, as well as time constraints (initial time, polygon time, and time limit)
    vertices = []
    time_range = []
    fr = open(f, 'r')
    
    for line in fr:
        if line[:2] == 'V:':
            vertex_line = line.split()
            lat = float(vertex_line[1])
            lon = float(vertex_line[2])
            vertices.append((lon,lat))

        elif line[:2] == 'T:':
            time_line = line.split()
            time_range.append(float(time_line[1]))
    
    fr.close()
    return vertices, time_range[0], time_range[1], time_range[2] # min_time, max_time, time_limit            

def parse_all(traj_file, query_file, poi_file):
    
    # Trajectory dictionary will separate each line into its associated path.
    traj_dict = trajectory_parser(traj_file)
    
    # Vertices of the polygon. Represents region the trajectory must pass through.
    vertices, start_time, end_time, time_limit = query_parser(query_file)
    
    # List of arbitrary points of interest
    poi_dict = poi_parser(poi_file)

    return traj_dict, poi_dict, vertices, start_time, end_time, time_limit

'''
4. Speeds
'''

highway_dict = {'motorway':65, 
                'trunk':45,
                'primary':45,
                'secondary':30,
                'tertiary':30,
                'unclassified':25,
                'service':20,
                'residential':25,
                'motorway_link':60,
                'trunk_link':25,
                'primary_link':40,
                'secondary_link':25,
                'tertiary_link':25}

'''
5. Graph generation
'''

def generate_graph(vertices):
    # Initialize graph as directed graph. This allows for one way streets.
    G=nx.DiGraph()

    # Open the retrieved osm files containing a region of Chicago.
    for i in range(14):
        #text = './/map.osm-' +str(i+1)+'.xml'
        text = './data/map.osm-' +str(i+1)+'.xml'
        tree = ET.parse(text)
        root = tree.getroot()
        make_graph_with_poly(G,root,vertices) # With or without poly
        print(str(i+1))

    # This set of nodes is the largest connected cluster of node. This allows for better handling when it comes to pathing, so use G_big for analysis
    big_cluster_nodes = sorted(nx.weakly_connected_components(G),key=len)[-1]
    G_big = G.subgraph(big_cluster_nodes)
    
    # Create a position dictionary which takes the latitudes and longitudes. This allows for proper plotting later
    pos = {}
    for node in G_big.nodes():
        pos[node]=(G_big.node[node]['lon'], G_big.node[node]['lat'])
    
    return G_big, pos

def make_graph_with_poly(G, root, poly_vertices):
    # Given latitude/longitude vertices of the polygon,
    poly = Polygon(poly_vertices)
    
    # Check each note to determine if it lies within the polygon region.
    for node in root.iter('node'):
        lat_= float(node.attrib['lat'])
        lon_ = float(node.attrib['lon'])
        p_test = Point(lon_,lat_)
        inside = poly.contains(p_test)
        G.add_node(int(node.attrib['id']), lat=lat_, lon=lon_, is_inside = inside)
    
    # For each 'way' (a connected series of points in the osm format), check if it's a street (highway or oneway). Tags will be stored in the edge
    for way in root.iter('way'):
        highway = False
        oneway = False
        highway_tag = None
        for tag in way.iter('tag'):
            if 'highway' in tag.attrib['k']:
                highway = True
                highway_tag = tag.attrib.get('v',None)

            if 'oneway' in tag.attrib['k']:
                oneway = True

        # If it's a highway, we also store the speed in the edge, per highway type listed in highway_dict
        if highway:
            pairs = pairwise(way.iter('nd'))
            
            # Calculate distance
            if highway_tag in highway_dict:
                speed = highway_dict[highway_tag]
            else:
                speed = 30

            # Here we actually make the edge with the stored variables.       
            for a, b in pairs:
                aa, bb = int(a.attrib['ref']), int(b.attrib['ref'])
                dist = node_dist(G, aa, bb)
                time = float(dist)/(speed*1.60934)
                G.add_edge(aa, bb, tag=highway_tag, speed_tag=speed, dist_tag=dist, time_tag=time)

            # If twoway, must store reverse distance
            if not oneway:
                pairs = pairwise(way.iter('nd'))
                for a, b in pairs:
                    aa, bb = int(a.attrib['ref']), int(b.attrib['ref'])
                    dist = node_dist(G, bb, aa)
                    time = float(dist)/(speed*1.60934)
                    G.add_edge(bb, aa, tag=highway_tag, speed_tag=speed, dist_tag=dist, time_tag=time)

    return G

'''
6. Node matching
'''

def closest_node(new_lat,new_lon, G):
    # Given latitude and longitude, find the nearest node in graph G
    min_node = 0
    min_dist = float('inf')
    for node in G.nodes():
        lat = G.node[node]['lat']
        lon = G.node[node]['lon']

        dist = haversine_dist(new_lat, new_lon,lat,lon)
        if dist < min_dist:
            min_node = node
            min_dist = dist
    return (min_node, min_dist)

'''
7. Points-of-Interest (POI) injection
'''

def inject_multiple_pois(G, poi_dict, pos, poly_vertices):
    # Add point-of-interest nodes to graph
    # This now properly inserts a new poi node if the closest point is on the line, or it adds poi attributes to the node if an endpoint is the closest point
    # Output poi_node indicates either the new node id or the already-placed node id
    
    poi_node_dict = {}
    poly = Polygon(poly_vertices)
    for poi in poi_dict:# sample_pois:
        poi_name = poi
        poi_lat = poi_dict[poi]['lat']
        poi_lon = poi_dict[poi]['lon']
        poi_time = poi_dict[poi]['wait_time']
        poi_test = Point(poi_lon, poi_lat)
        inside = poly.contains(poi_test)
        
        # find the two nearest nodes and make a subgraph
        nearest_node1, nearest_node2, edge_dist = check_dist_to_edges(G, poi_lat, poi_lon)
        relevant_nodes = G.subgraph([nearest_node1,nearest_node2])

        # Create poi node

        x1 = G.node[nearest_node1]['lon']
        y1 = G.node[nearest_node1]['lat']
        x2 = G.node[nearest_node2]['lon']
        y2 = G.node[nearest_node2]['lat']

        # Use shapely.LineString to determine where the closest point on the line between the nodes lies
        edge_line = LineString([(x1,y1),(x2,y2)])
        poi_p = Point(poi_lon, poi_lat)

        # Will either interpolate onto the line or be one of the two nodes
        new_lon, new_lat = list(edge_line.interpolate(edge_line.project(poi_p)).coords)[0]

        # If node is already there (i.e. closest point on line segment is an endpoint), then just add POI tags to that node
        if (new_lon,new_lat) == (x1,y1) or (new_lon,new_lat) == (x2,y2):
            nearest_node, nearest_node_dist = closest_node(new_lat,new_lon, relevant_nodes)
            G.node[nearest_node]['poi_name'] = poi_name
            G.node[nearest_node]['wait_time'] = poi_time
            poi_node = nearest_node

        else: # Otherwise, need to add a node and adjust the edges
            pos[poi_name] = (new_lon,new_lat)
            G.add_node(poi_name, lat = new_lat, lon = new_lon, poi_name=poi_name, is_inside = inside, wait_time = poi_time)

            # Remove old edge
            two_way = False
            edge_dict = G.get_edge_data(nearest_node1, nearest_node2)
            highway_tag = edge_dict['tag']
            speed = edge_dict['speed_tag']
            dist1 = node_dist(G, nearest_node1, poi_name)
            time1 = float(dist1)/(speed*1.60934)
            dist2 = node_dist(G, poi_name, nearest_node2)
            time2 = float(dist2)/(speed*1.60934)

            # Remove edge (or two edges if twoway) between the initial node pair, to be replaced by two edges or pairs of edges
            G.remove_edge(nearest_node1, nearest_node2)
            if (nearest_node2, nearest_node1) in G.edges():
                G.remove_edge(nearest_node2, nearest_node1)
                two_way = True

            # Create new edge, calculate edge characteristics
            G.add_edge(nearest_node1, poi_name, tag=highway_tag, speed_tag=speed, dist_tag=dist1, time_tag=time1)
            G.add_edge(poi_name, nearest_node2, tag=highway_tag, speed_tag=speed, dist_tag=dist2, time_tag=time2)

            if two_way:
                dist1 = node_dist(G, nearest_node2, poi_name)
                time1 = float(dist1)/(speed*1.60934)
                dist2 = node_dist(G, poi_name, nearest_node1)
                time2 = float(dist2)/(speed*1.60934)
                G.add_edge(nearest_node2, poi_name, tag=highway_tag, speed_tag=speed, dist_tag=dist1, time_tag=time1)
                G.add_edge(poi_name, nearest_node1, tag=highway_tag, speed_tag=speed, dist_tag=dist2, time_tag=time2)

            poi_node = poi_name

        poi_node_dict[poi_name] = poi_node

    return G, poi_node_dict, pos

def check_dist_to_edges(G, lat, lon):
    # Calculate the distance to edges. Involved with selecting the correct point in POI injection
    edge_dist = float('inf')
    edge_node1 = -1
    edge_node2 = -1
    for edge in G.edges():
        node1 = edge[0]
        node2 = edge[1]
        lat1, lon1 = G.node[node1]['lat'], G.node[node1]['lon']
        lat2, lon2 = G.node[node2]['lat'], G.node[node2]['lon']
        check_dist = shapely_dist_point_to_edge(lon, lat, lon1, lat1, lon2, lat2)
        if check_dist < edge_dist:
            edge_dist = check_dist
            edge_node1 = node1
            edge_node2 = node2
    G.node[node1]['lat']
            
    return edge_node1, edge_node2, edge_dist

def shapely_dist_point_to_edge(poi_lon, poi_lat, x1, y1, x2, y2):
    # Uses shapely classes to determine distance from point to edge
    edge_line = LineString([(x1,y1),(x2,y2)])
    poi_p = Point(poi_lon, poi_lat)
    dist = edge_line.distance(poi_p)
    return dist    

'''
8. Path searching
'''

def find_path(G, source, target, initial_time, min_time, max_time, target_times, time_limit):
    # Given source node, perform pathfinding until find a node which carries 'True' for 'is_inside', meaning the node lies within the input polygon region
    # This is closer to a depth-first search. Due to the fact that the path must find the node within a certain time,
    # an optimized pathfinder algorithm like A* doesn't perform correctly. However, this should later be modified for use with
    # a greedier algorithm to improve runtime.
    runtime = [initial_time] # track runtime at each node
    visited = [source] # track current path
    stack = [iter(G[source])] # put children onto stack
    count = 0 # used for monitoring

    while stack:
        children = stack[-1]
        child = next(children, None) # iterate over children

        if child is None:
            # no children left
            stack.pop()
            visited.pop()
            runtime.pop()
        elif runtime[-1] <= max_time:
            # continue on path if under cutoff time
            time_from_prev = G.get_edge_data(visited[-1],child)['time_tag'] # time from previous to current node
            current_time = runtime[-1] + time_from_prev # total time to current node
            wait_time = G.node[child].get('wait_time', None)
            if wait_time:
                current_time += wait_time

            if child in target_times:
                time_to_target = target_times[child]
            else:
                time_to_target = 0
                
            inside_polygon = G.node[child].get('is_inside', False)


            # stop trying this path if can't reach target under max time or if reached time limit
            if (current_time <= max_time) and (current_time + time_to_target) <= time_limit:
                count += 1
                if count%1000000==0:
                    # monitor progress
                    print(count)
                
                if inside_polygon and current_time >= min_time:
                    # inside polygon with valid time
                    print(current_time)
                    print(count)
                    return visited + [child]
                elif child not in visited and child != target:
                    # go further if child has not been visited and is not target
                    stack.append(iter(G[child]))    
                    visited.append(child)
                    runtime.append(current_time)
        else:
            # stop if runtime greater than cutoff time
            print(runtime[-1])
            stack.pop()
            visited.pop()
            runtime.pop()

    return visited # did not find a valid path

'''
9. Utilities
'''

def pairwise(way_iter):
    # Simple utility for getting all adjacent element pairs in a list. For use with steps
    a, b = itertools.tee(way_iter)
    next(b, None)
    return itertools.izip(a,b)

def node_dist(G,node1, node2):
    # Haversine distance calculation.
    lon1,lat1 = G.node[node1]['lon'], G.node[node1]['lat']
    lon2,lat2 = G.node[node2]['lon'], G.node[node2]['lat']
    dist = haversine_dist(lat1,lon1,lat2,lon2)
    return dist

def haversine_dist(lat1,lon1,lat2,lon2):
    # In kilometers. Haversine is used for latitude and longitude distance calculations
    radius = 6372.8
    dLat = radians(lat2-lat1)
    dLon = radians(lon2-lon1)
    lat1 = radians(lat1)
    lat2 = radians(lat2)
    
    a = sin(dLat/2)**2 + cos(lat1)*cos(lat2)*sin(dLon/2)**2
    c = 2*asin(sqrt(a))
    return radius*c

def nodes_in_poly(G):
    # Simply search graph to find nodes that are inside the polygon
    poly_nodes = []
    for node in G.nodes():
        try:
            if G.node[node]['is_inside'] == True:
                poly_nodes.append(node)
        except:
            print(node)
    K = G.subgraph(poly_nodes)
    return K

def f7(seq):
    # Used for plotting verification to prevent excess overlaying
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]

def retrieve_traj_coords(G,traj_dict):
    # Get the positions of trajectory points. For plotting purposes
    coords_dict = {}
    
    for traj in traj_dict:
        coords_list = []
        for node in range(len(traj_dict[traj])):
            lat1,lon1 = traj_dict[traj][node].lat_start, traj_dict[traj][node].lon_start
            lat2,lon2 = traj_dict[traj][node].lat_end, traj_dict[traj][node].lon_end
            coords_list.append((lat1,lon1))
            coords_list.append((lat2,lon2))
        coords_dict[traj] = f7(coords_list)
    return coords_dict

'''
10. Plotting
'''

def graph_plotter(G, pos, valid_traj, traj_dict,xlim, ylim):
    # Plot the graph. Includes plotting the map of the region, as well as the trajectory points and polygon points
    P = nodes_in_poly(G)
    plt.figure(figsize=(10,10))
    nx.draw_networkx_edges(G, pos, edge_color = '0.75', alpha = 0.9)
    nx.draw(P,pos=pos,node_color='r', node_size=20, alpha=0.7)
    #fig, ax = plt.subplots()
    coords_dict = retrieve_traj_coords(G,traj_dict)
    #coords = list(set(retrieve_traj_coords(G,traj_dict)))
    
    cols = 'gmycw'
    markers = 'osv*D^'
    col_index = 0
    marker_index = 0
    print(coords_dict)
    for ind in coords_dict:
        
        traj_nodes = list(set(coords_dict[ind]))
        x, y = zip(*traj_nodes)
        print(traj_nodes)
        if ind in valid_traj:
            K = G.subgraph(valid_traj[ind])
            nx.draw(K,pos=pos,node_color=cols[col_index], node_size=50, alpha=0.2, marker = 'o')
            col_index+=1
        for j in range(len(y)):
            plt.scatter(y[j], x[j], color=cols[col_index], marker = markers[marker_index], s=200)  
            col_index+=1
        #plt.scatter(y,x,color = cols[col_index], marker = markers[marker_index], s=80)
        #plt.text(y,x,0)
        label = 1
        
        for coord in traj_nodes:
            v, w = coord[0], coord[1]
            print(v,w)
            #plt.text(v ,w, label) # This breaks plotting
            label+=1
        
        col_index +=1
        if col_index >=5:
            col_index = 0
            marker_index +=1
        
        if marker_index >=6:
            marker_index = 0
    plt.xlim(xlim)
    plt.ylim(ylim)

#if __name__ == "__main__":
#    traj_file = './files/trajectory.txt'
#    query_file = './files/query.txt'
#    poi_file = './files/POI.txt'
#    run_solution(traj_file, query_file, poi_file)

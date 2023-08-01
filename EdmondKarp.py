class Vertex:
    def __init__(self, label):
        """
        A vertex class which contains various properties necessary for breadth first search.

        Input:
            label: The identity of the vertex

        Space complexity: O(E) since it is possible for a single vertex to be associated with every edge
        """
        self.label = label
        self.edges = []
        self.previous = None
        self.visited = False
        self.discovered = False

    """
    Add an edge associated with this vertex 
    Time complexity: O(1), since append is O(1)
    """
    def add_edge(self, vertex_one, vertex_two, available_flow, index, exist=True):
        self.edges.append(Edge(vertex_one, vertex_two, available_flow, index, exist))


class Edge:
    def __init__(self, vertex_one: Vertex, vertex_two: Vertex, available_flow, index, exist=True):
        """
        An edge class

        Input:
            vertex_one: The vertex the edge starts from
            vertex_two: The vertex the edge ends at
            capacity: The max flow capacity of the edge
            reverse: True if it is a backward flow, False if it is a forward flow
            index: Index of the reverse edge in the edges array of another vertex
            exist: True if the edge exists, False if edge does not
        """
        self.vertex_one = vertex_one
        self.vertex_two = vertex_two
        self.available_flow = available_flow
        self.exist = exist
        self.index = index

    """
    Switch the boolean value for edge.exist. Used when available flow becomes 0 from more than 0 or becomes 
    more than 0 from 0.
    
    Time complexity: O(1)
    """
    def toggle_exist(self):
        self.exist = not self.exist

class ResidualNetwork:
    """
    Based on connections, maxIn and maxOut, initialise a residual network. For every forward edge flowing between
    any pair of vertices in a residual network, there would be another edge flowing in the opposite direction. The
    available flow of that edge would be equal to the current flow of its corresponding forward edge.

    Input:
        connections: The edges connecting any pair of data centres
        maxIn: A list of positive integers representing the maximum amount of data that can enter a data centre
               at once.
        maxOut: A list of positive integers representing the maximum amount of data that can leave a data centre
                at once.
        origin: An integer signifying the starting data centre where all data travels from
        targets: A list of integers signifying a list of data centres that are available for data storage

    Output: None, since nothing is returned

    Postcondition: A residual network is created

    Time complexity: O(|C| + |D|)

    Auxiliary space complexity: O(|C| + |D|)
    """
    def __init__(self, connections, maxIn, maxOut, origin, targets):
        self.sink = Vertex("T")
        self.connections = connections
        self.targets = targets
        self.max_vertex = len(maxIn) - 1

        # Generate a graph with vertex(0....3*max_vertex_number), where max_vertex_number is the vertex with\
        # the largest numerical label
        # Vertex(0...max_vertex_number) represents a data centre
        # Vertex(max_vertex_number + 1...2*max_vertex_number) represents the bottleneck going in of a particular data centre (based on maxIn)
        # Vertex(2*max_vertex_number + 1...3*max_vertex_number) represents the bottleneck going out of a particular data centre (based on maxOut)
        # Time complexity: O(|D|)
        self.residual_network = list(Vertex(num) for num in range(3 * (self.max_vertex + 1)))
        self.origin = self.residual_network[origin]

        # Last vertex in the residual network list is the super sink
        self.residual_network.append(self.sink)

        # For each vertex representing a bottleneck, connect them to their corresponding vertex. Also, connect a reverse edge
        # travelling in the opposite direction with 0 available flow for every edge.
        # Example for vertex(1), o--maxIn--o--maxOut--o, where the middle o represents vertex(1) and the left o
        # and right o represents the bottlenecks of vertex(1) determined by the maxIn and maxOut of vertex(1)
        # Store the index of each reverse edge for O(1) access
        # Time complexity: O(|D|), since there are D amount of maxIn and maxOut, where D is the number of data centres
        for index in range(len(maxIn)):
            self.residual_network[index + self.max_vertex + 1].add_edge(self.residual_network[index + self.max_vertex + 1], self.residual_network[index], maxIn[index], len(self.residual_network[index].edges))
            self.residual_network[index].add_edge(self.residual_network[index], self.residual_network[index + self.max_vertex + 1], 0, len(self.residual_network[index + self.max_vertex + 1].edges) - 1, False)
        for index in range(len(maxOut)):
            self.residual_network[index].add_edge(self.residual_network[index], self.residual_network[index + 2 * (self.max_vertex + 1)], maxOut[index], len(self.residual_network[index + 2 * (self.max_vertex + 1)].edges))
            self.residual_network[index + 2 * (self.max_vertex + 1)].add_edge(self.residual_network[index + 2 * (self.max_vertex + 1)], self.residual_network[index], 0, len(self.residual_network[index].edges) - 1, False)

        # Add edges to the graph, each edge representing a connection between two data centres.
        # Also, connect a reverse edge travelling in the opposite direction with 0 available flow for every edge
        # Store the index of each reverse edge for O(1) access
        # Time complexity: O(|C|)
        for connection in connections:
            self.residual_network[connection[0] + 2 * (self.max_vertex + 1)].add_edge(
                self.residual_network[connection[0] + 2 * (self.max_vertex + 1)],
                self.residual_network[connection[1] + self.max_vertex + 1], connection[2], len(self.residual_network[connection[1] + self.max_vertex + 1].edges))
            self.residual_network[connection[1] + self.max_vertex + 1].add_edge(
                self.residual_network[connection[1] + self.max_vertex + 1],
                self.residual_network[connection[0] + 2 * (self.max_vertex + 1)], 0, len(self.residual_network[connection[0] + 2 * (self.max_vertex + 1)].edges) - 1, False)

        # For each target vertex, connect them to the super sink. The capacity of this edge is based on the maxIn of the vertex
        # Also, connect a reverse edge travelling in the opposite direction with 0 available for every edge
        # Store the index of each reverse edge for O(1) access
        # Time complexity: O(T), where T is the number of targets. T < D as T cannot include origin
        for target in targets:
            self.residual_network[target].add_edge(self.residual_network[target], self.sink, maxIn[target], len(self.sink.edges))
            self.sink.add_edge(self.sink, self.residual_network[target], 0, len(self.residual_network[target].edges) - 1, False)

    """
    Used to find the shortest path from the start vertex to every other vertex, but in this case we are only
    interested in the shortest from the source vertex to the sink vertex. BFS is only applicable to graphs with
    no weighted edges.
    
    Input: 
       start: The starting vertex where the flow ends  
       end: The last vertex where the flow begins  
       
    Output: None since the function does not return anything
        
    Postcondition: All the vertex.previous is updated with an edge. Reversely follow the path from the end 
                   to start would give us the shortest path from start to end.
                   
    Time complexity: O(V + E) since worst case scenario requires all vertices and edges to be traversed
    
    Space complexity: O(V + E) since the graph is stored in an adjacency list 
    """
    def breadth_first_search(self, start: Vertex, end: Vertex):
        discovered = deque([])
        discovered.append(start)

        while len(discovered) > 0:
            # Start on the start vertex
            current_vertex = discovered.popleft()
            current_vertex.visited = True

            # If target is found, end the method
            if current_vertex == end:
                return

            for edge in current_vertex.edges:
                # Check if edge.available flow is 0 or not, if true ignore the edge
                if not edge.exist:
                    continue

                visiting_vertex = edge.vertex_two

                # If vertex that can be traversed to from current edge is not discovered yet, add it into
                # the discovered queue
                if not visiting_vertex.discovered:
                    discovered.append(visiting_vertex)
                    visiting_vertex.discovered = True
                    visiting_vertex.previous = edge

        return

    """
    Reversely generate the shortest path from start to end

    Input: 
       start: The starting vertex where the flow begins  
       end: The last vertex where the flow ends  

    Output: An array containing the shortest path represented by edges from start to end in reverse order

    Time complexity: O(V) since the path may potentially require traversing all of the vertices in the graph

    Auxiliary space complexity: O(V), an array containing the path which may potentially be all of the vertices in the graph
    """
    def backtracking(self, start: Vertex, end: Vertex):
        # Create a list with the first edge from the destination in it
        shortest_path = [end.previous]

        # Initialise the next vertex in the path
        next_vertex = end.previous.vertex_one

        # Keep appending edges that are part of the path from start to end until the start vertex is reached
        while next_vertex != start:
            if next_vertex is None:
                return shortest_path

            shortest_path.append(next_vertex.previous)
            next_vertex = next_vertex.previous.vertex_one

        return shortest_path

    """
    Check whether there is a path from source to sink
    
    Output: True, if there is path to the sink. Otherwise, False
    
    Time complexity: O(|C|) since it's possible that the path contains all the connections in the network
    """
    def has_augmenting_path(self) -> bool:
        self.breadth_first_search(self.origin, self.sink)
        return not (self.sink.previous is None)

    """
    Return the augmenting path. (In this case, the shortest path from source to sink since BFS is used)
    
    Output: A list of edge forming the shortest path from origin to sink 
    
    Time complexity: O(|C|) since it's possible that the path contains all the connections in the network
    """
    def get_augmenting_path(self) -> list:
        path = []
        inverted_path = self.backtracking(self.origin, self.sink)

        for i in range(len(inverted_path) - 1, -1, -1):
            path.append(inverted_path[i])

        return path

    """
    Return the critical edge's flow. An edge is critical if it has the smallest available flow among all the
    edges in the augmenting path.
    
    Input: 
        path: The shortest path from origin to sink
        
    Output: The residual capacity of the path, which is the lowest available flow among all edges forming the path
    
    Time complexity: O(|C|) since it's possible that the path contains all the connections in the network
    """
    def get_critical_flow(self, path) -> int:
        critical_edge_flow = math.inf

        for edge in path:
            if edge.available_flow < critical_edge_flow:
                critical_edge_flow = edge.available_flow

        return critical_edge_flow

    """
    Get the edge travelling in the opposite direction
    
    Input:
        vertex_one: the starting vertex where the edge begins from
        index: the index of the reverse edge stored in the edges array of vertex_one
    
    Output: An edge travelling in the opposite direction
    
    Time complexity: O(1)
    """
    def get_reverse_edge(self, vertex_one: Vertex, index):
        return vertex_one.edges[index]

    """
    Update the available flow of each edge, "removing" each edge with available flow of zero and "adding" edges
    with available flow more than zero after the augment. Reset the network for the next iteration of BFS as well.
    
    Input: 
        path: The shortest path from source to sink
        min_flow: The residual capacity of the path
        
    Output: None, since the method does not return anything
    
    Postcondition: An updated residual network 
    
    Time complexity: O(|C|)
    """
    def augment_flow(self, path, min_flow):
        # For each edge, update its available flow and its reverse counterpart available flow. Toggle the
        # edge exist attribute according to the result after update
        # Time complexity: O(|C|)
        for edge in path:
            edge.available_flow -= min_flow
            if edge.available_flow == 0:
                edge.toggle_exist()
            reverse_edge = self.get_reverse_edge(edge.vertex_two, edge.index)
            reverse_edge.available_flow += min_flow
            if reverse_edge.available_flow > 0:
                reverse_edge.toggle_exist()

        # Reset each vertex
        # Time complexity: O(|D|)
        for vertex in self.residual_network:
            vertex.previous = None
            vertex.visited = False
            vertex.discovered = False

"""
Ford fulkerson with BFS (edmond karp algorithm). The edmond karp algorithm seeks to improve on the time complexity
of ford fulkerson as ford fulkerson has a time complexity of O(|E|*f), where |E| is the number of edges and f
is the maximum flow the network. Edmond karp specifies using BFS to find the available path from source to sink, provided
a path exists. This is because BFS promises to always find the shortest path from source to sink at each iteration as 
opposed to other path finding algorithms like DFS. This changes the overall time complexity to O(|V||E|^2), decoupling
it from the maximum flow of the network. The proof of the time complexity involves two steps: 

1. Proving the shortest path from source to sink is monotonically increasing after every flow augmentation

2. Proving each edge has a maximum number of V times to be critical. An edge becomes critical if it has the 
   smallest available flow among all the edges in the augmenting path. 
   
Since there are E edges that can become critical V times, there would be VE iterations of edmond karp. Each iteration takes
O(V+E) time due to BFS. Hence the overall time complexity is O(V^2E + VE^2), which can simplify to O(VE^2) as V is mostly
smaller or equal to E.


Input:
    connections: The edges connecting any pair of data centres
    maxIn: A list of positive integers representing the maximum amount of data that can enter a data centre
           at once. 
    maxOut: A list of positive integers representing the maximum amount of data that can leave a data centre
            at once.
    origin: An integer signifying the starting data centre where all data travels from
    targets: A list of integers signifying a list of data centres that are available for data storage
    
Output: The maximum flow of the graph

Time complexity: O(|D||C|^2), this time complexity can be proven through the edmond karp algorithm 
                 (ford fulkerson with BFS), where |D| is the number of data centres and |C| is the 
                 number of connections. 

Auxiliary space complexity: O(V + E), the space complexity of the residual network
"""
def ford_fulkerson(connections, maxIn, maxOut, origin, targets):
    flow = 0

    residual_network = ResidualNetwork(connections, maxIn, maxOut, origin, targets)
    while residual_network.has_augmenting_path():
        path = residual_network.get_augmenting_path()
        min_flow = residual_network.get_critical_flow(path)
        flow += min_flow
        residual_network.augment_flow(path, min_flow)

    return flow

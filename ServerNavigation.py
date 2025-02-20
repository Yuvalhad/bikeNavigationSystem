import requests
import networkx as nx
import matplotlib.pyplot as plt

def fetch_osm_data_bbox(min_lat, min_lon, max_lat, max_lon):
    """
    Fetches road network data using a bounding box.
    """
    query = f"""
    [out:json];
    (
      way["highway"]({min_lat},{min_lon},{max_lat},{max_lon});
    );
    out geom;
    """
    response = requests.get("http://overpass-api.de/api/interpreter", params={"data": query})

    if response.status_code == 200:
        return response.json()
    else:
        print("❌ Failed to fetch data:", response.text)
        return None


def osm_to_graph(osm_data):
    """
    Converts OSM road data into a NetworkX graph.
    """
    G = nx.Graph()

    if "elements" not in osm_data:
        print("❌ No road data found!")
        return G

    for element in osm_data["elements"]:
        if element["type"] == "way" and "geometry" in element:
            street_name = element.get("tags",{}).get("name","unknown")
            nodes = []
            #need to create new dict and add to there the node and the name of the street of these nodes
            #check how can i use the dict in simplify_graph
            for coord in element["geometry"]:
                node = (coord["lon"], coord["lat"])
                if node not in G:
                  G.add_node(node, streets = set())
                G.nodes[node]["streets"].add(street_name)
                nodes.append(node)

            # Create edges between consecutive nodes
            for i in range(len(nodes) - 1):
                G.add_edge(nodes[i], nodes[i + 1])

            for i in range(len(nodes) - 1):
                if G.has_edge(nodes[i], nodes[i + 1]):
                    # בדיקה אם מאפיין "streets" קיים עבור הקשת, ואם לא – יצירתו
                    if "streets" not in G.edges[nodes[i], nodes[i + 1]]:
                        G.edges[nodes[i], nodes[i + 1]]["streets"] = set()
                    G.edges[nodes[i], nodes[i + 1]]["streets"].add(street_name)
                else:
                    G.add_edge(nodes[i], nodes[i + 1], streets={street_name})
                    G.edges[nodes[i], nodes[i + 1]]["streets"].add(street_name)

    return G


def simplify_graph(G):
    """
    Removes redundant nodes that are not intersections or dead ends.
    """
    to_remove = []

    for node in list(G.nodes):
        if G.degree(node) == 2:  # Node is not an intersection
            neighbors = list(G.neighbors(node))
            # i will choose the street name of the neighbors
            edge1_streets = G.edges[node, neighbors[0]]["streets"]
            edge2_streets = G.edges[node, neighbors[1]]["streets"]

            #check if there is a intersection (חיתוך) and if yes then i will need to delete one
            common_streets = edge1_streets.intersection(edge2_streets)
            if common_streets:
                #Merge segment
                G.add_edge(neighbors[0], neighbors[1], streets=common_streets)
                to_remove.append(node)
        elif G.degree(node) == 1:
            #need to add a new edge from the back neighbor to the front neighbor befor deletening
            to_remove.append(node)

    G.remove_nodes_from(to_remove)  # Remove redundant nodes
    return G


def plot_graph(G):
    """
    Plots the road network graph.
    """
    if len(G.nodes) == 0:
        print("❌ No nodes to plot.")
        return

    plt.figure(figsize=(10, 10))
    pos = {node: (node[0], node[1]) for node in G.nodes}
    nx.draw(G, pos, node_size=10, edge_color="blue", with_labels=False)
    plt.title("Road Network Graph")
    plt.show()

# Example: Ashkelon, Israel bounding box
min_lat, min_lon = 31.65, 34.55
max_lat, max_lon = 31.70, 34.60

# Example: Tel Aviv, Israel bounding box
#min_lat, min_lon = 32.0150, 34.7450
#max_lat, max_lon = 32.1250, 34.8550

# Example: Israel bounding box
#min_lat, min_lon = 29.5, 34.2
#max_lat, max_lon = 33.3, 35.9

osm_data = fetch_osm_data_bbox(min_lat, min_lon, max_lat, max_lon)
print(osm_data)  # Check if roads exist

if osm_data:
    # Convert OSM data to a graph
    road_graph = osm_to_graph(osm_data)

    # Simplify the graph
    simplified_graph = simplify_graph(road_graph)

    # Plot the graph
    plot_graph(simplified_graph)


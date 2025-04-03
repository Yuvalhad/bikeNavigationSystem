import requests
import networkx as nx
import matplotlib.pyplot as plt
from math import sqrt


def fetch_osm_data_bbox(min_lat, min_lon, max_lat, max_lon):
    """
    Fetches road network data using a bounding box.
    Focus on highways suitable for cycling (e.g., secondary, tertiary, cycleway, residential).
    """
    query = f"""
    [out:json];
    (
      way["highway"~"secondary|tertiary|cycleway|residential|unclassified"]({min_lat},{min_lon},{max_lat},{max_lon});
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
            street_name = element.get("tags", {}).get("name", "unknown")
            highway_type = element.get("tags", {}).get("highway", "unknown")
            nodes = []
            for coord in element["geometry"]:
                node = (coord["lon"], coord["lat"])
                if node not in G:
                    G.add_node(node, streets=set(), highway=highway_type)
                G.nodes[node]["streets"].add(street_name)
                nodes.append(node)

            for i in range(len(nodes) - 1):
                G.add_edge(nodes[i], nodes[i + 1])

            for i in range(len(nodes) - 1):
                if G.has_edge(nodes[i], nodes[i + 1]):
                    if "streets" not in G.edges[nodes[i], nodes[i + 1]]:
                        G.edges[nodes[i], nodes[i + 1]]["streets"] = set()
                    G.edges[nodes[i], nodes[i + 1]]["streets"].add(street_name)
                else:
                    G.add_edge(nodes[i], nodes[i + 1], streets={street_name})

    return G


def simplify_graph(G, target_nodes=300):
    """
    Simplifies the graph to reach a target number of nodes (e.g., 300), keeping only key junctions (degree > 2)
    or endpoints (degree = 1), and merging nodes on the same road based on distance.
    """
    to_remove = []
    remaining_nodes = G.number_of_nodes()
    distance_threshold = 0.0001  # סף מרחק נמוך מאוד להתחלה, להתאמה דלילה

    while remaining_nodes > target_nodes:
        print(f"Current nodes: {remaining_nodes}, targeting {target_nodes}")
        removed_in_this_pass = False

        for node in list(G.nodes):
            degree = G.degree(node)
            if degree == 2:  # נודים באמצע כביש
                neighbors = list(G.neighbors(node))  # המרת השכנים לרשימה
                if len(neighbors) != 2:
                    continue  # דילוג אם יש בעיה במבנה
                edge1_streets = G.edges[node, neighbors[0]]["streets"]
                edge2_streets = G.edges[node, neighbors[1]]["streets"]
                common_streets = edge1_streets.intersection(edge2_streets)

                if common_streets:
                    # בדיקת מרחק בין השכנים
                    pos1 = neighbors[0]
                    pos2 = neighbors[1]
                    dist = sqrt((pos2[0] - pos1[0]) ** 2 + (pos2[1] - pos1[1]) ** 2)
                    if dist < distance_threshold:
                        G.add_edge(pos1, pos2, streets=common_streets)
                        to_remove.append(node)
                        removed_in_this_pass = True

        if not removed_in_this_pass:
            print("⚠️ No more nodes to remove with current distance threshold. Increasing threshold...")
            distance_threshold *= 2  # מגדיל את הסף אם אין עוד נודים למחיקה
            if distance_threshold > 0.05:  # הגבלה על הגדלת הסף
                print("⚠️ Maximum distance threshold reached. Stopping simplification.")
                break

        G.remove_nodes_from(to_remove)
        remaining_nodes = G.number_of_nodes()
        print(f"Removed {len(to_remove)} nodes. Remaining: {remaining_nodes}")
        to_remove = []  # איפוס הרשימה למחזור הבא

    # וידוא שהגרף נשאר מחובר (אם לא, נוסיף חיבורים נוספים)
    if not nx.is_connected(G):
        largest_component = max(nx.connected_components(G), key=len)
        G = G.subgraph(largest_component).copy()

    # הסרת קצוות לא חיוניים (degree = 1) אם אין להם חיבור משמעותי
    to_remove_endpoints = []
    for node in G.nodes:
        if G.degree(node) == 1:
            neighbors = list(G.neighbors(node))  # המרת השכנים לרשימה
            if len(neighbors) == 1:
                neighbor = neighbors[0]
                if G.degree(neighbor) > 2:  # אם השכן הוא צומת, נשאיר
                    continue
                to_remove_endpoints.append(node)
    G.remove_nodes_from(to_remove_endpoints)

    print(
        f"✅ Final simplification: Removed {G.number_of_nodes() - remaining_nodes} nodes. Remaining nodes: {G.number_of_nodes()}")
    return G


def plot_graph(G):
    """
    Plots the road network graph with improved visualization, similar to the provided image.
    """
    if len(G.nodes) == 0:
        print("❌ No nodes to plot.")
        return

    plt.figure(figsize=(12, 12))  # גודל גדול יותר לראות טוב יותר
    pos = {node: (node[0], node[1]) for node in G.nodes}
    nx.draw(G, pos, node_size=10, edge_color="blue", node_color="blue", with_labels=False)
    plt.title("Road Network Graph (Minimal Nodes for Cycling Navigation)")
    plt.show()


# Example: Ashkelon, Israel bounding box
min_lat, min_lon = 31.65, 34.55
max_lat, max_lon = 31.70, 34.60

osm_data = fetch_osm_data_bbox(min_lat, min_lon, max_lat, max_lon)
print(osm_data)  # Check if roads exist

if osm_data:
    # Convert OSM data to a graph
    road_graph = osm_to_graph(osm_data)

    # Simplify the graph to reach ~300 nodes (adjustable)
    simplified_graph = simplify_graph(road_graph, target_nodes=300)

    # Plot the graph
    plot_graph(simplified_graph)
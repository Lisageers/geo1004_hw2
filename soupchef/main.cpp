#include <iostream>
#include <fstream>
#include <sstream>
#include <list>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include "DCEL.hpp"
#include <stack>    // Constructs a stack container adaptor object


// forward declarations; these functions are given below main()
void printDCEL(DCEL & D);

typedef std::pair<int,int> pair;
struct pair_hash
{
    template <class T1, class T2>
    std::size_t operator() (const std::pair<T1, T2> &pair) const
    {
        return std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
    }
};

std::vector<double> cross_product(std::vector<double> v0, std::vector<double> v1)
{
    double vx = v0[1] * v1[2] - v0[2] * v1[1];
    double vy = v0[2] * v1[0] - v0[0] * v1[2];
    double vz = v0[0] * v1[1] - v0[1] * v1[0];
    return std::vector<double> {vx, vy, vz};
}
float dot_product(std::vector<double> v0, std::vector<double> v1)
{
    return (v0[0] * v1[0] + v0[1] * v1[1] + v0[2] * v1[2]);
}

float signed_volume(Vertex* &a, Vertex* &b, Vertex* &c, std::vector<double> &d)
{
    std::vector<double> ad{a->x-d[0], a->y-d[1], a->z-d[2]};
    std::vector<double> bd{b->x-d[0], b->y-d[1], b->z-d[2]};
    std::vector<double> cd{c->x-d[0], c->y-d[1], c->z-d[2]};
    std::vector<double> crossproduct = cross_product(bd, cd);
    float volume =  dot_product(ad, crossproduct);
    return volume;
}

// computes the distance of a point to a plane
float compute_distance(Vertex* &a, Vertex* &b, Vertex* &c, std::vector<double> &p)
{
    // calculate plane
    std::vector<double> ba{b->x-a->x, b->y-a->y, b->z-a->z};
    std::vector<double> ca{c->x-a->x, c->y-a->y, c->z-a->z};
    std::vector<double> a_double{a->x, a->y, a->z};
    std::vector<double> n = cross_product(ba, ca);
    std::vector<double> plane {n[0], n[1], n[2], -dot_product(n, a_double)};

    // calculate distance
    float d = (plane[0] * p[0]) + (plane[1]* p[1]) + (plane[2] * p[2]) + plane[3];
    float e = sqrt(pow(plane[0], 2) + pow(plane[1], 2) + pow(plane[2], 2));
    return abs(d/e);
}

void orient_faces(HalfEdge* e)
{
    Face* f = e->incidentFace;
    HalfEdge* e0 = f->exteriorEdge;
    HalfEdge* e1 = f->exteriorEdge->next;
    HalfEdge* e2 = f->exteriorEdge->next->next;

    // if the twin edges have the same origin, orientation of the twin is wrong
    if (e->origin == e->twin->origin)
    {
        // switch vertices first edge
        Vertex* originholder = e0->origin;
        e0->origin = e0->destination;
        e0->destination = originholder;

        // switch vertices second edge
        originholder = e1->origin;
        e1->origin = e1->destination;
        e1->destination = originholder;;

        // switch vertices third edge
        originholder = e2->origin;
        e2->origin = e2->destination;
        e2->destination = originholder;

        // switch order first edge
        HalfEdge* orderholder = e0->prev;
        e0->prev = e0->next;
        e0->next = orderholder;

        // switch order second edge
        orderholder = e1->prev;
        e1->prev = e1->next;
        e1->next = orderholder;

        // switch order first edge
        orderholder = e2->prev;
        e2->prev = e2->next;
        e2->next = orderholder;
    }
    // check holes
    if (f->holes.size() > 1)
    {
        const auto & holes = f->holes;
        for (const auto & h : holes )
        {
            std::vector<double> point{0,0,0};
            float volume_face = signed_volume(e0->origin, e0->destination, e0->next->destination, point);
            float volume_hole = signed_volume(h->origin, h->destination, h->next->destination, point);

            // orientation hole is incorrect if same as face
            if ((volume_face > 0 and volume_hole > 0) or (volume_face < 0 and volume_hole < 0))
            {
                HalfEdge* h0 = h;
                HalfEdge* h1 = h->next;
                HalfEdge* h2 = h->next->next;
                // switch vertices first edge
                Vertex* originholder = h0->origin;
                h0->origin = h0->destination;
                h0->destination = originholder;

                // switch vertices second edge
                originholder = h1->origin;
                h1->origin = h1->destination;
                h1->destination = originholder;;

                // switch vertices third edge
                originholder = h2->origin;
                h2->origin = h2->destination;
                h2->destination = originholder;

                // switch order first edge
                HalfEdge* orderholder = h0->prev;
                h0->prev = h0->next;
                h0->next = orderholder;

                // switch order second edge
                orderholder = h1->prev;
                h1->prev = h1->next;
                h1->next = orderholder;

                // switch order first edge
                orderholder = h2->prev;
                h2->prev = h2->next;
                h2->next = orderholder;
            }
        }
    }

}

/* 
  Example functions that you could implement. But you are 
  free to organise/modify the code however you want.
  After each function you should have a DCEL without invalid elements!
*/
// 1.
std::unordered_map< Vertex*, int> importOBJ(DCEL & D, const char *file_in) {
    std::string line;
    std::ifstream input(file_in);
    double x, y, z;
    int v0, v1, v2;
    int linecounter = 0;
    std::unordered_map<int, Vertex*> vertices;
    std::unordered_map<Vertex*, int> vertices_reverse;
    std::unordered_map<pair,std::vector<HalfEdge*>,pair_hash> edge_points;

    // read lines
    while (std::getline(input, line)) {
        linecounter++;
        // if line starts with v
        if (line[0] == 'v') {
            // split line at spaces
            std::istringstream ss(line);
            // move to second element
            std::string word;
            ss >> word;
            ss >> word;
            // create x with string converted to double
            x = std::stof(word);
            // move one further to get y
            ss >> word;
            y = std::stof(word);
            // move one further to get z
            ss >> word;
            z = std::stof(word);
            // put vertex point in vector
            Vertex* v = D.createVertex(x, y, z);
            vertices[linecounter] = v;
            vertices_reverse[v] = linecounter;
        }
            // if line starts with f
        else if (line[0] == 'f')
        {
            // split line at spaces
            std::istringstream ss(line);
            // initialise
            std::string word;
            // move two to get first vertex
            ss >> word;
            ss >> word;
            v0 = std::stoi(word);
            // move one further to get second
            ss >> word;
            v1 = std::stoi(word);
            // move one further to get third
            ss >> word;
            v2 = std::stoi(word);

            // dcel initialise
            HalfEdge* e0 = D.createHalfEdge();
            HalfEdge* e1 = D.createHalfEdge();
            HalfEdge* e2 = D.createHalfEdge();
            Face* f = D.createFace();
            f->exteriorEdge = e0;

            // create dcel for first edge
            e0->origin = vertices[v0];
            e0->destination = vertices[v1];
            e0->next = e1;
            e0->prev = e2;
            e0->incidentFace = f;
            // put edge and points in dict
            pair  edge_pair = {std::min(v0, v1), std::max(v0, v1)};
            // if there is an edge in dict
            if (!edge_points[edge_pair].empty())
            {
                edge_points[edge_pair].push_back(e0);
            }
            else // if there is no edge yet
            {
                edge_points.emplace(edge_pair, std::vector<HalfEdge*>());
                edge_points[edge_pair].push_back(e0);
            }

            // create dcel for second edge
            e1->origin = vertices[v1];
            e1->destination = vertices[v2];
            e1->next = e2;
            e1->prev = e0;
            e1->incidentFace = f;
            // put edge and points in dict
            edge_pair = {std::min(v1, v2), std::max(v1, v2)};
            // if there is an edge in dict
            if (!edge_points[edge_pair].empty())
            {
                edge_points[edge_pair].push_back(e1);
            }
            else // if there is no edge yet
            {
                edge_points.emplace(edge_pair, std::vector<HalfEdge*>());
                edge_points[edge_pair].push_back(e1);
            }

            // create dcel for third edge
            e2->origin = vertices[v2];
            e2->destination = vertices[v0];
            e2->next = e0;
            e2->prev = e1;
            e2->incidentFace = f;
            // put edge and points in dict
            edge_pair = {std::min(v0, v2), std::max(v0, v2)};
            // if there is an edge in dict
            if (!edge_points[edge_pair].empty())
            {
                edge_points[edge_pair].push_back(e2);
            }
            else // if there is no edge yet
            {
                edge_points.emplace(edge_pair, std::vector<HalfEdge*>());
                edge_points[edge_pair].push_back(e2);
            }
        }
    }
    input.close();
    // link twin edges
    for (auto const &values: edge_points)
    {
        values.second[0]->twin = values.second[1];
        values.second[1]->twin = values.second[0];
    }
    return vertices_reverse;
}

// 2.
std::unordered_map< Face*, int> groupTriangles(DCEL & D) {
    std::stack<Face*> faceStack;
    std::list<Face*> allfacesList;
    std::list<Face*> meshList;
    std::unordered_map<Face*, int> mesh_list;
    int feature = 1;
    bool map_complete = false;
    for (auto const &current_face:D.faces()) {
        if (std::find(allfacesList.begin(), allfacesList.end(), current_face.get()) != allfacesList.end()) {
            continue;
        }
        else {
            faceStack.push(current_face.get());
        }
        meshList.clear();

        while (!faceStack.empty()) {
            // Accessing the next element in the stack of faces
            Face* new_face = faceStack.top();
            // Removing the top element from the stack
            faceStack.pop();
            // Constructing a vector for the current_face's three edges by accessing the exteriorEdge and next pointers
            std::vector<HalfEdge*> face_vector = {new_face->exteriorEdge, new_face->exteriorEdge->next, new_face->exteriorEdge->next->next};
            for(auto const &current_edge: face_vector) {
                if (std::find(allfacesList.begin(), allfacesList.end(), (current_edge->twin->incidentFace)) != allfacesList.end()) {
                    continue;
                }
                else {
                    // Add the current edge to the stack of faces, the list of all faces, and the mesh list.
                    faceStack.push(current_edge->twin->incidentFace);
                    allfacesList.push_back(current_edge->twin->incidentFace);
                    meshList.push_back(current_edge->twin->incidentFace);
                }
            }
        }
        // infiniteFace stores all holes in the mesh
        // If current_face stack is empty, then file the remaining interior holes by adding one of the edges of each interior hole, to infiniteFace
        D.infiniteFace()->holes.push_back(current_face->exteriorEdge);
    }
    return mesh_list;
}


// 3.
void orientMeshes(DCEL & D) {
    // find lower boundary mesh
    // initialise min values
    float minx = INFINITY;
    float miny = INFINITY;
    float minz = INFINITY;
    const auto & faces = D.faces();
    for (const auto & f : faces) {
        // take vertices from face
        Vertex* v0 = f->exteriorEdge->origin;
        Vertex* v1 = f->exteriorEdge->destination;
        Vertex* v2 = f->exteriorEdge->next->destination;
        // update min if one of the vertices are smaller
        if (std::min({v0->x, v1->x, v2->x}) < minx)
        {
            minx = std::min({v0->x, v1->x, v2->x});
        }
        if (std::min({v0->y, v1->y, v2->y}) < miny)
        {
            miny = std::min({v0->y, v1->y, v2->y});
        }
        if (std::min({v0->z, v1->z, v2->z}) < minz)
        {
            minz = std::min({v0->z, v1->z, v2->z});
        }
    }

    // calculate point that is outside lower bounding box
    std::vector<double> exterior_point{minx-5, miny-5, minz-5};

    // find face that is closest to exterior point
    float min_distance = INFINITY;
    Face* closest_face;
    // find triangle with smallest distance to point
    for (const auto &f: faces) {
        // take vertices from face
        Vertex* v0 = f->exteriorEdge->origin;
        Vertex* v1 = f->exteriorEdge->destination;
        Vertex* v2 = f->exteriorEdge->next->destination;
        float distance = compute_distance(v0, v1, v2, exterior_point);
        if (distance < min_distance)
        {
            min_distance = distance;
            closest_face = f->exteriorEdge->incidentFace;
        }
    }

    // take vertices from closest face
    Vertex* v0 = closest_face->exteriorEdge->origin;
    Vertex* v1 = closest_face->exteriorEdge->destination;
    Vertex* v2 = closest_face->exteriorEdge->next->destination;
    HalfEdge* e0 = closest_face->exteriorEdge;
    HalfEdge* e1 = closest_face->exteriorEdge->next;
    HalfEdge* e2 = closest_face->exteriorEdge->next->next;
    float volume = signed_volume(v0, v1, v2, exterior_point);
    // if orientation is not correct
    if (volume < 0) {
        // switch vertices first edge
        Vertex* originholder = e0->origin;
        e0->origin = e0->destination;
        e0->destination = originholder;

        // switch vertices second edge
        originholder = e1->origin;
        e1->origin = e1->destination;
        e1->destination = originholder;;

        // switch vertices third edge
        originholder = e2->origin;
        e2->origin = e2->destination;
        e2->destination = originholder;

        // switch order first edge
        HalfEdge* orderholder = e0->prev;
        e0->prev = e0->next;
        e0->next = orderholder;

        // switch order second edge
        orderholder = e1->prev;
        e1->prev = e1->next;
        e1->next = orderholder;

        // switch order first edge
        orderholder = e2->prev;
        e2->prev = e2->next;
        e2->next = orderholder;
    }

    std::stack<HalfEdge*> edgeStack;
    std::vector<Face*> visited_faces;
    edgeStack.push(e0->twin);
    edgeStack.push(e1->twin);
    edgeStack.push(e2->twin);
    visited_faces.push_back(e0->incidentFace);
    visited_faces.push_back(e0->twin->incidentFace);
    visited_faces.push_back(e1->twin->incidentFace);
    visited_faces.push_back(e2->twin->incidentFace);
    while(!edgeStack.empty()) {
        HalfEdge* e = edgeStack.top();
        edgeStack.pop();
        orient_faces(e);
        bool twin0_visited = false;
        bool twin1_visited = false;
        bool twin2_visited = false;
        // check if twin faces are visited
        for (int i=0; i<visited_faces.size();i++) {
            if (e->twin->incidentFace == visited_faces[i])
            {
                twin0_visited = true;
            }
            else if (e->next->twin->incidentFace == visited_faces[i])
            {
                twin1_visited = true;
            }
            else if (e->next->next->twin->incidentFace == visited_faces[i])
            {
                twin2_visited = true;
            }
        }
        // update stack if face needs to be visited
        if (!twin0_visited) {
            visited_faces.push_back(e->twin->incidentFace);
            edgeStack.push(e->twin);
        }
        if (!twin1_visited) {
            visited_faces.push_back(e->next->twin->incidentFace);
            edgeStack.push(e->next->twin);
        }
        if (!twin2_visited) {
            visited_faces.push_back(e->next->next->twin->incidentFace);
            edgeStack.push(e->next->next->twin);
        }
    }
}

// 4.
void mergeCoPlanarFaces(DCEL & D) {
  // to do
  // merge adjacent triangles that are co-planar into larger polygonal faces

  // Approach: Where the normals are the same and there is a twin edge

}
// 5.
void exportCityJSON(DCEL & D, const char *file_out, std::unordered_map<Vertex*, int> verticesdict) {
    std::ofstream myfile(file_out);
    myfile << "{\n";
    myfile << "\t\"type\"" << ":" << " \"CityJSON\"" << ",\n";
    myfile << "\t\"version\"" << ":" << " \"1.0\"" << ",\n";

    myfile << "\t\"CityObjects\"" << ":" << " {\n";
    myfile << "\t\t\"id-1\"" << ":" << " {\n";
    myfile << "\t\t\t\"type\"" << ":" << " \"Building\"" << ",\n";
    myfile << "\t\t\t\"attributes\"" << ":" << " {},\n";
    myfile << "\t\t\t\"geometry\"" << ":" <<" [],\n";
    myfile << "\t\t\t\"children\"" << ":" <<" [\n";
    for (int i=0;i<D.faces().size();i++)
    {
        myfile << "\t\t\t\t\"id-" << (i+2) << "\"";
        if (i < D.faces().size()-1)
        {
            myfile << ",\n";
        }
        else
        {
            myfile << "]\n";
        }
    }
    myfile << "\t\t},\n";

    // write faces
    int loopcounterf = 0;
    const auto & faces = D.faces();
    for ( const auto & f : faces )
    {
        loopcounterf++;
        myfile << "\t\t\"id-" << (loopcounterf+1) << "\"" << ":" << " {\n";
        myfile << "\t\t\t\"type\"" << ":" << " \"BuildingPart\"" << ",\n";
        myfile << "\t\t\t\"parents\"" << ":" << " \"id-1\"" << ",\n";

        myfile << "\t\t\t\"geometry\"" << ":" <<" [{\n";
        myfile << "\t\t\t\t\"type\"" << ":" << " \"MultiSurface\"" << ",\n";
        myfile << "\t\t\t\t\"lod\"" << ":" << " \"2\"" << ",\n";
        myfile << "\t\t\t\t\"boundaries\"" << ":" << "[\n";
        Vertex* origin = f->exteriorEdge->origin;
        Vertex* next_vertex = f->exteriorEdge->destination;
        HalfEdge* next_edge = f->exteriorEdge;
        myfile << "\t\t\t\t[[" << (verticesdict[origin]-1);
        myfile << ", " << (verticesdict[next_vertex]-1);
        // loop through edges
        while (true)
        {
            next_edge = next_edge->next;
            next_vertex = next_edge->destination;
            if (next_vertex != origin)
            {
                myfile << ", " << (verticesdict[next_vertex]-1);
            }
            else // stop when loop is made
            {
                myfile << "]";
                break;
            }
        }
        // write holes
        const auto & holes = f->holes;
        for (const auto & h : holes ) {
            Vertex *origin = h->origin;
            Vertex *next_vertex = h->destination;
            HalfEdge *next_edge = h;
            myfile << ", [" << (verticesdict[origin] - 1);
            myfile << ", " << (verticesdict[next_vertex] - 1);
            while (true) {
                next_edge = next_edge->next;
                next_vertex = next_edge->destination;
                if (next_vertex != origin)
                {
                    myfile << ", " << (verticesdict[next_vertex] - 1);
                }
                else // stop when loop is made
                {
                    myfile << "]";
                    break;
                }
            }
        }
        myfile << "]\n\t\t\t\t]\n" << "\t\t\t}]\n";

        if (loopcounterf < faces.size())
        {
            myfile << "\t\t},\n";
        }
        else
        {
            myfile << "\t\t}\n";
        }
    }

    // write vertices
    myfile << "\t},\n";
    myfile << "\t\"vertices\"" << ":" << " [\n";
    int loopcounterv = 0;
    const auto &vertices = D.vertices();
    for ( const auto & v : vertices ) {
        loopcounterv++;
        myfile << "\t\t[" << v->x << ", " << v->y << ", " << v->z;
        if (loopcounterv < vertices.size())
        {
            myfile << "],\n";
        }
        else
        {
            myfile << "]\n";
        }
    }
    myfile << "\t]\n";
    myfile << "}\n";

    myfile.close();
}

int main(int argc, const char * argv[]) {
  const char *file_in = "../../bk_soup.obj";
  const char *file_out = "../../bk.json";

  // create an empty DCEL
  DCEL D;

  // 1. read the triangle soup from the OBJ input file and convert it to the DCEL,
    std::unordered_map<Vertex*, int> verticesdict = importOBJ(D, file_in);

  // 2. group the triangles into meshes,
    std::unordered_map<Face*, int> meshes_map = groupTriangles(D);

    // 3. determine the correct orientation for each mesh and ensure all its triangles
  //    are consistent with this correct orientation (ie. all the triangle normals 
  //    are pointing outwards).
    orientMeshes(D);

  // 4. merge adjacent triangles that are co-planar into larger polygonal faces.
  
  // 5. write the meshes with their faces to a valid CityJSON output file.
    exportCityJSON(D, file_out, verticesdict);

  return 0;
}


void printDCEL(DCEL & D) {

  // Quick check if there is an invalid element
  auto element = D.findInValid();
  if ( element == nullptr ) {
    // Beware that a 'valid' DCEL here only means there are no dangling links and no elimated elements.
    // There could still be problems like links that point to the wrong element.
    std::cout << "DCEL is valid\n";
  } else {
    std::cout << "DCEL is NOT valid ---> ";
    std::cout << *element << "\n";
  }

  // iterate all elements of the DCEL and print the info for each element
  const auto & vertices = D.vertices();
  const auto & halfEdges = D.halfEdges();
  const auto & faces = D.faces();
  std::cout << "DCEL has:\n";
  std::cout << " " << vertices.size() << " vertices:\n";
  for ( const auto & v : vertices ) {
    std::cout << "  * " << *v << "\n";
  }
  std::cout << " " << halfEdges.size() << " half-edges:\n";
  for ( const auto & e : halfEdges ) {
    std::cout << "  * " << *e << "\n";
  }
  std::cout << " " << faces.size() << " faces:\n";
  for ( const auto & f : faces ) {
    std::cout << "  * " << *f << "\n";
  }

}
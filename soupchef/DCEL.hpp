#pragma once

#include "DCELElements.hpp"

using VertexList    = std::list< std::unique_ptr<Vertex> >;
using HalfEdgeList  = std::list< std::unique_ptr<HalfEdge> >;
using FaceList      = std::list< std::unique_ptr<Face> >;

/*
  The class DCEL is a complete DCEL implementation.
*/
class DCEL {

  /* 
    maintain a list for each element type. All the Vertices, HalfEdges and
    Faces (except the infiniteFace) that are created are put in these lists.
  */
  VertexList _vertices;
  HalfEdgeList _halfEdges;
  FaceList _faces;

  /*
    The infiniteFace can be used
    1) as an incidentFace to any half-edge that does not have a real incident face
       (eg. in the triangle example in main.cpp:DemoDCEL()). In this way we prevent
       a dangling link. Notice that in a 2-manifold surface there should be no half-edge 
       that links to the infiniteFace (otherwise it would not be 2-manifold).
    2) To store the meshes contained in this DCEL as its holes. Simply push one half-edge
       for each mesh to the holes list of the infiniteFace. If you do this you can easily
       find your meshes by simply iterating the infiniteFace's holes list.

    The infinite Face is automatically created with a new DCEL instance and is not part of 
    the _faces list. 
  */
  Face _infiniteFace = Face();

  public:

  // accessor for this DCEL's vertices list
  const VertexList & vertices() {
    return _vertices;
  }

  // accessor for this DCEL's half-edges list
  const HalfEdgeList & halfEdges() {
    return _halfEdges;
  }

  // accessor for this DCEL's faces list
  const FaceList & faces() {
    return _faces;
  }

  // accessor for this DCEL's infiniteFace ()
  Face * infiniteFace() {
    return &_infiniteFace;
  }

  /* 
  Function to create a new vertex. It is created and put in the _vertices list. 
  Always use this function to create a new vertex!!

  Returns a pointer to the newly created vertex.
  */
  Vertex * createVertex(double x, double y, double z) {
    _vertices.emplace_back(std::make_unique<Vertex>(x,y,z));
    // return pointer to the vertex we just created
    return _vertices.back().get();
  }

  /* 
  Function to create a new half-edge. It is created and put in the _halfEdges list. 
  Always use this function to create a new half-edge!!

  Returns a pointer to the newly created half-edge.
  */
  HalfEdge * createHalfEdge() {
    // create a new HalfEdge instance in the half-edge list
    _halfEdges.emplace_back(std::make_unique<HalfEdge>());
   
    // return the handle to the half-edge we just created
    return _halfEdges.back().get();
  }

  /* 
  Function to create a new face. It is created and put in the _faces list. 
  Always use this function to create a new face!!

  Returns a pointer to the newly created face.
  */
  Face * createFace() {
    // create a new Face instance in the face list
    _faces.emplace_back(std::make_unique<Face>());
    return _faces.back().get();
  }

  /*
  Function to remove all the DCEL elements that were elminated from the element lists.

  !!! This wil invalidate all the pointers (any pointer you created but also the links of the 
  DCEL elements) to any eliminated elements !!!

  !!! So only do this when you are sure there are no more links to the eliminated elements !!!
  */
  void cleanup() {
    for ( auto it = _vertices.begin(); it != _vertices.end(); ) {
      if ( (*it)->isEliminated() ) {
        it = _vertices.erase(it);
      } else {
        ++it;
      }
    }
    for ( auto it = _halfEdges.begin(); it != _halfEdges.end(); ) {
      if ( (*it)->isEliminated() ) {
        it = _halfEdges.erase(it);
      } else {
        ++it;
      }
    }
    for ( auto it = _faces.begin(); it != _faces.end(); ) {
      if ( (*it)->isEliminated() ) {
        it = _faces.erase(it);
      } else {
        ++it;
      }
    }
  }

  /* 
    Function to check if all elements are valid. See DCELElements.hpp for what valid means.
    It is recommended that you use this function after major changes to the DCEL as to ensure
    that all the links in the DCEL are pointing to a valid Element and that all the elements
    are not eliminated themselves.
    
    Returns 
      - a pointer to first encountered invalid element if at least one element is invalid
      - nullptr if all elements are valid
  */
  DCELElement * findInValid() {
    for ( auto &v : _vertices ) {
      if ( !v->isValid() ) {
        return v.get();
      }
    }
    for ( auto &e : _halfEdges ) {
      if ( !e->isValid() ) {
        return e.get();
      }
    }
    for ( auto &f : _faces ) {
      if ( !f->isValid() ) {
        return f.get();
      }
    }
    // if all of the above passed it means there are no invalid elements, thus we return the nullptr
    return nullptr;
  }
};

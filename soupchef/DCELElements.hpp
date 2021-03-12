#pragma once

#include <memory>

// forward declarations
struct Vertex;
struct HalfEdge;
struct Face;

/* 
  The base class for all DCEL elements. The Vertex, HalfEdge, and Face class inherit from this base.
  All of the methods defined here therefore are also work for the inheriting classes
*/
class DCELElement {

  bool _eliminated = false;
  
  public:
  
  /*
    Helper function for printing the Element to an output stream
  */
  virtual std::ostream& toStream(std::ostream& os) const = 0;

  /*
    By eliminating an element we indicate that this element should no longer be part of the DCEL 
    and that we intent to remove it later (using the cleanup() function of the DCEL class).
  */
  void eliminate() {
    _eliminated = true;
  }

  /*
    Returns whether this element was elminated before:
      
      true  - this element was eliminated
      false - this element was not eliminated
  */
  bool isEliminated() const {
    return _eliminated;
  }

  /*
    This function (which has a separate implementation for the Vertex, HalfEdge and Face structs)
    checks if all the links of the element are pointing to another element that was not elminated
    earlier. After an operation (eg removing an edge, or constructing a DCEL) there should be no 
    dangling links in any element!
  */
  virtual bool hasDanglingLink() const = 0;

  /*
    A DCEL Element is valid if none of its links are dangling and if it was not eliminated
  */
  bool isValid() const {
    return !_eliminated && !hasDanglingLink();
  }

};

/*
  Helper function for printing the Element to an output stream
*/
std::ostream& operator<<(std::ostream& os, const DCELElement& f) {
  f.toStream(os) << " at " << &f;
  
  if ( f.isEliminated() ) {
    os << " [was eliminated]";
  }
  if ( f.hasDanglingLink() ) {
    os << " [has dangling link]";
  }
  return os;
} 

/*
  The Vertex element is the simplest one. It just holds the x,y,z coordinates.
  There are no links to other elements.
*/
struct Vertex : DCELElement {

  double x;
  double y;
  double z;
  
  Vertex(double x, double y, double z) : x(x), y(y), z(z) {}
  Vertex(const Vertex&) = delete;

  bool hasDanglingLink() const override;

  std::ostream& toStream(std::ostream& os) const override {
    os << "Vertex(" << x << ", " << y << ", " << z << ")";
    return os;
  };

};

/*
  The HalfEdge element is the most important one in a DCEL. There are 6 links to other elements:
  2 vertices, 3 other halfedges and 1 face.
*/
struct HalfEdge : DCELElement {
  
  // link to the vertex at the origin of this half-edge
  // this should be the destination of this half-edge's twin!
  Vertex* origin = nullptr;

  // link to the vertex at the origin of this half-edge
  // this should be the origin of this half-edge's twin!
  Vertex* destination = nullptr;
  
  // link to this half-edge's twin
  // The twin should always go in the opposite direction!
  HalfEdge* twin = nullptr;

  // link to this half-edge's next half-edge
  // The next half-edge should always have the same incidident face!
  HalfEdge* next = nullptr;

  // link to this half-edge's previous half-edge
  // The previous half-edge should always have the same incidident face!
  HalfEdge* prev = nullptr;

  // link to the face that is incident to this half-edge
  Face* incidentFace = nullptr;

  HalfEdge() {};
  HalfEdge(const HalfEdge&) = delete;

  bool hasDanglingLink() const override;

  std::ostream& toStream(std::ostream& os) const override {
    os << "HalfEdge";
    return os;
  };

};

/*
  Finally there is the Face element. There are 6 links to other elements:
  2 vertices, 3 other halfedges and 1 face.
*/
struct Face : DCELElement {

  // link to one of the half-edges on the exterior chain of this face
  HalfEdge* exteriorEdge = nullptr;

  // a list containing one half-edge for each interior chain (a hole) of this face
  std::list<HalfEdge*> holes;
  
  Face() {};
  Face(const Face&) = delete;

  bool hasDanglingLink() const override;
  
  std::ostream& toStream(std::ostream& os) const override {
    os << "Face";
    return os;
  };

};

/*
  Following are the implementations of hasDanglingLink() for each DCEL element type.
*/
bool Vertex::hasDanglingLink() const {
    // A vertex does not have any links to other DCEL elements, so it can never have dangling links either
    return false;
  }

bool HalfEdge::hasDanglingLink() const {
  // check if all the members of this HalfEdge are linking to a non-eliminated element
  if ( origin == nullptr ) return true;
  if ( origin->isEliminated() ) return true;
  
  if ( destination == nullptr ) return true;
  if ( destination->isEliminated() ) return true;
  
  if ( twin == nullptr ) return true;
  if ( twin->isEliminated() ) return true;
  
  if ( next == nullptr ) return true;
  if ( next->isEliminated() ) return true;
  
  if ( prev == nullptr ) return true;
  if ( prev->isEliminated() ) return true;
  
  if ( incidentFace == nullptr ) return true;
  if ( incidentFace->isEliminated() ) return true;

  // if we get to this point in the function we know there are no dangling links
  return false;
}

bool Face::hasDanglingLink() const {
  // check if the exteriorEdge of this Face is linking to a non-eliminated element
  if ( exteriorEdge == nullptr ) return true;
  if ( exteriorEdge->isEliminated() ) return true;
  
  // check if all the interior half-edges are linking to a non-eliminated element
  for ( auto he : holes ) {
    if ( he == nullptr ) return true;
    if ( he->isEliminated() ) return true;
  }

  // if we get to this point in the function we know there are no dangling links
  return false;
}

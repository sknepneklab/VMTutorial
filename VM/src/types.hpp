/*!
 * \file types.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 30-Nov-2023
 * \brief Data types for the Mesh 
 */ 

#ifndef __TYPES_HPP__
#define __TYPES_HPP__

//#define VERSION

#include <vector>
#include <iterator>

#include <pybind11/pybind11.h>

#include "vec.hpp"


namespace py = pybind11;
using std::vector;



namespace VMTutorial
{
  
  template<typename Property> class HalfEdge;
  template<typename Property> class Vertex;
  template<typename Property> class Edge;
  template<typename Property> class Face;
  template<typename Property> class Mesh;
  template<typename Property> class VertexCirculator; 
  template<typename Property> class VertexCCirculator; 
  template<typename Property> class FaceCirculator;
  template<typename Property> class FaceCCirculator;

  template<typename Property> using HEHandle     = typename vector<HalfEdge<Property>>::iterator;
  template<typename Property> using VertexHandle = typename vector<Vertex<Property>>::iterator;
  template<typename Property> using EdgeHandle   = typename vector<Edge<Property>>::iterator;
  template<typename Property> using FaceHandle   = typename vector<Face<Property>>::iterator;

  template<typename Property> using HECHandle     = typename vector<HalfEdge<Property>>::const_iterator;
  template<typename Property> using VertexCHandle = typename vector<Vertex<Property>>::const_iterator;
  template<typename Property> using EdgeCHandle   = typename vector<Edge<Property>>::const_iterator;
  template<typename Property> using FaceCHandle   = typename vector<Face<Property>>::const_iterator;
  
  //! HalfEdge class
  template<typename Property>
  class HalfEdge
  {

    public:

      // Constructors
      HalfEdge(Mesh<Property>& mesh) :  _mesh{mesh},
                                        _idx{-1}, 
                                        _property{typename Property::HEProperty()}, 
                                        _from{-1},
                                        _to{-1},
                                        _edge{-1},
                                        _face{-1},
                                        _pair{-1},
                                        _next{-1},
                                        _prev{-1},
                                        erased{false}  
                                        { 

                                        }
      
      // Member functions
      int idx() { return _idx; }  
      void set_idx(int idx) { _idx = idx; }    

      typename Property::HEProperty& data() { return _property; }
      typename Property::HEProperty  data() const { return _property; }

      VertexHandle<Property>  from() {  return _mesh.get_mesh_vertex(_from); }
      VertexHandle<Property>  to()   {  return _mesh.get_mesh_vertex(_to); }

      VertexCHandle<Property> from() const  {  return _mesh.get_mesh_vertex(_from); }
      VertexCHandle<Property> to()   const  {  return _mesh.get_mesh_vertex(_to);  }
 
      EdgeHandle<Property>   edge()  {  return _mesh.get_mesh_edge(_edge);  }
      EdgeCHandle<Property>  edge() const  {  return _mesh.get_mesh_edge(_edge);   }

      FaceHandle<Property>   face()  {  return _mesh.get_mesh_face(_face); }
      FaceCHandle<Property>  face()  const {  return _mesh.get_mesh_face(_face); }

      HEHandle<Property>    pair()  {  return _mesh.get_mesh_he(_pair); }
      HEHandle<Property>    next()  {  return _mesh.get_mesh_he(_next); }
      HEHandle<Property>    prev()  {  return _mesh.get_mesh_he(_prev); }

      HECHandle<Property>     pair() const {  return _mesh.get_mesh_he(_pair); }
      HECHandle<Property>     next() const {  return _mesh.get_mesh_he(_next); }
      HECHandle<Property>     prev() const {  return _mesh.get_mesh_he(_prev); }


      Vec direction() { return this->to()->r - this->from()->r; }

      // Public members
      bool erased;

      friend class Mesh<Property>;
      
    private:

      int _idx;                               // Half-edge index (used for debugging)
      typename Property::HEProperty          _property;

      int  _from;         // vertex it starts from (that vertex will have this halfedge as its he)
      int  _to;           // vertex it points to (its pair will have this vertex as its he)

      int   _edge;          // edge this he is part of

      int    _face;         // face to the left of it, when looking in the direction this he points to

      int   _pair;            // its pair half edge (undefined for boundary edges)
      int   _next;            // next he in the same face
      int   _prev;            // previous he in the same face

      Mesh<Property>& _mesh;
      
  };

  
  //!< Vertex class
  template<typename Property>
  class Vertex
  {
    public:

      // Constructors    
      Vertex(Mesh<Property>& mesh)  : _mesh{mesh},
                                      id{0}, 
                                      r{0.0,0.0}, 
                                      _he{-1},
                                      _property{typename Property::VertexProperty()}, 
                                      erased{false}, 
                                      boundary{false} 
      { 

      }
      Vertex(int id, const Vec& r, Mesh<Property>& mesh) : _mesh{mesh},
                                                           id{id}, 
                                                           r{r}, 
                                                           _property{typename Property::VertexProperty()}, 
                                                           _he{-1},
                                                           erased{false}, 
                                                           boundary{false} 
      {

      }
      Vertex(int id, const Vec& r, bool bnd, Mesh<Property>& mesh) : _mesh{mesh},
                                                                     id{id}, 
                                                                     r{r},
                                                                     _he{-1},
                                                                     _property{typename Property::VertexProperty()}, 
                                                                     erased{false}, 
                                                                     boundary{bnd} 
      {

      }
      
      
      // Member functions
      typename Property::VertexProperty& data() { return _property; }
      typename Property::VertexProperty data() const { return _property; }

      HEHandle<Property>  he()  { return _mesh.get_mesh_he(_he); }
      HECHandle<Property> he() const { return _mesh.get_mesh_he(_he); }

      VertexCirculator<Property> circulator() { return VertexCirculator<Property>(this->he()); }
      VertexCCirculator<Property> circulator() const { return VertexCCirculator<Property>(this->he()); }

      // Public members 
      Vec   r;              // position
      int   id;             // unique id
      bool  erased;         // marks vertices that are not connected to the rest of the mesh, but are still in memory 
      bool  boundary;       // if true, vertex is on boundary 
      int   coordination;   // number of neighbours this vertex has

      friend class Mesh<Property>;

    private:

      typename Property::VertexProperty  _property;
      int   _he;          // outgoing half edge
      Mesh<Property>& _mesh;
         
  };

  
  //!< Edge class
  template<typename Property>
  class Edge
  {
    public:

      // Constructors    
      Edge(Mesh<Property>& mesh)  : _mesh{mesh},
                                    _idx{0}, 
                                    i{0}, 
                                    j{0},
                                    _he{-1}, 
                                    _property{typename Property::EdgeProperty()}, 
                                    boundary{true}, 
                                    erased{false} 
      {

      }
      Edge(int i, int j, Mesh<Property>& mesh) : _mesh{mesh},
                                                 _idx{0}, 
                                                 i{i}, 
                                                 j{j}, 
                                                 _property{typename Property::EdgeProperty()}, 
                                                 boundary{true}, 
                                                 erased{false} 
      { 

      }
      
      // Member functions
      int idx() const { return _idx;  }
      void set_idx(int idx) { _idx = idx;  }
      typename Property::EdgeProperty& data() { return _property; }
      typename Property::EdgeProperty  data() const { return _property; }

      HEHandle<Property>  he()  { return _mesh.get_mesh_he(_he); }
      HECHandle<Property> he() const { return _mesh.get_mesh_he(_he); }

      // Public members
      int i, j;               // indices of two vertices
      bool boundary;          // if true, edge is a boundary edge
      bool erased;            // marks all erased edges that are still in memory

      friend class Mesh<Property>;

    private:

      int _idx;  // Unique edge index
      typename Property::EdgeProperty _property;
      int   _he;     // one of the two half edges
      Mesh<Property>& _mesh;
      
  };

  //!< Face class
  template<typename Property>
  class Face
  {
    public:

      // Constructors    
      Face(Mesh<Property>& mesh)  : _mesh{mesh},
                                    id{0}, 
                                    _he{-1},
                                    _property{typename Property::FaceProperty()}, 
                                    outer{false}, 
                                    erased{false} 
      {

      }
      Face(int id, Mesh<Property>& mesh)  : _mesh{mesh},
                                            id{id}, 
                                            _he{-1},
                                            _property{typename Property::FaceProperty()}, 
                                            outer{false}, 
                                            erased{false} 
      { 

      }
      Face(int id, bool erased, Mesh<Property>& mesh)  : _mesh{mesh},
                                                         id{id}, 
                                                         _property{typename Property::FaceProperty()}, 
                                                         outer{false}, 
                                                         erased{erased} 
      {

      }
      
      // Member functions
      typename Property::FaceProperty& data() { return _property; }
      typename Property::FaceProperty  data() const { return _property; }

      HEHandle<Property>  he()  { return _mesh.get_mesh_he(_he); }
      HECHandle<Property> he() const { return _mesh.get_mesh_he(_he); }

      FaceCirculator<Property> circulator() { return FaceCirculator<Property>(this->he()); }
      FaceCCirculator<Property> circulator() const { return FaceCCirculator<Property>(this->he()); }
    
      // public members 
      int id;        // face id
      bool outer;    // if true, face is a ghost outer face
      int nsides;     // number of sides face has
      bool erased;    // if true, face is marked as erased

      friend class Mesh<Property>;

    private:

      typename Property::FaceProperty      _property;
      int   _he;           // one of its half edges
      Mesh<Property>& _mesh;
      
  };

// Vertex circulator
template<typename Property>
class VertexCirculator 
{

public:

    using iterator_category = std::forward_iterator_tag;
    using value_type = HalfEdge<Property>;
    using difference_type = std::ptrdiff_t;
    using pointer = HalfEdge<Property>*;
    using reference = HalfEdge<Property>&;

    VertexCirculator() : _start{}, 
                         _current{},
                         _isEnd{true}
    {
        
    }

    explicit VertexCirculator(HEHandle<Property> he) : _start{he}, 
                                                       _current{he},
                                                       _isEnd{false}
    {

    }

    VertexCirculator& operator++() 
    {
        _current = _current->next()->pair();
        if (_current == _start) 
        {
            _isEnd = true; // Completed full circle
        }
        return *this;
    }

    VertexCirculator operator++(int) 
    {
        VertexCirculator temp = *this;
        ++(*this);
        return *temp;
    }

    reference operator*() 
    {
        return *_current;
    }

    pointer operator->() 
    {
        return &(*_current);
    }

    bool operator==(const VertexCirculator& other) const 
    {
        return (_isEnd && other._isEnd) || (_current == other._current);
    }

    bool operator!=(const VertexCirculator& other) const 
    {
        return !(*this == other);
    }

    // Begin and end methods for range-based for loop
    VertexCirculator begin()  { return *this; }
    VertexCirculator cbegin() const  { return *this; }
    VertexCirculator end()  { return VertexCirculator(); }
    VertexCirculator cend() const { return VertexCirculator(); }

private:

    HEHandle<Property> _start;
    HEHandle<Property> _current;
    bool _isEnd;

};

// constant Vertex circulator
template<typename Property>
class VertexCCirculator 
{

public:

    using iterator_category = std::forward_iterator_tag;
    using value_type = HalfEdge<Property>;
    using difference_type = std::ptrdiff_t;
    using pointer = const HalfEdge<Property>*;
    using reference = const HalfEdge<Property>&;

    VertexCCirculator() : _start{}, 
                         _current{},
                         _isEnd{true}
    {
        
    }

    explicit VertexCCirculator(HECHandle<Property> he) : _start{he}, 
                                                       _current{he},
                                                       _isEnd{false}
    {

    }

    VertexCCirculator& operator++() 
    {
        _current = _current->next()->pair();
        if (_current == _start) 
        {
            _isEnd = true; // Completed full circle
        }
        return *this;
    }

    VertexCCirculator operator++(int) 
    {
        VertexCCirculator temp = *this;
        ++(*this);
        return *temp;
    }

    const reference operator*() const
    {
        return *_current;
    }

    pointer operator->() 
    {
        return &(*_current);
    }

    bool operator==(const VertexCCirculator& other) const 
    {
        return (_isEnd && other._isEnd) || (_current == other._current);
    }

    bool operator!=(const VertexCCirculator& other) const 
    {
        return !(*this == other);
    }

    // Begin and end methods for range-based for loop
    VertexCCirculator begin()  { return *this; }
    VertexCCirculator cbegin() const  { return *this; }
    VertexCCirculator end()  { return VertexCCirculator(); }
    VertexCCirculator cend() const { return VertexCCirculator(); }

private:

    HECHandle<Property> _start;
    HECHandle<Property> _current;
    bool _isEnd;

};

// Face Circulator
template<typename Property>
class FaceCirculator 
{

public:

    using iterator_category = std::forward_iterator_tag;
    using value_type = HalfEdge<Property>;
    using difference_type = std::ptrdiff_t;
    using pointer = const HalfEdge<Property>*;
    using reference = const HalfEdge<Property>&;

    FaceCirculator() : _start{}, 
                       _current{}
    {
        
    }
    explicit FaceCirculator(HEHandle<Property> he) : _start{he}, 
                                                     _current{he},
                                                     _isEnd{false}
                                 
    {
       
    }

    FaceCirculator& operator++() 
    {
        _current = _current->next();
        if (_current == _start) 
        {
            _isEnd = true; // Completed full circle
        }
        return *this;
    }

    FaceCirculator operator++(int) 
    {
        FaceCirculator temp = *this;
        ++(*this);
        return *temp;
    }

    reference operator*() 
    {
        return *_current;
    }

    pointer operator->() 
    {
        return &(*_current);
    }

    bool operator==(const FaceCirculator& other) const 
    {
        return (_isEnd && other._isEnd) || (_current == other._current);
    }

    bool operator!=(const FaceCirculator& other) const 
    {
        return !(*this == other);
    }

    // Begin and end methods for range-based for loop
    FaceCirculator begin()  { return *this; }
    FaceCirculator cbegin() const { return *this; }
    FaceCirculator end()  { return FaceCirculator(); }
    FaceCirculator cend() const { return FaceCirculator(); }

private:

    HEHandle<Property> _start;
    HEHandle<Property> _current;
    bool _isEnd;

};

// Constant Face Circulator
template<typename Property>
class FaceCCirculator 
{

public:

    using iterator_category = std::forward_iterator_tag;
    using value_type = HalfEdge<Property>;
    using difference_type = std::ptrdiff_t;
    using pointer = const HalfEdge<Property>*;
    using reference = const HalfEdge<Property>&;

    FaceCCirculator() : _start{}, 
                       _current{}
    {
        
    }
    explicit FaceCCirculator(HECHandle<Property> he) : _start{he}, 
                                                     _current{he},
                                                     _isEnd{false}
                                 
    {
       
    }

    FaceCCirculator& operator++() 
    {
        _current = _current->next();
        if (_current == _start) 
        {
            _isEnd = true; // Completed full circle
        }
        return *this;
    }

    FaceCCirculator operator++(int) 
    {
        FaceCCirculator temp = *this;
        ++(*this);
        return *temp;
    }

    const reference operator*() const
    {
        return *_current;
    }

    pointer operator->() 
    {
        return &(*_current);
    }

    bool operator==(const FaceCCirculator& other) const 
    {
        return (_isEnd && other._isEnd) || (_current == other._current);
    }

    bool operator!=(const FaceCCirculator& other) const 
    {
        return !(*this == other);
    }

    // Begin and end methods for range-based for loop
    FaceCCirculator begin()  { return *this; }
    FaceCCirculator cbegin() const { return *this; }
    FaceCCirculator end()  { return FaceCCirculator(); }
    FaceCCirculator cend() const { return FaceCCirculator(); }

private:

    HECHandle<Property> _start;
    HECHandle<Property> _current;
    bool _isEnd;

};

  
}

#endif

#ifndef SCREAM_FIELD_LAYOUT_HPP
#define SCREAM_FIELD_LAYOUT_HPP

#include "field_tag.hpp"
#include <ekat/std_meta/ekat_std_utils.hpp>
#include <ekat/ekat_assert.hpp>

#include <string>
#include <vector>

namespace scream
{

// The type of the layout, that is, the kind of field it represent.
enum class LayoutType {
  Invalid,
  Scalar2D,
  Vector2D,
  Tensor2D,
  Scalar3D,
  Vector3D,
  Tensor3D
};

inline std::string e2str (const LayoutType lt) {
  std::string name;
  switch (lt) {
    case LayoutType::Scalar2D: name = "Scalar2D"; break;
    case LayoutType::Vector2D: name = "Vector2D"; break;
    case LayoutType::Tensor2D: name = "Tensor2D"; break;
    case LayoutType::Scalar3D: name = "Scalar3D"; break;
    case LayoutType::Vector3D: name = "Vector3D"; break;
    case LayoutType::Tensor3D: name = "Tensor3D"; break;
    case LayoutType::Invalid:  name = "INVALID" ; break;
    default:
      EKAT_ERROR_MSG ("Error! Unrecognized LayoutType.\n");
  }
  return name;
}

/*
 *  A small class to hold basic info about a field layout
 *
 */

class FieldLayout {
public:

  // Constructor(s)
  FieldLayout () = delete;
  FieldLayout (const FieldLayout&) = default;
  FieldLayout (const std::initializer_list<FieldTag>& tags);
  FieldLayout (const std::vector<FieldTag>& tags);
  FieldLayout (const std::vector<FieldTag>& tags,
               const std::vector<int>& dims);

  // Assignment (defaulted)
  FieldLayout& operator= (const FieldLayout&) = default;

  // Create invalid layout
  static FieldLayout invalid () { return FieldLayout({}); }

  // ----- Getters ----- //

  // Name and layout informations
  const std::vector<FieldTag>& tags () const { return m_tags; }
  FieldTag tag  (const int idim) const;
  bool has_tag (const FieldTag t) const { return ekat::contains(m_tags,t); }

  // The rank is the number of tags associated to this field.
  int     rank () const  { return m_rank; }

  int dim (const FieldTag tag) const;
  int dim (const int idim) const;
  const std::vector<int>& dims () const { return m_dims; }

  int      size ()               const;

  bool is_dimension_set  (const int idim) const;
  bool are_dimensions_set () const;

  // Check if this layout is that of a vector fielt
  bool is_vector_layout () const;

  // If this is the layout of a vector field, get the idx of the vector dimension
  // Note: throws if is_vector_layout()==false.
  int get_vector_dim () const;

  // ----- Setters ----- //

  // Note: as soon as a dimension is set, it cannot be changed.
  void set_dimension  (const int idim, const int dimension);
  void set_dimensions (const std::vector<int>& dims);

protected:

  int                   m_rank;
  std::vector<FieldTag> m_tags;
  std::vector<int>      m_dims;
};

bool operator== (const FieldLayout& fl1, const FieldLayout& fl2);
LayoutType get_layout_type (const std::vector<FieldTag>& field_tags);
std::string to_string (const FieldLayout& l);

// ========================== IMPLEMENTATION ======================= //

inline int FieldLayout::dim (const FieldTag t) const {
  auto it = ekat::find(m_tags,t);

  // Check if found
  EKAT_REQUIRE_MSG(it!=m_tags.end(), "Error! Tag '" + e2str(t) + "' not found.\n");

  // Check only one tag (no ambiguity)
  EKAT_REQUIRE_MSG(ekat::count(m_tags,t)==1,
                     "Error! Tag '" + e2str(t) + "' appears multiple times.\n"
                     "       You must inspect tags() and dims() manually.\n");

  return m_dims[std::distance(m_tags.begin(),it)];
}

inline int FieldLayout::dim (const int idim) const {
  ekat::error::runtime_check(idim>=0 && idim<m_rank, "Error! Index out of bounds.", -1);
  return m_dims[idim];
}

inline int FieldLayout::size () const {
  ekat::error::runtime_check(are_dimensions_set(), "Error! Field dimensions not yet set.\n",-1);
  int prod = m_rank>0 ? 1 : 0;
  for (int idim=0; idim<m_rank; ++idim) {
    prod *= m_dims[idim];
  }
  return prod;
}

inline FieldTag FieldLayout::tag (const int idim) const { 
  ekat::error::runtime_check(idim>=0 && idim<m_rank, "Error! Index out of bounds.", -1);
  return m_tags[idim];
} 

inline bool FieldLayout::is_dimension_set (const int idim) const {
  ekat::error::runtime_check(idim>=0 && idim<m_rank, "Error! Index out of bounds.", -1);
  return m_dims[idim]>=0;
}

inline bool FieldLayout::are_dimensions_set () const {
  for (int idim=0; idim<m_rank; ++idim) {
    if (m_dims[idim]<0) {
      return false;
    }
  }
  return true;
}

inline bool operator== (const FieldLayout& fl1, const FieldLayout& fl2) {
  return fl1.rank()==fl2.rank() &&
         fl1.tags()==fl2.tags() &&
         fl1.dims()==fl2.dims();
}

} // namespace scream

#endif // SCREAM_FIELD_LAYOUT_HPP


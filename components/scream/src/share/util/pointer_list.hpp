#ifndef SCREAM_POINTER_LIST
#define SCREAM_POINTER_LIST

#include <vector>

namespace scream {

// A PointerList is a linearly traversible list of pointers that provides
// iterators that automatically double-dereference their referents. Under the
// covers, a PointerList is just a vector of the given pointer type.
template <typename PointerType, typename ValueType>
class PointerList final {

public:

  class iterator final {
  public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = ValueType;
    using difference_type = std::ptrdiff_t;
    using pointer = PointerType;
    using reference = ValueType&;
    using base_iter_type = typename std::vector<PointerType>::iterator;

    explicit iterator(base_iter_type iter) : m_iter(iter) {}
    iterator(const iterator& iter) : m_iter(iter.m_iter) {}
    iterator& operator=(const iterator& iter) {
      m_iter = iter.m_iter;
      return *this;
    }

    iterator& operator++() { m_iter++; return *this;}
    iterator operator++(int) { auto retval = *this; m_iter++; return retval; }
    bool operator==(iterator other) const {return m_iter == other.m_iter;}
    bool operator!=(iterator other) const {return m_iter != other.m_iter;}
    pointer operator->() const {return *m_iter;}
    reference operator*() const {return **m_iter;}
  private:
    base_iter_type m_iter;
  };

  class const_iterator final {
  public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = const ValueType;
    using difference_type = std::ptrdiff_t;
    using pointer = const PointerType;
    using reference = const ValueType&;
    using base_iter_type = typename std::vector<PointerType>::const_iterator;

    explicit const_iterator(base_iter_type iter) : m_iter(iter) {}
    const_iterator(const const_iterator& iter) : m_iter(iter.m_iter) {}
    const_iterator& operator=(const const_iterator& iter) {
      m_iter = iter.m_iter;
      return *this;
    }

    const_iterator& operator++() { m_iter++; return *this;}
    const_iterator operator++(int) { auto retval = *this; m_iter++; return retval; }
    bool operator==(const_iterator other) const {return m_iter == other.m_iter;}
    bool operator!=(const_iterator other) const {return m_iter != other.m_iter;}
    pointer operator->() const {return *m_iter;}
    reference operator*() const {return **m_iter;}
  private:
    base_iter_type m_iter;
  };

  // These iterators provide access to the list's contents.
  iterator begin() { return iterator(m_list.begin()); }
  const_iterator begin() const { return const_iterator(m_list.begin()); }
  iterator end() { return iterator(m_list.end()); }
  const_iterator end() const { return const_iterator(m_list.end()); }

  // Returns the number of elements in the list.
  size_t size() const { return m_list.size(); }

  // Adds a pointer to the end of the list.
  void append(PointerType ptr) { m_list.push_back(ptr); }

private:

  std::vector<PointerType> m_list;
};

} // end namespace scream

#endif

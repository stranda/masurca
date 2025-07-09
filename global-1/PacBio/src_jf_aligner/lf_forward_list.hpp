/******************************************
Copyright University of Maryland 2015
******************************************/
#ifndef _LF_FORWARD_LIST_HPP_
#define _LF_FORWARD_LIST_HPP_

#include <memory>

template<typename T, bool LF = true>
struct lf_forward_list_base {
  struct head_node {
    mutable head_node* next_;
  };
  struct node : head_node {
    T val_;
  };

  lf_forward_list_base(head_node* h)  { head_.next_ = h; }

  class iterator;
  // Const iterator
  class const_iterator : public std::iterator<std::forward_iterator_tag, T> {
  public:
    friend class iterator;
    friend struct lf_forward_list_base<T, LF>;
    explicit const_iterator(const head_node* head = 0) : head_(head) { }
    const_iterator(const const_iterator& rhs) : head_(rhs.head_) { }
    const_iterator(const iterator& rhs) : head_(rhs.head_) { }
    const_iterator& operator=(const const_iterator& rhs) {
      head_ = rhs.head_;
      return *this;
    }
    bool operator==(const const_iterator& rhs) const { return head_ == rhs.head_; }
    bool operator!=(const const_iterator& rhs) const { return head_ != rhs.head_; }
    const T& operator*() const { return static_cast<const node*>(head_)->val_; }
    const T* operator->() const { return &static_cast<const node*>(head_)->val_; }
    const_iterator& operator++() { head_ = (*const_cast<head_node* volatile*>(&head_))->next_; return *this; }
    const_iterator operator++(int) {
      const_iterator res(*this);
      ++*this;
      return res;
    }
  private:
    const head_node* head_;
  };

  // Iterator
  class iterator : public std::iterator<std::forward_iterator_tag, T> {
  public:
    friend class const_iterator;
    friend struct lf_forward_list_base<T, LF>;
    explicit iterator(head_node* head = 0) : head_(head) { }
    iterator(const iterator& rhs) : head_(rhs.head_) { }
    iterator& operator=(const iterator& rhs) {
      head_ = rhs.head_;
      return *this;
    }
    bool operator==(const iterator& rhs) const { return head_ == rhs.head_; }
    bool operator!=(const iterator& rhs) const { return head_ != rhs.head_; }
    bool operator==(const const_iterator& rhs) const { return head_ == rhs.head_; }
    bool operator!=(const const_iterator& rhs) const { return head_ != rhs.head_; }
    T& operator*() { return static_cast<node*>(head_)->val_; }
    const T& operator*() const { return static_cast<const node*>(head_)->val_; }
    T* operator->() { return &static_cast<node*>(head_)->val_; }
    const T* operator->() const { return &static_cast<const node*>(head_)->val_; }
    iterator& operator++() { head_ = (*const_cast<head_node* volatile*>(&head_))->next_; return *this; }
    iterator operator++(int) {
      iterator res(*this);
      ++*this;
      return res;
    }
  private:
    head_node* head_;
  };

  /** Get an iterator */
  iterator begin() { return iterator(head_.next_); }
  const_iterator begin() const { return const_iterator(head_.next_); }
  const_iterator cbegin() const { return const_iterator(head_.next_); }
  iterator before_begin() { return iterator(&head_); }
  const_iterator before_begin() const { return const_iterator(&head_); }
  const_iterator cbefore_begin() const { return const_iterator(&head_); }

  iterator end() { return iterator(); }
  const_iterator end() const { return const_iterator(); }
  const_iterator cend() const { return const_iterator(); }

  static void insert_after_(const head_node* head, head_node* n) {
    if(LF) {
      head_node* hn = *const_cast<head_node* volatile*>(&(head->next_));
      do {
        n->next_ = hn;
        hn = __sync_val_compare_and_swap(&head->next_, hn, n);
      } while(hn != n->next_);
    } else {
      n->next_    = head->next_;
      head->next_ = n;
    }
  }

  static void insert_after_(const_iterator pos, head_node* n) {
    insert_after_(pos.head_, n);
  }

  iterator push_front_(head_node* n) {
    insert_after_(&head_, n);
    return iterator(n);
  }

  // Remove elements after position up to last. Somewhat thread
  // safe. It can be mixed with insert operation but not with erase
  // operations on ranges that overlap. The base class does not
  // allocate/deallocate the elements, hence it does not control the
  // life span of the element removed.
  // Returns a pointer to the first node removed from the list (position->next_).
  static head_node* erase_after_(head_node* position, head_node* last) {
    head_node* nn;
    if(LF) {
      head_node* onn;
      nn = *const_cast<head_node* volatile*>(&position->next_);
      do {
        onn = nn;
        nn = __sync_val_compare_and_swap(&position->next_, nn, last);
      } while(nn != onn);
      return nn;
    } else {
      nn              = position->next_;
      position->next_ = last;
      return nn;
    }
  }
  static head_node* erase_after_(iterator position, iterator last) {
    return erase_after_(position.head_, last.head_);
  }

  head_node head_;
};

template<typename T, bool LF = true, class Alloc = std::allocator<T> >
class lf_forward_list : public lf_forward_list_base<T, LF>
{
  typedef lf_forward_list_base<T, LF> super;
  typedef typename super::head_node   head_node;
  typedef typename super::node        node;

public:
  typedef Alloc                          allocator_type;
  typedef T                              value_type;
  typedef typename super::iterator       iterator;
  typedef typename super::const_iterator const_iterator;

  /** Construct an empty list */
  explicit lf_forward_list (const allocator_type& alloc = allocator_type()) :
    super(0), T_alloc_(alloc), node_alloc_(alloc)
  { }

  ~lf_forward_list() {
    clear();
  }

  /** Remove all element from the list. Not safe if other threads are
      accessing elements in the list. */
  void clear() {
    erase_after(this->before_begin(), this->end());
  }

  iterator erase_after(iterator position, iterator last) {
    head_node* hn = super::erase_after_(position, last);
    head_node* ohn;

    while(hn) {
      ohn = hn->next_;
      T_alloc_.destroy(&static_cast<node*>(hn)->val_);
      node_alloc_.deallocate(static_cast<node*>(hn), 1);
      hn = ohn;
    }
    return last;
  }

  /** Push an element to the front of the list (copy). Returns an
      iterator to the newly inserted element. */
  iterator push_front(const value_type& val) {
    node* n = node_alloc_.allocate(1);
    T_alloc_.construct(&n->val_, val);
    return super::push_front_(n);
  }

  /** Push an element to the front of the list (move). Returns an
      iterator to the newly inserted element. */
  iterator push_front(value_type&& val) {
    node* n = node_alloc_.allocate(1);
    T_alloc_.construct(&n->val_, std::move(val));
    return super::push_front_(n);
  }


private:
  typedef typename Alloc::template rebind<T>::other    T_alloc_type;
  typedef typename Alloc::template rebind<node>::other node_alloc_type;
  T_alloc_type    T_alloc_;
  node_alloc_type node_alloc_;
};

#endif /* _LF_FORWARD_LIST_HPP_ */

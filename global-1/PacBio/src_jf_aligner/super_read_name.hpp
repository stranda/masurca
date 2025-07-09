/******************************************
Copyright University of Maryland 2015
******************************************/
#ifndef __SUPER_READ_NAME_H__
#define __SUPER_READ_NAME_H__

#include <limits>
#include <algorithm>
#include <cstdint>
#include <string>
#include <vector>
#include <debug.hpp>

class super_read_name {
public:
  struct u_id_ori {
    union {
      struct {
        uint32_t ori_:1; // 0 -> F, 1 -> R
        uint32_t id_:31;
      } split;
      uint32_t raw;
    };

    u_id_ori() : raw(0) { }
    u_id_ori(uint32_t i, char c) : split{ c == 'R', i} {}
    u_id_ori(uint32_t o, uint32_t n) : split{ o, n} {}

    bool orib() const { return split.ori_; }
    char ori() const { return split.ori_ ? 'R' : 'F'; }
    uint32_t id() const { return split.id_; }
    std::string name() const { return std::to_string(split.id_) + ori(); }
    void reverse() { split.ori_ = ~split.ori_; }
    u_id_ori reversed() const { return u_id_ori((uint32_t)1 - split.ori_, split.id_); }
    bool operator==(const u_id_ori& rhs) const { return raw == rhs.raw;}
    bool operator!=(const u_id_ori& rhs) const { return raw != rhs.raw; }
  };
  typedef std::vector<u_id_ori> unitigs_list;
  static const uint32_t invalid_id = std::numeric_limits<uint32_t>::max() >> 1;
  //  static const u_id_ori invalid_unitig;

private:
  unitigs_list unitigs_;

public:
  super_read_name() = default;
  explicit super_read_name(size_t n) : unitigs_(n) { }
  explicit super_read_name(const std::string& name) : unitigs_(parse(name)) { }
  //  explicit super_read_name(std::vector<long>&& unitigs) : unitigs_(std::move(unitigs)) { }
  // super_read_name(super_read_name&& rhs) : unitigs_(std::move(rhs.unitigs_)) { }
  super_read_name(const super_read_name& rhs) : unitigs_(rhs.unitigs_) { }

  super_read_name& operator=(const std::string& name) {
    unitigs_ = parse(name);
    return *this;
  }
  super_read_name& operator=(const super_read_name& rhs) {
    unitigs_ = rhs.unitigs_;
    return *this;
  }
  bool operator==(const super_read_name& rhs) const { return unitigs_ == rhs.unitigs_; }

  const unitigs_list& unitigs() const { return unitigs_; }
  size_t nb_unitigs() const { return unitigs_.size(); }
  size_t size() const { return nb_unitigs(); }
  std::string name() const;

  void append(const super_read_name& rhs, size_t skip = 0);

  // Prepend the unitigs from rhs, in the range [first, last] (note
  // closed interval on both sides), starting at offset. If offset ==
  // this.size(), it prepends from the end. Returns the next offset
  // where it is free to prepend to.
  size_t prepend(size_t offset, const super_read_name& rhs, size_t first, size_t last);

  size_t prepend(const super_read_name& rhs, size_t first, size_t last) {
    return prepend(size(), rhs, first, last);
  }

  u_id_ori operator[](size_t i) const { return i < nb_unitigs() ? unitigs_[i] : u_id_ori(invalid_id, 'F'); }

  uint32_t unitig_id(size_t i) const { return i < nb_unitigs() ? unitigs_[i].id() : invalid_id; }

  void reverse();

  super_read_name get_reverse() const {
    super_read_name res(*this);
    res.reverse();
    return res;
  }

  // Return the length of the longest overlap by unitigs between two
  // super reads, in a dovetail fashion. The return value is the
  // largest integer m such that the last m k-unitigs of *this are
  // equal to the first m k-unitigs of rhs. Note that this
  // relationship is NOT symmetrical, i.e. if *this overlaps with rhs
  // then it comes before rhs. If the return is 0, there is no
  // overlap.
  int overlap(const super_read_name& rhs) const;

  void print_sequence(std::ostream& os, const std::vector<std::string>& unitigs_sequences, int k_len,
                      size_t start_unitig = 0, ssize_t nb_unitigs = -1) const;

protected:
  static unitigs_list parse(const std::string& name);

  friend std::ostream& operator<<(std::ostream& os, const super_read_name& sr);
};

std::ostream& operator<<(std::ostream& os, const super_read_name& sr);
#endif /* __SUPER_READ_NAME_H__ */

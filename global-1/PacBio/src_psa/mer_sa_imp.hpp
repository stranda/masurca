/******************************************
Copyright University of Maryland 2015
******************************************/
#ifndef __MER_SA_IMP_H__
#define __MER_SA_IMP_H__

#include <cmath>
#include <iterator>
#include <atomic>
#include <algorithm>

#include "barrier.hpp"
#include "slice.hpp"
#include "global_timer.hpp"
#include "said_traits.hpp"

template<typename T>
bool lexicographical_compare_n(T* first1, const size_t len1,
                               T* first2, const size_t len2) {
  return std::lexicographical_compare(first1, first1 + len1, first2, first2 + len2);
}

template<typename T>
bool lexicographical_compare_n(const T* first1, const size_t len1,
                               const T* first2, const size_t len2) {
  return std::lexicographical_compare(first1, first1 + len1, first2, first2 + len2);
}

namespace mer_sa_imp {

inline  uint64_t base_to_code(char c) {
  switch(c) {
  case 'a': case 'A': return 0;
  case 'c': case 'C': return 1;
  case 'g': case 'G': return 2;
  case 't': case 'T': return 3;
  default: return 0;
  }
}

template<typename CHARPTR>
uint64_t str_to_mer(CHARPTR p, unsigned int mer_size) {
  uint64_t m = 0;
  for(unsigned int i = 0; i < mer_size; ++i, ++p)
      m = (m << 2) | base_to_code(*p);
  return m;
}

template<typename U, int len, int l = sizeof(U) * 8 / (2 * len)>
struct cmask {
  static const U v =
    (cmask<U, len, l - 1>::v << (2 * len)) | (((U)1 << len) - 1);
};
template<typename U, int len>
struct cmask<U, len, 0> {
  static const U v = 0;
};

// Fast reverse complement of one word through bit tweedling.
inline uint64_t word_reverse(uint64_t w) {
  typedef uint64_t U;
  w = ((w >> 2)  & cmask<U, 2 >::v) | ((w & cmask<U, 2 >::v) << 2);
  w = ((w >> 4)  & cmask<U, 4 >::v) | ((w & cmask<U, 4 >::v) << 4);
  w = ((w >> 8)  & cmask<U, 8 >::v) | ((w & cmask<U, 8 >::v) << 8);
  w = ((w >> 16) & cmask<U, 16>::v) | ((w & cmask<U, 16>::v) << 16);
  w = ( w >> 32                   ) | ( w                    << 32);
  return w;
}


template<>
inline uint64_t str_to_mer<const_compact_iterator<uint8_t, uint64_t>>(const_compact_iterator<uint8_t, uint64_t> p, unsigned int mer_size) {
  return mer_size ? word_reverse(p.get_bits(2 * mer_size)) >> (sizeof(uint64_t) * 8 - 2 * mer_size) : 0;
}

template<>
inline uint64_t str_to_mer<compact_iterator<uint8_t, uint64_t>>(compact_iterator<uint8_t, uint64_t> p, unsigned int mer_size) {
  return mer_size ? p.get_bits(2 * mer_size) : 0;
}

template<typename CHARPTR>
struct mer_iterator : std::iterator<std::input_iterator_tag, const uint64_t>{
  CHARPTR          m_p;
  mutable uint64_t m_m;
  const uint64_t   m_mask;
  mer_iterator() = default;
  mer_iterator(CHARPTR start, unsigned int mer_size) : m_p(start), m_m(0), m_mask(~(uint64_t)0 >> (8 * sizeof(uint64_t) - 2 * mer_size)) {
    for(unsigned int i = 0; i < mer_size - 1; ++i, ++m_p)
      m_m = (m_m << 2) | base_to_code(*m_p);
    m_m <<= 2;
  }
  // Iterator used as a guard. It is not deferencable
  mer_iterator(CHARPTR pos) : m_p(pos), m_m(0), m_mask(0) { }

  mer_iterator(const mer_iterator& rhs) : m_p(rhs.m_p), m_m(rhs.m_m), m_mask(rhs.m_mask) { }

  bool operator==(const mer_iterator& rhs) const { return m_p == rhs.m_p; }
  bool operator!=(const mer_iterator& rhs) const { return m_p != rhs.m_p; }
  bool operator==(const CHARPTR& rhs) const { return m_p == rhs; }
  bool operator!=(const CHARPTR& rhs) const { return m_p != rhs; }


  const uint64_t operator*() const {
    m_m = m_m | base_to_code(*m_p);
    return m_m & m_mask;
  }

  mer_iterator& operator++() {
    ++m_p;
    m_m <<= 2;
    return *this;
  }

  mer_iterator operator++(int) {
    mer_iterator res(*this);
    ++*this;
    return res;
  }
};

template<>
struct mer_iterator<const_compact_iterator<uint8_t, uint64_t>> : public std::iterator<std::input_iterator_tag, const uint64_t> {
  typedef const_compact_iterator<uint8_t, uint64_t> CHARPTR;
  CHARPTR      m_p;
  unsigned int m_mer_size;

  mer_iterator() = default;
  mer_iterator(CHARPTR start, unsigned int mer_size) : m_p(start), m_mer_size(mer_size) { }
  // Iterator used as a guard. pos is not deferencable
  mer_iterator(CHARPTR pos, unsigned int mer_size, bool guard) : m_p(pos - mer_size + 1), m_mer_size(mer_size) { }

  mer_iterator(const mer_iterator& rhs) : m_p(rhs.m_p), m_mer_size(rhs.m_mer_size) { }

  bool operator==(const mer_iterator& rhs) const { return m_p == rhs.m_p; }
  bool operator!=(const mer_iterator& rhs) const { return m_p != rhs.m_p; }
  bool operator==(const CHARPTR& rhs) const { return m_p + m_mer_size - 1 == rhs; }
  bool operator!=(const CHARPTR& rhs) const { return m_p + m_mer_size - 1 != rhs; }

  const uint64_t operator*() const { return str_to_mer(m_p, m_mer_size); }
  mer_iterator& operator++() { ++m_p; return *this; }

  mer_iterator operator++(int) {
    mer_iterator res(*this);
    ++*this;
    return res;
  }
};

template<>
struct mer_iterator<compact_iterator<uint8_t, uint64_t>> : public std::iterator<std::input_iterator_tag, const uint64_t> {
  typedef compact_iterator<uint8_t, uint64_t> CHARPTR;
  CHARPTR      m_p;
  unsigned int m_mer_size;

  mer_iterator() = default;
  mer_iterator(CHARPTR start, unsigned int mer_size) : m_p(start), m_mer_size(mer_size) { }
  // Iterator used as a guard. pos is not deferencable
  mer_iterator(CHARPTR pos, unsigned int mer_size, bool guard) : m_p(pos - mer_size + 1), m_mer_size(mer_size) { }

  mer_iterator(const mer_iterator& rhs) : m_p(rhs.m_p), m_mer_size(rhs.m_mer_size) { }

  bool operator==(const mer_iterator& rhs) const { return m_p == rhs.m_p; }
  bool operator!=(const mer_iterator& rhs) const { return m_p != rhs.m_p; }
  bool operator==(const CHARPTR& rhs) const { return m_p + m_mer_size - 1 == rhs; }
  bool operator!=(const CHARPTR& rhs) const { return m_p + m_mer_size - 1 != rhs; }

  const uint64_t operator*() const { return str_to_mer(m_p, m_mer_size); }
  mer_iterator& operator++() { ++m_p; return *this; }

  mer_iterator operator++(int) {
    mer_iterator res(*this);
    ++*this;
    return res;
  }
};

template<typename CHARPTR, typename SAIDPTR>
struct SA {
  typedef typename compactsufsort_imp::type_traits<SAIDPTR>::SAIDX          SAIDX;
  typedef typename compactsufsort_imp::const_iterator_traits<SAIDPTR>::type CSAIDPTR;
  static const size_t ALPHABET_SIZE = compactsufsort_imp::alphabet_traits<CHARPTR>::size;

  static SAIDX
  create(CHARPTR T, SAIDPTR SA, SAIDX n, uint64_t* mer_counts, unsigned int mer_size, unsigned int max_size) {
    const size_t nb_counts = (1 << (2 * mer_size)) + 1;
    std::fill_n(mer_counts, nb_counts, (uint64_t)0);

    count_mers(T, n - mer_size + 1, mer_counts, mer_size);
    SAIDX total_mers = partial_sums(mer_counts, nb_counts);
    fill_mers(T, 0, n, SA, mer_counts, mer_size);
    for(size_t i = 0; i < nb_counts - 1; ++i)
      sort_one_mer(T, n, SA, mer_counts[i], mer_counts[i + 1], mer_size, max_size);

    return total_mers;
  }

  static void
  create_thread(int thid, int nb_threads, barrier<std::mutex>* thread_barrier,
                slice_for<SAIDX, barrier<std::mutex>>* slicer,
                CHARPTR T, SAIDPTR SA, SAIDX n, uint64_t* mer_counts,
                unsigned int mer_size, unsigned int max_size) {
    const size_t nb_counts    = (1 << (2 * mer_size)) + 1;
    const size_t counts_mask  = ~(size_t)0 >> (8 * sizeof(size_t) - 2 * mer_size);
    const size_t counts_start = (nb_counts - 1) * thid / nb_threads;
    const size_t counts_end   = (nb_counts - 1) * (thid + 1) / nb_threads;
    const size_t n_step       = std::min((SAIDX)(1024 * 1024), std::max((SAIDX)1, n / nb_threads));

    // Add 1 to length of zeroing if thid == nb_threads because
    // counts_end is at most nb_counts - 1, and need to zero up to
    // nb_counts.
    global_timer.start(*thread_barrier, "zero mer counts");
    std::fill_n(mer_counts + counts_start, counts_end - counts_start + (thid == nb_threads), (uint64_t)0);

    global_timer.start(*thread_barrier, "count mers");
    {
      std::unique_ptr<uint64_t[]> tmp_counts(new uint64_t[nb_counts - 1]);
      std::fill_n(tmp_counts.get(), nb_counts - 1, (uint64_t)0);
      slicer->chunk(0, n, n_step, [=,&tmp_counts](SAIDX s, SAIDX e) {
          e = std::max(s, std::min((SAIDX)(n - mer_size + 1), e));
          count_mers(T + s, e - s, tmp_counts.get(), mer_size);
        });
      for(size_t i = 0, j = counts_start; i < nb_counts - 1; ++i, j = ((j + 1) & counts_mask))
        __sync_fetch_and_add(mer_counts + j, tmp_counts[j]);
    }

    global_timer.start(*thread_barrier, "partial sums");
    if(thread_barrier->wait()) {
#ifndef NDEBUG
      const SAIDX total =
#endif
        partial_sums(mer_counts, nb_counts);
      assert(total == (SAIDX)(n - mer_size + 1));
    }

    typedef typename compactsufsort_imp::parallel_iterator_traits<SAIDPTR>::type PSAIDPTR;
    PSAIDPTR PSA(SA);
    global_timer.start(*thread_barrier, "fill_mers");
    slicer->chunk(0, n, n_step, [=](SAIDX s, SAIDX e) {
        e = std::min(n, (SAIDX)(e + mer_size - 1));
        fill_mers(T, s, e, PSA, mer_counts, mer_size);
      });
#ifndef NDEBUG
    if(thread_barrier->wait())
      assert(mer_counts[nb_counts - 1] == (SAIDX)(n - mer_size + 1));
#endif
    global_timer.start(*thread_barrier, "sorting");
    if(mer_size < max_size) {
      slicer->loop(0, nb_counts - 1, 100, [=](SAIDX i) {
          sort_one_mer(T, n, PSA, mer_counts[i], mer_counts[i + 1], mer_size, max_size);
        });
    }
    global_timer.stop(*thread_barrier);
  }

  SAIDX
  static create_mt(CHARPTR T, SAIDPTR SA, SAIDX n, unsigned int t, uint64_t* mer_counts, unsigned int mer_size, unsigned int max_size) {
    barrier<std::mutex>      SA_barrier(t);
    slice_for<SAIDX, barrier<std::mutex>> slicer(SA_barrier);
    std::vector<std::thread> threads;

    for(unsigned int i = 0; i < t; ++i)
      threads.push_back(std::thread(SA::create_thread, i, t, &SA_barrier, &slicer,
                                    T, SA, n, mer_counts, mer_size, max_size));
    for(unsigned int i = 0; i < t; ++i)
      threads[i].join();
    return 0;
  }


  // Mer counting minimizing cache misses
  static void count_mers(CHARPTR T, SAIDX n, uint64_t* counts, unsigned int mer_size) {

    static const unsigned int cache_size = 16;
    uint64_t mers[cache_size];
    mer_iterator<CHARPTR> mit(T, mer_size);

    if(n >= (SAIDX)cache_size) {
      // Warm up cache
      SAIDX i;
      for(i = 0; i < cache_size; ++mit, ++i) {
        mers[i] = *mit;
        __builtin_prefetch(counts + mers[i], 1, 0);
      }

      // Count with a delay
      for( ; i < n; ++mit, ++i) {
        const int j = i % cache_size;
        ++counts[mers[j]];
        const auto m = mers[j] = *mit;
        __builtin_prefetch(counts + m);
      }

      // Finish up
      for(unsigned int j = 0; j < cache_size; ++j, ++i)
        ++counts[mers[i % cache_size]];
    } else { // n < cache_size
      for(SAIDX i = 0; i < n; ++mit, ++i)
        ++counts[*mit];
    }
  }

  // static void count_mers(CHARPTR T, SAIDX n, uint64_t* counts, unsigned int mer_size) {
  //   const CHARPTR end(T + n);
  //   for(mer_iterator<CHARPTR> mit(T, mer_size); mit != end; ++mit)
  //     ++counts[*mit];
  // }


  static void count_mers_thread(CHARPTR T, SAIDX n, uint64_t* counts, unsigned int mer_size) {
    if(n < mer_size) return;
    const CHARPTR end(T + n);
    for(mer_iterator<CHARPTR> mit(T, mer_size); mit != end; ++mit)
      __sync_fetch_and_add(counts + *mit, (uint64_t)1);
  }

  // Sum up
  static size_t partial_sums(uint64_t* counts, const SAIDX nb_counts) {
    uint64_t tmp1 = counts[0];
    uint64_t tmp2;
    uint64_t sum = 0;

    counts[0] = 0;
    for(SAIDX i = 1; i < nb_counts; ++i) {
      tmp2       = counts[i];
      counts[i]  = sum;
      sum       += tmp1;
      tmp1       = tmp2;
    }
    return sum;
  }

  // Fill up the SA array with the index of the prefixes, partially
  // sorted by their prefix mer.
  template<typename SA_TYPE>
  static void fill_mers(CHARPTR T, SAIDX start, const SAIDX end,
                        SA_TYPE SA,
                        uint64_t* counts, unsigned int mer_size) {
    if(end - start < mer_size) return;
    uint64_t*     offsets = counts + 1;
#ifndef NDEBUG
    const CHARPTR endptr(T + end);
#endif
    
    for(mer_iterator<CHARPTR> mit(T + start, mer_size); start < end - mer_size + 1; ++mit, ++start) {
      assert(mit.m_p < endptr);
      const uint64_t m   = *mit;
      const uint64_t off = __sync_fetch_and_add(offsets + m, (uint64_t)1);
      SA[off] = start;
    }
  }

  template<typename SA_TYPE>
  static void sort_one_mer(CHARPTR T, SAIDX n, SA_TYPE SA,
                           uint64_t start, uint64_t end,
                           unsigned int mer_size, unsigned int max_size) {
    assert(start <= end);
    const SAIDX extra_size = max_size - mer_size;
    auto comp = [=](SAIDX i, SAIDX j) -> bool {
      auto si = std::min(n, i + (SAIDX)mer_size);
      auto ei = std::min(n, si + extra_size);
      auto sj = std::min(n, j + (SAIDX)mer_size);
      auto ej = std::min(n, sj + extra_size);
      return lexicographical_compare_n(T + si, ei - si, T + sj, ej - sj) || (!lexicographical_compare_n(T + sj, ej - sj, T + si, ei - si) && si > sj);
    };
    std::sort(SA + start, SA + end, comp);
  }

  // Assumed that Psize is <= max_size
  template<typename CHARPTR2>
  static std::pair<SAIDX, SAIDX>
  search(const CHARPTR T, SAIDX Tsize, const CSAIDPTR SA, SAIDX SAsize,
         const uint64_t* counts, unsigned int mer_size, unsigned int max_size,
         const CHARPTR2 P, SAIDX Psize) {
    const auto  min_size  = std::min((unsigned int)Psize, mer_size);
    uint64_t    mer       = str_to_mer(P, min_size);
    if(Psize <= mer_size) {
      // XXX: should the 2 here be changed by the number of bits per letter of the alphabet?
      const SAIDX shift   = 2 * (mer_size - Psize);
      uint64_t    lmer    = ~(~mer << shift);
      mer               <<= shift;
      return std::make_pair(counts[lmer + 1] - counts[mer], counts[mer]);
    }
    const SAIDX left  = Psize - mer_size;
    SAIDX       count = counts[mer + 1] - counts[mer];
    if(count == 0) return std::make_pair((SAIDX)0, (SAIDX)0); // Not in SA

    const uint64_t left_mer = str_to_mer(P + min_size, left);
    CSAIDPTR       first    = SA + counts[mer];
    SAIDX          step     = 0;

    // Return mer and lenght at it
    auto get_mer = [=](CSAIDPTR it) {
      const SAIDX pos = *it + mer_size;
      const SAIDX cs  = pos < Tsize ? std::min(left, Tsize - pos) : 0;
      return std::make_pair(str_to_mer(T + pos, cs), cs);
    };
    // Compare mer returned by get_mer with left_mer. Return -1, 0, or
    // 1 if mer is less, equal or greater than left_mer.
    auto compare_mer = [left, left_mer](const std::pair<uint64_t, SAIDX>& mer) {
      if(left == mer.second)
        return mer.first < left_mer ? -1 : (mer.first > left_mer ? 1 : 0);
      // Only possibility is mer.second < left_mer. Then mer and
      // left_mer cannot be equal.
      const uint64_t lmer = left_mer >> (2 * (left - mer.second));
      return mer.first <= lmer ? -1 : 1;
    };

    uint64_t start_mer = get_mer(first).first;
    uint64_t end_mer   = get_mer(first + count - 1).first;

    bool found = false;
    while(count > 0 && !found) {
      step            = std::min(count - 1,
                                 (SAIDX)std::lrint(count * (double)(2 * (left_mer - start_mer) + 1) / (double)(2 * (end_mer - start_mer + 1))));
      auto       it   = first + step;
      const auto cmer = get_mer(it);
      switch(compare_mer(cmer)) {
      case 0:
        found = true;
        break;
      case -1:
        first      = ++it;
        count     -= step + 1;
        start_mer  = cmer.first + 1;
        break;
      case 1:
        count   = step;
        end_mer = cmer.first;
        break;
      }
    }

    if(count == 0) return std::make_pair((SAIDX)0, (SAIDX)0); // Not in SA

    // Find lower range
    CSAIDPTR lower = first;
    {
      SAIDX lstep;
      uint64_t lend_mer = left_mer;
      for(SAIDX ncount = step; ncount > 0; ) {
        lstep           = std::min(ncount - 1,
                                   (SAIDX)std::rint(ncount * (double)(left_mer - start_mer) / (double)(lend_mer - start_mer + 1)));
        auto       it   = lower + lstep;
        const auto cmer = get_mer(it);
        if(compare_mer(cmer) < 0) {
          lower      = ++it;
          ncount    -= lstep + 1;
          start_mer  = cmer.first;
        } else {
          ncount   = lstep;
          lend_mer = cmer.first;
        }
      }
    }

    // Find upper range
    CSAIDPTR upper = first + step + 1;
    {
      SAIDX ustep;
      uint64_t ustart_mer = left_mer;
      for(SAIDX ncount = count - step - 1; ncount > 0; ) {
        ustep           = std::min(ncount - 1,
                                   (SAIDX)std::lrint(ncount * (double)(left_mer + 1 - ustart_mer) / (double)(end_mer - ustart_mer + 1)));
        ustep           = std::min(ncount - 1, ustep);
        auto       it   = upper + ustep;
        const auto cmer = get_mer(it);
        if(compare_mer(cmer) <= 0) {
          upper       = ++it;
          ncount     -= ustep + 1;
          ustart_mer  = cmer.first;
        } else {
          ncount  = ustep;
          end_mer = cmer.first;
        }
      }
    }

    return std::make_pair(upper - lower, lower - SA);
  }

  template<typename T>
  static typename std::enable_if<std::is_signed<T>::value, bool>::type
  is_negative(T x) { return x < 0; }
  template<typename T>
  static typename std::enable_if<!std::is_signed<T>::value, bool>::type
  is_negative(T x) { return false; }

  static bool
  check(const CHARPTR T, SAIDX Tsize, const CSAIDPTR SA, SAIDX SAsize,
        const uint64_t* const counts, unsigned int mer_size, unsigned int max_size,
        std::ostream& out = std::cerr) {
    const size_t nb_counts   = (1 << (2 * mer_size)) + 1;

    if(SAsize > Tsize) {
      out << "SAsize=" << SAsize << '>' << Tsize << "=Tsize" << std::endl;
      return false;
    }
    if(SAsize != counts[nb_counts - 1]) {
      out << "SAsize=" << SAsize << "!=counts[-1]=" << counts[nb_counts - 1] << std::endl;
      return false;
    }
    global_timer.start("Check SA offsets are in range");
    for(SAIDX i = 0; i < SAsize; ++i) {
      if(SA[i] >= Tsize || is_negative(SA[i])) {
        out << "Offset  SA[" << i << "] = " << SA[i] << " is out of range [0, " << Tsize << ")" << std::endl;
        return false;
      }
    }
    global_timer.start("Check mers offsets");
    for(size_t i = 0; i < nb_counts; ++i) {
      if(counts[i] >= Tsize || is_negative(counts[i])) {
        out << "Mer count[" << i << "]=" << counts[i] << " is out of range [0," << Tsize << ")" << std::endl;
        return false;
      }
      if(i < nb_counts - 1) {
        if(counts[i + 1] < counts[i]) {
          out << "Mer counts are not ordered: counts[" << (i + 1) << "]=" << counts[i+1] << '<'
              << counts[i] << "=counts[" << i << ']' << std::endl;
          return false;
        }
        for(size_t j = counts[i]; j < counts[i + 1]; ++j) {
          auto pos = SA[j];
          if(pos + mer_size < Tsize) {
            uint64_t m = str_to_mer(T + pos, mer_size);
            if(m != i) {
              out << "Incorrect mer at index " << j << " pos " << pos << ": " << std::hex << i << "!=" << std::hex << m << std::endl;
              return false;
            }
          }
        }
      }
    }
    global_timer.start("Check partial order");
    for(SAIDX i = 0; i < SAsize - 1; ++i) {
      if(lexicographical_compare_n(T + SA[i + 1], std::min((SAIDX)max_size, Tsize - SA[i+1]),
                                   T + SA[i], std::min((SAIDX)max_size, Tsize - SA[i]))) {
        out << "Not ordered: SA[" << (i + 1) << "]=" << SA[i+1] << "->" << "..." << '<' << "..." << "<-" << SA[i] << "=SA[" << i << ']' << std::endl;
        return false;
      }
    }
    global_timer.stop();
    return true;
  }

  static bool
  check_suffixes(const CHARPTR T, SAIDX Tsize, const CSAIDPTR SA, SAIDX SAsize,
                const uint64_t* const counts,
                unsigned int mer_size, unsigned int max_size,
        std::ostream& out = std::cerr) {
    if(SAsize > Tsize) {
      out << "SAsize=" << SAsize << '>' << Tsize << "=Tsize" << std::endl;
      return false;
    }

    // occurrences contains for each suffix size, the number of time we
    // saw a suffix occurring k times found at position j (with 1<=k<=M
    // and 0<=j<M).
    global_timer.start("Check all suffixes present in SA");
    const unsigned int M = 10;
    size_t occurences[max_size - mer_size + 1][M][M];
    memset(occurences, '\0', sizeof(size_t) * (max_size - mer_size + 1) * M * M);

    for(SAIDX i = 0; i < Tsize - mer_size + 1; ++i) {
      for(unsigned int j = mer_size; j <= std::min((SAIDX)max_size, Tsize - i); ++j) {
        auto res = search(T, Tsize, SA, SAsize, counts, mer_size, max_size,
                          T + i, j);
        if(res.first == 0) {
          out << "Suffix at position " << i << " and length " << j << " not found"
              << std::endl;
          return false;
        }
        SAIDX k = 0;
        for(k = 0; k < res.first; ++k)
          if(SA[res.second + k] == i) break;
        if(k >= res.first) {
          out << "Suffix at position " << i << " and length " << j
              << " has no match at that position" << std::endl;
          return false;
        }
        if(res.first <= M)
          ++occurences[j - mer_size][res.first - 1][k];
      }
    }

    // Check that the suffixes with k occurrences have been seen at all
    // positions (0<=j<k) the same number of times.
    global_timer.start("Check aggregate number of occurences of suffixes");
    for(unsigned int i = 0; i <= max_size - mer_size; ++i) {
      for(unsigned int k = 1; k <= M; ++k) {
        for(unsigned int j = 1; j < k; ++j) {
          if(occurences[i][k - 1][j] != occurences[i][k - 1][0]) {
            // out << "For length " << (i + mer_size) << " occuring " << k << " times, " << j << " not seen the right number of times: "
            //     << occurences[i][k - 1][j] << " != " << occurences[i][k - 1][0] << '\n';
            //            return false;
          }
        }
      }
    }
    global_timer.stop();
    return true;
  }

};
} // namespace mer_sa_imp

#endif /* __MER_SA_IPM_H__ */

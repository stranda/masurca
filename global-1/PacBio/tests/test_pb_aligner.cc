#include <gtest/gtest.h>
#include <src_jf_aligner/jf_aligner.hpp>
#include <src_jf_aligner/coarse_aligner.hpp>
#include <src_jf_aligner/superread_parser.hpp>
#include <tests/misc.hpp>

namespace {
using misc::rev_comp;
using misc::remove_file;

// Create a fake sequence of 1000 bp. Then a PacBio read of length 500
// starting at base 300 of the sequence. The PacBio reads has an error
// a base 100 and an insertion at base 400.
//
// A bunch of super-reads are created, all of length 100. R1 and R2
// spans the left boundary of the PacBio read, R3 and R4 are over the
// sequencing error at base 100, R5 and R6 are over the insertion of
// length 4 at base 400, R7 and R8 span the right boundary of the
// PacBio read, finally R9 and R10 are in the error free area (between
// base 101 and 399 of the PacBio read).
//
// An odd numbered read is forward, and even numbered read is reversed.
//
// This is the picture. X representing the sequencing error on I the
// insertion point.
//
//    Full sequence
// 0     300          400                         700                     1000
// -----------------------------------------------------------------------
// PacBio --------------X---------------------------I------------
//     |----> R1       |----> R3  |----> R9      |----> R5    |----> R7
//      <----| R2    <----| R4      <----| R10    <----| R6    <----| R8

std::string generate_sequences(const char* superreads) {
  std::string sequence(1000, 'A');
  static const char base[4] = {'A', 'C', 'G', 'T' };
  for(size_t i = 0; i < sequence.size(); ++i)
    sequence[i] = base[random() % 4];

  std::string pacbio_sequence = sequence.substr(300, 500);
  char error = pacbio_sequence[100];
  do {
    pacbio_sequence[100] = base[random() % 4];
  } while(error == pacbio_sequence[100]);
  pacbio_sequence[399] = 'C';
  pacbio_sequence[400] = 'G';
  pacbio_sequence.insert(400, "ACGT");

  {
    std::ofstream super(superreads);
    if(!super.good())
      throw std::runtime_error("Can't open superreads temp file");
    super << ">R1\n" << sequence.substr(250, 100) << "\n"
          << ">R2\n" << rev_comp(sequence.substr(275, 100)) << "\n"
          << ">R3\n" << sequence.substr(390, 100) << "\n"
          << ">R4\n" << rev_comp(sequence.substr(380, 100)) << "\n"
          << ">R5\n" << sequence.substr(660, 100) << "\n"
          << ">R6\n" << rev_comp(sequence.substr(670, 100)) << "\n"
          << ">R7\n" << sequence.substr(770, 100) << "\n"
          << ">R8\n" << rev_comp(sequence.substr(780, 100)) << "\n"
          << ">R9\n" << sequence.substr(500, 100) << "\n"
          << ">R10\n" << rev_comp(sequence.substr(550, 100)) << "\n";
  }

  return pacbio_sequence;
}

TEST(PbAligner, FakeSequences) {
  mer_dna::k(15);
  remove_file sr_file;
  std::string pacbio_sequence = generate_sequences(sr_file.path);

  auto psa = superread_parse(sr_file.path);
  EXPECT_EQ((size_t)10, psa.nb_sequences());

  parse_sequence           parser(pacbio_sequence);
  align_pb::frags_pos_type frags_pos;
  // lis_align::affine_capped accept_mer(10, 2, 1e9);
  // lis_align::linear        accept_sequence(10);
  align_pb::fetch_super_reads(psa, parser, frags_pos);
  align_pb::do_all_LIS(frags_pos, lis_align::accept_all(), lis_align::accept_all(), 1);

  EXPECT_EQ((size_t)10, frags_pos.size());
  for(auto it = frags_pos.cbegin(); it != frags_pos.cend(); ++it) {
    const align_pb::mer_lists& ml = it->second;
    int read_id = std::atoi(it->first + 1);
    ASSERT_TRUE(read_id >= 1 && read_id <= 10); // Read id is valid
    SCOPED_TRACE(::testing::Message() << "Read:" << read_id);
    EXPECT_TRUE(std::is_sorted(ml.fwd.offsets.cbegin(), ml.fwd.offsets.cend())); // mers offsets must be sorted
    EXPECT_TRUE(std::is_sorted(ml.bwd.offsets.cbegin(), ml.bwd.offsets.cend())); // mers offsets must be sorted
    EXPECT_EQ(ml.fwd.offsets.size(), ml.fwd.lis.size()); // lis has same size
    {
      auto iit = ml.fwd.lis.cbegin();
      for(auto it = ml.fwd.offsets.cbegin(); it != ml.fwd.offsets.cend(); ++it, ++iit) {
        EXPECT_EQ(it->second, ml.fwd.offsets[*iit].second); // and is equivalent
        EXPECT_EQ(it->first, ml.fwd.offsets[*iit].first); // and is equivalent
      }
    }

    EXPECT_EQ(ml.bwd.offsets.size(), ml.bwd.lis.size()); // lis has same size
    {
      auto iit = ml.bwd.lis.cbegin();
      for(auto it = ml.bwd.offsets.cbegin(); it != ml.bwd.offsets.cend(); ++it, ++iit) {
        EXPECT_EQ(it->second, ml.bwd.offsets[*iit].second); // and is equivalent
        EXPECT_EQ(it->first, ml.bwd.offsets[*iit].first); // and is equivalent
      }
    }

    switch(read_id) {
    case 1:
      ASSERT_EQ((size_t)36, ml.fwd.offsets.size());
      EXPECT_EQ(51, ml.fwd.offsets.front().second);
      EXPECT_EQ(86, ml.fwd.offsets.back().second);
      break;

    case 2:
      ASSERT_EQ((size_t)61, ml.bwd.offsets.size());
      EXPECT_EQ(-61, ml.bwd.offsets.front().second);
      EXPECT_EQ(-1, ml.bwd.offsets.back().second);
      break;

    case 3:
      ASSERT_EQ((size_t)75, ml.fwd.offsets.size());
      EXPECT_EQ(12, ml.fwd.offsets.front().second);
      EXPECT_EQ(86, ml.fwd.offsets.back().second);
      break;

    case 4:
      ASSERT_EQ((size_t)71, ml.bwd.offsets.size());
      EXPECT_EQ(-86, ml.bwd.offsets.front().second);
      EXPECT_EQ(-1, ml.bwd.offsets.back().second);
      break;

    case 5:
      ASSERT_EQ((size_t)70, ml.fwd.offsets.size());
      EXPECT_EQ(1, ml.fwd.offsets.front().second);
      EXPECT_EQ(86, ml.fwd.offsets.back().second);
      break;

    case 6:
      ASSERT_EQ((size_t)70, ml.bwd.offsets.size());
      EXPECT_EQ(-86, ml.bwd.offsets.front().second);
      EXPECT_EQ(-1, ml.bwd.offsets.back().second);
      break;

    case 7:
      ASSERT_EQ((size_t)16, ml.fwd.offsets.size());
      EXPECT_EQ(1, ml.fwd.offsets.front().second);
      EXPECT_EQ(16, ml.fwd.offsets.back().second);
      break;

    case 8:
      ASSERT_EQ((size_t)6, ml.bwd.offsets.size());
      EXPECT_EQ(-86, ml.bwd.offsets.front().second);
      EXPECT_EQ(-81, ml.bwd.offsets.back().second);
      break;

    case 9:
      ASSERT_EQ((size_t)86, ml.fwd.offsets.size());
      EXPECT_EQ(1, ml.fwd.offsets.front().second);
      EXPECT_EQ(86, ml.fwd.offsets.back().second);
      break;

    case 10:
      ASSERT_EQ((size_t)86, ml.bwd.offsets.size());
      EXPECT_EQ(-86, ml.bwd.offsets.front().second);
      EXPECT_EQ(-1, ml.bwd.offsets.back().second);
      break;
    }
  }
}

TEST(NTH_MIN, Random) {
  static const int n = 10;
  static const int N = 1000;
  std::vector<int> nbs;
  for(int i = 0; i < N; ++i)
    nbs.push_back(i);
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::shuffle(nbs.begin(), nbs.end(), std::default_random_engine(seed));

  align_pb::nth_min<int> l(n);
  for(int i : nbs)
    l.add(i);
  EXPECT_EQ(n - 1, l.highest());
} // NTH_MIN.Random

} // namespace {

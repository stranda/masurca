/******************************************
Copyright University of Maryland 2015
******************************************/
#include <memory>
#include <ios>
#include <thread>

#include <src_jf_aligner/superread_parser.hpp>
#include <src_jf_aligner/coarse_aligner.hpp>
#include <src_jf_aligner/fine_aligner.hpp>
#include <src_jf_aligner/overlap_graph.hpp>
#include <src_jf_aligner/misc.hpp>
#include <src_jf_aligner/create_mega_reads_cmdline.hpp>

#include <src_psa/global_timer.hpp>

using align_pb::coarse_aligner;
using align_pb::fine_aligner;
using align_pb::coords_info_type;

typedef create_mega_reads_cmdline cmdline_args;
cmdline_args args;


void create_mega_reads(read_parser* reads, Multiplexer* output_m,
                       const coarse_aligner* align_data, const fine_aligner* short_align_data,
                       overlap_graph* graph_walker, const std::vector<std::string>* unitigs_sequences,
                       Multiplexer* dot_m) {
  parse_sequence                        parser;
  coarse_aligner::thread                aligner(*align_data);
  std::unique_ptr<fine_aligner::thread> short_aligner;
  std::unique_ptr<short_parse_sequence> short_parser;

  overlap_graph::thread                 graph(*graph_walker);
  std::string                           name;
  Multiplexer::ostream                  output(output_m);
  std::unique_ptr<Multiplexer::ostream> dot(dot_m ? new Multiplexer::ostream(dot_m) : 0);
  std::vector<int> sort_array;

  if(short_align_data) {
    short_aligner.reset(new fine_aligner::thread(*short_align_data));
    short_parser.reset(new short_parse_sequence);
  }

  if(dot)
    graph.dot(dot.get());
  switch(args.trim_arg) {
  case cmdline_args::trim::match: graph.trim_match(); break;
  }

  while(true) {
    read_parser::job job(*reads);
    if(job.is_empty()) break;

    for(size_t i = 0; i < job->nb_filled; ++i) { // Process each PB read
      auto name_end = job->data[i].header.find_first_of(" \t\n\v\f\r");
      name = job->data[i].header.substr(0, name_end);
      const size_t pb_size = job->data[i].seq.size();
      parser.reset(job->data[i].seq);
      aligner.align_sequence_max(parser, pb_size);

      const coords_info_type* coords = &aligner.coords();

      if(short_aligner) { // Lets recompute the alignment with smaller mers
        short_parser->reset(job->data[i].seq);
        short_aligner->align_sequence(*short_parser, pb_size, *coords);
        coords = &short_aligner->coords();
      }
      const int n = coords->size();
      if((int)sort_array.size() < n)
        sort_array.resize(n);
      for(int i = 0; i < n; ++i)
        sort_array[i] = i;
      std::sort(sort_array.begin(), sort_array.begin() + n, [coords] (int i, int j) { return (*coords)[i] < (*coords)[j]; });
      coords_info_type sorted_coords;
      for(int i = 0; i < n; ++i)
        sorted_coords.push_back((*coords)[sort_array[i]]);

      graph.reset(sorted_coords, name);
      graph.traverse();
      graph.term_node_per_comp(pb_size, args.density_arg, args.min_length_arg);
      switch(args.tiling_arg) {
      case cmdline_args::tiling::maximal: graph.tile_maximal(); break;
      case cmdline_args::tiling::greedy: graph.tile_greedy(); break;
      case cmdline_args::tiling::weighted: graph.tile_weighted(); break;
      }
      graph.print_mega_reads(output, name, unitigs_sequences);

      output.end_record();
      if(dot) dot->end_record();
    }
 }
}

int main(int argc, char *argv[])
{
  args.parse(argc, argv);
  mer_dna::k(args.mer_arg);
  std::ios::sync_with_stdio(false);

  // Open output file for early error reporting
  output_file output;
  if(args.output_given) {
    output.open(args.output_arg, args.threads_arg);
  } else {
    output.set(std::cout, args.threads_arg);
  }
  output_file dot;
  if(args.dot_given)
    dot.open(args.dot_arg, args.threads_arg);

  // Read k-unitig lengths
  std::vector<int> unitigs_lengths;
  std::vector<std::string> sequences;
  if(args.unitigs_lengths_given) { // File with lengths
    std::ifstream is(args.unitigs_lengths_arg);
    if(!is.good())
      cmdline_args::error() << "Failed to open unitig lengths map file '" << args.unitigs_lengths_arg << "'";
    read_unitigs_lengths(is, unitigs_lengths);
  } else { // Sequence in fasta file given
    std::ifstream is(args.unitigs_sequences_arg);
    if(!is.good())
      cmdline_args::error() << "Failed to open unitigs sequence file '" << args.unitigs_sequences_arg << "'";
    read_unitigs_sequences(is, unitigs_lengths, sequences);
  }

  // Read the super reads
  if(args.fine_mer_given)
    short_mer_type::k(args.fine_mer_arg);
  global_timer.start("Super read parse");
  auto psa = superread_parse(args.superreads_arg.cbegin(), args.superreads_arg.cend(),
                             std::min(short_mer_type::k(), args.psa_min_arg), mer_dna::k());
  global_timer.stop();

  // Prepare I/O
  stream_manager streams(args.pacbio_arg.cbegin(), args.pacbio_arg.cend());
  read_parser    reads(4 * args.threads_arg, 100, 1, streams);

  // Create aligners
  coarse_aligner align_data(psa, mer_dna::k(),
                            args.stretch_factor_arg, args.stretch_constant_arg, args.stretch_cap_arg, args.window_size_arg,
                            true /* forward */, args.max_match_flag,
                            args.max_count_arg ? args.max_count_arg : std::numeric_limits<int>::max(),
                            args.mers_matching_arg / 100.0, args.bases_matching_arg / 100.0);
  align_data.unitigs_lengths(&unitigs_lengths, args.k_mer_arg);
  std::unique_ptr<fine_aligner> short_align_data;
  if(args.fine_mer_given)
    short_align_data.reset(new fine_aligner(psa, args.fine_mer_arg, &unitigs_lengths, args.k_mer_arg));
  else
    short_mer_type::k(mer_dna::k());
  
  // Output candidate mega_reads
  //  std::cerr << args.density_arg << ' ' << args.min_length_arg << '\n';
  global_timer.start("create mega reads");
  overlap_graph graph_walker(args.overlap_play_arg, args.k_mer_arg, unitigs_lengths, args.errors_arg, args.bases_flag);
  std::vector<std::thread> threads;
  for(unsigned int i = 0; i < args.threads_arg; ++i)
    threads.push_back(std::thread(create_mega_reads, &reads, output.multiplexer(),
                                  &align_data, short_align_data.get(),
                                  &graph_walker, args.unitigs_sequences_given ? &sequences : nullptr,
                                  dot.multiplexer()));
  for(auto& th : threads)
    th.join();
  global_timer.stop();

  return 0;
}

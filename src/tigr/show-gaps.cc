//-----------------------------------------------------------------------------
//   Programmer: Adam M Phillippy, The Institute for Genomic Research
//         File: show-snps.cc
//         Date: 12 / 08 / 2004
//
//        Usage: show-snps [options] <deltafile>
//               Try 'show-snps -h' for more information
//
//  Description: For use in conjunction with the MUMmer package. "show-snps"
//              displays human readable (and machine parse-able) single
//             base-pair polymorphisms, including indels from the .delta output
//            of the "nucmer" program. Outputs SNP positions and relative
//          information to stdout.
//
//-----------------------------------------------------------------------------

#include "delta.hh"
#include "tigrinc.hh"
#include <string>
#include <cstdlib>
using namespace std;


//================================================================ Options ====
string  OPT_AlignName;                  // delta file name


//=========================================================== Declarations ====
struct EdgeletQCmp_t
//!< Compares query low coord
{
  bool operator( ) (const DeltaEdgelet_t * i, const DeltaEdgelet_t * j) const
  { return ( i->loQ < j->loQ ); }
};


struct EdgeletRCmp_t
//!< Compares reference low coord
{
  bool operator( ) (const DeltaEdgelet_t * i, const DeltaEdgelet_t * j) const
  { return ( i->loR < j->loR ); }
};


void PrintShuffle(DeltaGraph_t & graph);
void ParseArgs(int argc, char ** argv);
void PrintHelp(const char * s);
void PrintUsage(const char * s);



//============================================================ Definitions ====
//------------------------------------------------------------------- main ----
int main(int argc, char **argv)
{
  DeltaGraph_t graph;

  ParseArgs(argc, argv);

  graph.build(OPT_AlignName, false);

  //-- Keep the union of FlagRLIS and FlagQLIS
  graph.flagWGA();
  graph.clean();

  PrintShuffle(graph);

  return EXIT_SUCCESS;
}


//----------------------------------------------------------- PrintShuffle ----
void PrintShuffle(DeltaGraph_t & graph)
{
  long int nAligns, i;

  DeltaEdgelet_t *A, *PA;
  vector<DeltaEdgelet_t *> aligns;
  vector<long int> qindex;

  map<string, DeltaNode_t>::const_iterator mi;
  vector<DeltaEdge_t *>::const_iterator ei;
  vector<DeltaEdgelet_t *>::iterator eli;

  //-- For each reference sequence
  for ( mi = graph.refnodes.begin( ); mi != graph.refnodes.end( ); ++ mi )
    {
      printf(">%s\n", mi->first.c_str());

      //-- Collect all alignments for this reference sequence
      aligns.clear();
      for ( ei  = (mi->second).edges.begin();
            ei != (mi->second).edges.end(); ++ei )
        for ( eli  = (*ei)->edgelets.begin();
              eli != (*ei)->edgelets.end(); ++eli )
          aligns.push_back (*eli);

      nAligns = aligns.size();
      if ( !nAligns )
        {
          // do something clever
          printf("EMPTY!\n");
          continue;
        }

      //-- Override *stpc* with query ordering
      sort(aligns.begin(), aligns.end(), EdgeletQCmp_t());
      for ( i = 0; i != nAligns; ++i )
        aligns[i]->stpc = i;

      //-- Sort by reference order
      sort(aligns.begin(), aligns.end(), EdgeletRCmp_t());


      //-- Walk alignments, low to high in reference
      PA = NULL;
      long int size = aligns.size();
      for ( long int i = 0; i != size; ++i )
        {
          //-- Only looking at alignments in R's LIS
          if ( !aligns[i]->isRLIS )
            continue;

          A = aligns[i];

          //-- For beginning of sequence to first alignment
          if ( !PA )
            {
              PA = A;
              continue;
            }

          //-- Jump if query coords out of expected order
          if ( (PA->dirR == PA->dirQ && A->stpc - 1 != PA->stpc) ||
               (PA->dirR != PA->dirQ && A->stpc + 1 != PA->stpc) )
            {
              if ( PA->dirR == PA->dirQ )
                printf("JUMP %ld,%ld", PA->hiR, PA->hiQ);
              else
                printf("JUMP %ld,%ld", PA->hiR, PA->loQ);

              if ( A->dirR == A->dirQ )
                printf(" to %ld,%ld\n", A->loR, A->loQ);
              else
                printf(" to %ld,%ld\n", A->loR, A->hiQ);
            }

          PA = A;
        }

      //-- For last alignment to end of sequence

    }
}


//-------------------------------------------------------------- ParseArgs ----
void ParseArgs (int argc, char ** argv)
{
  int ch, errflg = 0;
  optarg = NULL;

  while ( !errflg  &&
          ((ch = getopt (argc, argv, "h")) != EOF) )
    switch (ch)
      {
      case 'h':
        PrintHelp (argv[0]);
        exit (EXIT_SUCCESS);
        break;

      default:
        errflg ++;
      }

  if ( errflg > 0  ||  optind != argc - 1 )
    {
      PrintUsage (argv[0]);
      cerr << "Try '" << argv[0] << " -h' for more information.\n";
      exit (EXIT_FAILURE);
    }

  OPT_AlignName = argv [optind ++];
}


//-------------------------------------------------------------- PrintHelp ----
void PrintHelp (const char * s)
{
  PrintUsage (s);
  cerr
    << "-h            Display help information\n"
    << endl;

  cerr
    << "  Description\n"
    << endl;

  return;
}


//------------------------------------------------------------- PrintUsage ----
void PrintUsage (const char * s)
{
  cerr
    << "\nUSAGE: " << s << "  [options]  <deltafile>\n\n";
  return;
}

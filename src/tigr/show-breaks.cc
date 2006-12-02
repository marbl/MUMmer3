//-----------------------------------------------------------------------------
//   Programmer: Adam M Phillippy, University of Maryland
//         File: show-breaks.cc
//         Date: 12 / 01 / 2006
//
//        Usage: show-breaks [options] <deltafile>
//               Try 'show-breaks -h' for more information
//
//  Description:
//
//-----------------------------------------------------------------------------

#include "delta.hh"
#include "tigrinc.hh"
#include <string>
#include <cstdlib>
#include <cassert>
#include <climits>
using namespace std;


//================================================================ Options ====
string  OPT_AlignName;              // delta file name


//=========================================================== Declarations ====
struct EdgeletLoQCmp_t
//!< Sorts query low coord, lo to hi
{
  bool operator( ) (const DeltaEdgelet_t * i, const DeltaEdgelet_t * j) const
  { return ( i->loQ < j->loQ ); }
};

struct EdgeletLoRCmp_t
//!< Sorts reference lo coord, lo to hi
{
  bool operator( ) (const DeltaEdgelet_t * i, const DeltaEdgelet_t * j) const
  { return ( i->loR < j->loR ); }
};


void PrintBreaks(DeltaGraph_t & graph);
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

  PrintBreaks(graph);

  return EXIT_SUCCESS;
}


//------------------------------------------------------------- PrintBreaks ----
void PrintBreaks(DeltaGraph_t & graph)
{
  long i,j;
  long nAligns, gapR, gapQ, diff;
  DeltaEdgelet_t lpad, rpad;
  lpad.isRLIS = rpad.isRLIS = true;
  lpad.isQLIS = rpad.isQLIS = true;
  lpad.loR = lpad.hiR = lpad.loQ = lpad.hiQ = 0;

  DeltaEdgelet_t *A, *PA, *PGA;
  vector<DeltaEdgelet_t *> aligns;

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
          aligns.push_back(*eli);

      rpad.loR = rpad.hiR = mi->second.len + 1;
      rpad.loQ = rpad.hiQ = LONG_MAX;
      aligns.push_back(&lpad);
      aligns.push_back(&rpad);

      nAligns = aligns.size();

      //-- OVERRIDE *stpc* with loQ QLIS ordering
      sort(aligns.begin(), aligns.end(), EdgeletLoQCmp_t());
      for ( i = 0, j = 0; i != nAligns; ++i )
        aligns[i]->stpc = aligns[i]->isQLIS ? j++ : -1;

      //-- Sort by reference order
      sort(aligns.begin(), aligns.end(), EdgeletLoRCmp_t());
      assert ( aligns[0] == &lpad && aligns[nAligns-1] == &rpad );

      //-- Walk reference coverr alignments, low to high
      PA = PGA = aligns[0];
      for ( i = 1; i != nAligns; ++i )
        {
          if ( !aligns[i]->isRLIS ) continue;

          A = aligns[i];

          gapR = A->loR - PA->hiR - 1;
          if ( gapR > 0 )
            printf("INS %ld\t%ld\t%ld\n", PA->hiR + 1, A->loR - 1, gapR);
     
          if ( A->isQLIS )
            {
              gapR = A->loR - PGA->hiR - 1;

              if ( A->stpc != PGA->stpc + PGA->slope() ||
                   A->slope() != PGA->slope() )
                {
                  if ( A->slope() == PGA->slope() )
                    {
                      //-- Relocation            
                      printf("BRK %ld\t%ld\t%ld\n", PGA->hiR, A->loR, gapR);
                    }
                  else
                    {
                      //-- Inversion
                      printf("IBR %ld\t%ld\t%ld\n", PGA->hiR, A->loR, gapR);
                    }
                }
              else
                {
                  //-- Unknown
                  printf("???? %ld\t%ld\t%ld\n", PGA->hiR, A->loR, gapR);
                }
            }
          else
            {
              //-- Duplication
              printf("DUP %ld\t%ld\t%ld\n",
                     A->loR, A->hiR, A->hiR - A->loR + 1);
            }

          PA = A;
          if ( A->isQLIS ) PGA = A;
        }
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

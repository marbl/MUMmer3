//------------------------------------------------------------------------------
//   Programmer: Adam M Phillippy, The Institute for Genomic Research
//         File: delta-filter.cc
//         Date: 09 / 22 / 2004
//
//        Usage: delta-filter  [options]  <deltafile>
//               Try 'show-coords -h' for more information
//
//  Description: For use in conjunction with the MUMmer package.
//              "delta-filter" cleans delta alignment files by filtering
//             alignments that fail to meet the specifications required by the
//            command line switches. Produces a new, filtered delta file on
//           stdout, and works for both nucleotide and amino-acid alignments.
//
//------------------------------------------------------------------------------

#include "delta.hh"
#include "tigrinc.hh"
#include <map>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>
using namespace std;




//=============================================================== Options ====//
string         OPT_AlignName;                // alignment name parameter

bool           OPT_QLIS         = false;     // do query based LIS
bool           OPT_RLIS         = false;     // do reference based LIS
bool           OPT_GLIS         = false;     // do global LIS
long int       OPT_MinLength    = 0;         // minimum alignment length
float          OPT_MinIdentity  = 0.0;       // minimum %identity
float          OPT_MinUnique    = 0.0;       // minimum %unique
float          OPT_MaxOverlap   = 100.0;     // maximum olap as % of align len
float          OPT_Epsilon      = -1.0;      // negligible alignment score




//========================================================== Fuction Decs ====//

//------------------------------------------------------------- ParseArgs ----//
void ParseArgs (int argc, char ** argv);


//------------------------------------------------------------- PrintHelp ----//
void PrintHelp (const char * s);


//------------------------------------------------------------ PrintUsage ----//
void PrintUsage (const char * s);




//========================================================= Function Defs ====//
int main (int argc, char ** argv)
{
  DeltaGraph_t graph;
  srand (1);

  //-- Command line parsing
  ParseArgs (argc, argv);

  //-- Build the alignment graph from the delta file
  graph . build (OPT_AlignName, true);


  //-- Identity requirements
  if ( OPT_MinIdentity > 0  ||  OPT_MinLength > 0 )
    graph . flagScore (OPT_MinLength, OPT_MinIdentity);

  //-- Uniqueness requirements
  if ( OPT_MinUnique > 0 )
    graph . flagUNIQ (OPT_MinUnique);

  //-- Query-based LIS
  if ( OPT_QLIS )
    graph . flagQLIS (OPT_Epsilon, OPT_MaxOverlap);

  //-- Reference-based LIS
  if ( OPT_RLIS )
    graph . flagRLIS (OPT_Epsilon, OPT_MaxOverlap);

  //-- Global LIS
  if ( OPT_GLIS )
    graph . flagGLIS (OPT_Epsilon);


  //-- Output the filtered delta file
  graph . outputDelta (cout);

  return EXIT_SUCCESS;
}




//------------------------------------------------------------- ParseArgs ----//
void ParseArgs (int argc, char ** argv)
{
  int ch, errflg = 0;
  optarg = NULL;
  
  while ( !errflg  &&
          ((ch = getopt (argc, argv, "e:ghi:l:o:qru:")) != EOF) )
    switch (ch)
      {
      case 'e':
        OPT_Epsilon = atof (optarg);
        break;

      case 'g':
        OPT_GLIS = true;
        break;

      case 'h':
        PrintHelp (argv[0]);
        exit (EXIT_SUCCESS);
        break;

      case 'i':
        OPT_MinIdentity = atof (optarg);
        break;

      case 'l':
        OPT_MinLength = atol (optarg);
        break;

      case 'o':
	OPT_MaxOverlap = atof (optarg);
	break;

      case 'q':
        OPT_QLIS = true;
        break;

      case 'r':
        OPT_RLIS = true;
        break;

      case 'u':
        OPT_MinUnique = atof (optarg);
        break;

      default:
        errflg ++;
      }

  if ( OPT_MinIdentity < 0.0  ||  OPT_MinIdentity > 100.0 )
    {
      cerr << "ERROR: Minimum identity must be within the range [0, 100]\n";
      errflg ++;
    }

  if ( OPT_MinLength < 0 )
    {
      cerr << "ERROR: Minimum length must be greater than or equal to zero\n";
      errflg ++;
    }

  if ( OPT_MinUnique < 0.0  ||  OPT_MinUnique > 100.0 )
    {
      cerr << "ERROR: Minimum uniqueness must be within the range [0, 100]\n";
      errflg ++;
    }

  if ( OPT_MaxOverlap < 0.0  ||  OPT_MaxOverlap > 100.0 )
    {
      cerr << "ERROR: Maximum overlap must be within the range [0, 100]\n";
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




//------------------------------------------------------------- PrintHelp ----//
void PrintHelp (const char * s)
{
  PrintUsage (s);
  cerr
    << "-e float      For switches -g -r -q, keep repeats within e percent\n"
    << "              of the best LIS score [0, 100], no repeats by default\n"
    << "-g            Global alignment using length*identity weighted LIS.\n"
    << "              For every reference-query pair, leave only the aligns\n"
    << "              which form the longest mutually consistent set\n"
    << "-h            Display help information\n"
    << "-i float      Set the minimum alignment identity [0, 100], default "
    << OPT_MinIdentity << endl
    << "-l int        Set the minimum alignment length, default "
    << OPT_MinLength << endl
    << "-q            Query alignment using length*identity weighted LIS.\n"
    << "              For each query, leave only the aligns which form the\n"
    << "              longest consistent set for the query\n"
    << "-r            Reference alignment using length*identity weighted LIS.\n"
    << "              For each reference, leave only the aligns which form\n"
    << "              the longest consistent set for the reference\n"
    << "-u float      Set the minimum alignment uniqueness, i.e. percent of\n"
    << "              the alignment matching to unique reference AND query\n"
    << "              sequence [0, 100], default "
    << OPT_MinUnique << endl
    << "-o float      Set the maximum alignment overlap for -r and -q options\n"
    << "              as a percent of the alignment length [0, 100], default "
    << OPT_MaxOverlap << endl
    << endl;

  cerr
    << "  Reads a delta alignment file from either nucmer or promer and\n"
    << "filters the alignments based on the command-line switches, leaving\n"
    << "only the desired alignments which are output to stdout in the same\n"
    << "delta format as the input. For multiple switches, order of operations\n"
    << "is as follows: -i -l -u -q -r -g. If an alignment is excluded by a\n"
    << "preceding operation, it will be ignored by the succeeding operations\n"
    << "  An important distinction between the -g option and the -r -q\n"
    << "options is that -g requires the alignments to be mutually consistent\n"
    << "in their order, while the -r -q options are not required to be\n"
    << "mutually consistent and therefore tolerate translocations,\n"
    << "inversions, etc. Thus, -r provides a one-to-many, -q a many-to-one,\n"
    << "-r -q a one-to-one local mapping, and -g a one-to-one global mapping\n"
    << "of reference and query bases respectively.\n"
    << endl;

  return;
}




//------------------------------------------------------------ PrintUsage ----//
void PrintUsage (const char * s)
{
  cerr
    << "\nUSAGE: " << s << "  [options]  <deltafile>\n\n";
  return;
}

//-----------------------------------------------------------------------------
//   Programmer: Adam M Phillippy, University of Maryland
//         File: show-diff.cc
//         Date: 12 / 01 / 2006
//
//        Usage: show-diff [options] <deltafile>
//               Try 'show-diff -h' for more information
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
string  OPT_AlignName;            // delta file name
bool    OPT_AMOS    = false;      // AMOS output
bool    OPT_RefDiff = true;       // output reference diff only
bool    OPT_QryDiff = true;       // output query diff only


//=========================================================== Declarations ====
struct EdgeletLoQCmp_t
//!< Sorts query by lo coord, lo to hi
{
  bool operator()(const DeltaEdgelet_t * i, const DeltaEdgelet_t * j) const
  { return ( i->loQ < j->loQ ); }
};

struct EdgeletIdQLoQCmp_t
//!< Sorts query by ID and lo coord, lo to hi
{
  bool operator()(const DeltaEdgelet_t * i, const DeltaEdgelet_t * j) const
  {
    if ( i->edge && j->edge )
      {
        if ( i->edge->qrynode < j->edge->qrynode )
          return true;
        else if ( i->edge->qrynode > j->edge->qrynode )
          return false;
      }
    else if ( !i->edge && j->edge )
      return true;
    else if ( i->edge && !j->edge )
      return false;
    return ( i->loQ < j->loQ );
  }
};

struct EdgeletLoRCmp_t
//!< Sorts reference by lo coord, lo to hi
{
  bool operator()(const DeltaEdgelet_t * i, const DeltaEdgelet_t * j) const
  { return ( i->loR < j->loR ); }
};

struct EdgeletIdRLoRCmp_t
//!< Sorts reference by ID and lo coord, lo to hi
{
  bool operator()(const DeltaEdgelet_t * i, const DeltaEdgelet_t * j) const
  {
    if ( i->edge && j->edge )
      {
        if ( i->edge->refnode < j->edge->refnode )
          return true;
        else if ( i->edge->refnode > j->edge->refnode )
          return false;
      }
    else if ( !i->edge && j->edge )
      return true;
    else if ( i->edge && !j->edge )
      return false;
    return ( i->loR < j->loR );
  }
};


void PrintDiff(DeltaGraph_t & graph);
void PrintNewSeq(const char* seq);
void PrintGap(const char* seq, long s, long e);
void PrintSeqJmp(const char* seq,
                 const char* seqp, const char* seqn,
                 long s, long e);
void PrintLisJmp(const char* seq, long s, long e);
void PrintInv(const char* seq, long s, long e);
void PrintIndel(const char* seq, long s, long e, long gap1, long gap2);
void PrintDup(const char* seq, long s, long e);
void ParseArgs(int argc, char ** argv);
void PrintHelp(const char * s);
void PrintUsage(const char * s);


//============================================================ Definitions ====
//------------------------------------------------------------------- main ----
int main(int argc, char **argv)
{
  DeltaGraph_t graph;

  ParseArgs(argc, argv);

  //-- Get M-to-M alignment, i.e. union of QLIS and RLIS
  graph.build(OPT_AlignName, false);
  graph.flagMtoM();
  graph.clean();

  //-- Output diff
  PrintDiff(graph);

  return EXIT_SUCCESS;
}


//-------------------------------------------------------------- PrintDiff ----
void PrintDiff(DeltaGraph_t & graph)
{
  const char* refid;
  const char* qryid;
  long i,j;
  long nAligns, gapR, gapQ;
  DeltaEdgelet_t lpad, rpad;          // padding for the alignment vector
  lpad.isRLIS = rpad.isRLIS = true;
  lpad.isQLIS = rpad.isQLIS = true;
  lpad.loR = lpad.hiR = lpad.loQ = lpad.hiQ = 0;

  DeltaEdgelet_t *A, *PA, *PGA;       // alignment, prev, prev global
  vector<DeltaEdgelet_t *> aligns;

  map<string, DeltaNode_t>::const_iterator mi;
  vector<DeltaEdge_t *>::const_iterator ei;
  vector<DeltaEdgelet_t *>::iterator eli;

  //-- For each reference sequence
  if ( OPT_RefDiff )
    for ( mi = graph.refnodes.begin(); mi != graph.refnodes.end(); ++ mi )
      {
        //-- Announce new reference sequence
        refid = mi->first.c_str();
        PrintNewSeq(refid);

        //-- Collect all alignments for this reference sequence
        aligns.clear();
        for ( ei  = (mi->second).edges.begin();
              ei != (mi->second).edges.end(); ++ei )
          for ( eli  = (*ei)->edgelets.begin();
                eli != (*ei)->edgelets.end(); ++eli )
            aligns.push_back(*eli);

        //-- Pad the front and back of the alignment vector
        rpad.loR = rpad.hiR = mi->second.len + 1;
        rpad.loQ = rpad.hiQ = LONG_MAX;
        aligns.push_back(&lpad);
        aligns.push_back(&rpad);

        nAligns = aligns.size();

        //-- OVERRIDE *stpc* value with loQ QLIS ordering
        sort(aligns.begin(), aligns.end(), EdgeletIdQLoQCmp_t());
        for ( i = 0, j = 0; i != nAligns; ++i )
          aligns[i]->stpc = aligns[i]->isQLIS ? j++ : -1;

        //-- Sort by reference order
        sort(aligns.begin(), aligns.end(), EdgeletLoRCmp_t());
        assert ( aligns[0] == &lpad && aligns[nAligns-1] == &rpad );

        //-- Walk reference cover alignments, low to high
        PA = PGA = aligns[0];
        for ( i = 1; i != nAligns; ++i )
          {
            //-- Only interested in reference covering alignments
            if ( !aligns[i]->isRLIS ) continue;

            A = aligns[i];
            gapR = A->loR - PA->hiR - 1;

            //-- Reached pad, end of alignments
            if ( A->edge == NULL )
              {
                PrintGap(refid, PA->hiR, A->loR);
              }
            //-- 1-to-1 alignment
            else if ( A->isQLIS && A->edge == PGA->edge )
              {
                //-- Jump within Q
                if ( A->slope() != PGA->slope() ||
                     A->stpc != PGA->stpc + PGA->slope() )
                  {
                    if ( A->slope() == PGA->slope() )
                      PrintLisJmp(refid, PA->hiR, A->loR);
                    else
                      PrintInv(refid, PA->hiR, A->loR);
                  }
                //-- Lined up, nothing in between
                else if ( PA == PGA )
                  {
                    gapQ = A->isPositive() ?
                      A->loQ - PGA->hiQ - 1 :
                      PGA->loQ - A->hiQ - 1;
                    PrintIndel(refid, PA->hiR, A->loR, gapR, gapQ);
                  }
                //-- Lined up, duplication in between
                else
                  {
                    PrintGap(refid, PA->hiR, A->loR);
                  }
              }
            //-- Not in QLIS? Must be a duplication in R
            else if ( !A->isQLIS )
              {
                PrintGap(refid, PA->hiR, A->loR);
                PrintDup(refid, A->loR, A->hiR);
              }
            //-- A->edge != PGA->edge? Jump to different query sequence
            else if ( PGA->edge != NULL )
              {
                PrintSeqJmp(refid,
                            PGA->edge->qrynode->id->c_str(),
                            A->edge->qrynode->id->c_str(),
                            PA->hiR, A->loR);
              }
            //-- Gap before first alignment
            else
              {
                PrintGap(refid, PA->hiR, A->loR);
              }

            if ( A->isQLIS )
              PGA = A;
            PA = A;
          }
      }


  //---------- WARNING! Same code as above but Q's for R's and R's for Q's
  if ( OPT_QryDiff )
    for ( mi = graph.qrynodes.begin(); mi != graph.qrynodes.end(); ++ mi )
      {
        qryid = mi->first.c_str();
        PrintNewSeq(qryid);

        aligns.clear();
        for ( ei  = (mi->second).edges.begin();
              ei != (mi->second).edges.end(); ++ei )
          for ( eli  = (*ei)->edgelets.begin();
                eli != (*ei)->edgelets.end(); ++eli )
            aligns.push_back(*eli);

        rpad.loQ = rpad.hiQ = mi->second.len + 1;
        rpad.loR = rpad.hiR = LONG_MAX;
        aligns.push_back(&lpad);
        aligns.push_back(&rpad);

        nAligns = aligns.size();

        sort(aligns.begin(), aligns.end(), EdgeletIdRLoRCmp_t());
        for ( i = 0, j = 0; i != nAligns; ++i )
          aligns[i]->stpc = aligns[i]->isRLIS ? j++ : -1;

        sort(aligns.begin(), aligns.end(), EdgeletLoQCmp_t());
        assert ( aligns[0] == &lpad && aligns[nAligns-1] == &rpad );

        PA = PGA = aligns[0];
        for ( i = 1; i != nAligns; ++i )
          {
            if ( !aligns[i]->isQLIS ) continue;

            A = aligns[i];
            gapQ = A->loQ - PA->hiQ - 1;

            if ( A->edge == NULL )
              {
                PrintGap(qryid, PA->hiQ, A->loQ);
              }
            else if ( A->isRLIS && A->edge == PGA->edge )
              {
                if ( A->slope() != PGA->slope() ||
                     A->stpc != PGA->stpc + PGA->slope() )
                  {
                    if ( A->slope() == PGA->slope() )
                      PrintLisJmp(qryid, PA->hiQ, A->loQ);
                    else
                      PrintInv(qryid, PA->hiQ, A->loQ);
                  }
                else if ( PA == PGA )
                  {
                    gapR = A->isPositive() ?
                      A->loR - PGA->hiR - 1 :
                      PGA->loR - A->hiR - 1;
                    PrintIndel(qryid, PA->hiQ, A->loQ, gapQ, gapR);
                  }
                else
                  {
                    PrintGap(qryid, PA->hiQ, A->loQ);
                  }
              }
            else if ( !A->isRLIS )
              {
                PrintGap(qryid, PA->hiQ, A->loQ);
                PrintDup(qryid, A->loQ, A->hiQ);
              }
            else if ( PGA->edge != NULL )
              {
                PrintSeqJmp(qryid,
                            PA->edge->refnode->id->c_str(),
                            A->edge->refnode->id->c_str(),
                            PA->hiQ, A->loQ);
              }
            else
              {
                PrintGap(qryid, PA->hiQ, A->loQ);
              }

            if ( A->isRLIS )
              PGA = A;
            PA = A;
          }
      }
}


void PrintNewSeq(const char* seq)
{
  if ( !OPT_AMOS )
    printf(">%s\n", seq);
}


void PrintGap(const char* seq, long s, long e)
{
  if ( e-s-1 == 0 ) return;

  if ( !OPT_AMOS )
    printf("GAP\t%ld\t%ld\t%ld\n",
           s, e, e-s-1);
  else
    printf("%s\tA\t%ld\t%ld\tGAP\t%ld\t%ld\t%ld\n",
           seq, s, e, s, e, e-s-1);
}

void PrintSeqJmp(const char* seq,
                 const char* seqp,
                 const char* seqn,
                 long s, long e)
{
  if ( !OPT_AMOS )
    printf("SEQ\t%ld\t%ld\t%ld\t%s\t%s\n",
           s, e, e-s-1, seqp, seqn);
  else
    printf("%s\tA\t%ld\t%ld\tSEQ\t%ld\t%ld\t%ld\t%s\t%s\n",
           seq, s, e, s, e, e-s-1, seqp, seqn);
}

void PrintLisJmp(const char* seq, long s, long e)
{
  if ( !OPT_AMOS )
    printf("JMP\t%ld\t%ld\t%ld\n",
           s, e, e-s-1);
  else
    printf("%s\tA\t%ld\t%ld\tJMP\t%ld\t%ld\t%ld\n",
           seq, s, e, s, e, e-s-1);
}

void PrintInv(const char* seq, long s, long e)
{
  if ( !OPT_AMOS )
    printf("INV\t%ld\t%ld\t%ld\n", s, e, e-s-1);
  else
    printf("%s\tA\t%ld\t%ld\tINV\t%ld\t%ld\t%ld\n",
           seq, s, e, s, e, e-s-1);
}

void PrintIndel(const char* seq, long s, long e, long gap1, long gap2)
{
  if ( !OPT_AMOS )
    printf("%s\t%ld\t%ld\t%ld\t%ld\t%ld\n",
           (gap1-gap2 > 0 ? "INS":"DEL"), s, e, gap1, gap2, gap1-gap2);
  else
    printf("%s\tA\t%ld\t%ld\t%s\t%ld\t%ld\t%ld\t%ld\t%ld\n",
           seq,s,e,(gap1-gap2 > 0 ? "INS":"DEL"), s, e, gap1, gap2, gap1-gap2);
}

void PrintDup(const char* seq, long s, long e)
{
  if ( !OPT_AMOS )
    printf("DUP\t%ld\t%ld\t%ld\n",
           s, e, e-s+1);
  else
    printf("%s\tA\t%ld\t%ld\tDUP\t%ld\t%ld\t%ld\n",
           seq, s, e, s, e, e-s+1);
}


//-------------------------------------------------------------- ParseArgs ----
void ParseArgs (int argc, char ** argv)
{
  int ch, errflg = 0;
  optarg = NULL;

  while ( !errflg  &&
          ((ch = getopt (argc, argv, "fhqr")) != EOF) )
    switch (ch)
      {
      case 'f':
        OPT_AMOS = true;
        break;

      case 'h':
        PrintHelp (argv[0]);
        exit (EXIT_SUCCESS);
        break;

      case 'q':
        OPT_RefDiff = false;
        break;
        
      case 'r':
        OPT_QryDiff = false;
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
    << "-f            Output diff information as AMOS features\n"
    << "-h            Display help information\n"
    << "-q            Show break information for queries only\n"
    << "-r            Show break information for references only\n"
    << endl;

  cerr
    << "  Outputs a list of structural differences for each sequence in\n"
    << "the reference and query, sorted by position. For a reference\n"
    << "sequence R, and its matching query sequence Q, differences are\n"
    << "categorized as SEQ (jump to new Q), GAP (unaligned gap in R),\n"
    << "JMP (rearrangement), INV (inversion and rearrangement), INS\n"
    << "(insertion into R), DEL (deletion from R), and DUP (duplicate\n"
    << "insertion into R). The first three columns of output are feature\n"
    << "type, R feature start, and R feature end. Additional columns are\n"
    << "added depending on the feature type. Negative feature lengths\n"
    << "indicate adjacent alignment blocks overlap.\n"
    << "  SEQ gap-start gap-end gap-length new-sequence\n"
    << "  GAP gap-start gap-end gap-length\n"
    << "  JMP gap-start gap-end gap-length\n"
    << "  INV gap-start gap-end gap-length\n"
    << "  INS gap-start gap-end gap-lengthR gap-lengthQ gap-diff\n"
    << "  DEL gap-start gap-end gap-lengthR gap-lengthQ gap-diff\n"
    << "  DUP dup-start dup-end dup-length\n"
    << "Feature lists are headed by a > line with the ID of the sequence\n"
    << "being described. Positions always reference this sequence. The\n"
    << "sum of the third column (ignoring negative values) for a sequence\n"
    << "R is the total amount of sequence inserted into R with respect\n"
    << "to Q. Summing the third column after removing DUP features is the\n"
    << "total amount of unique sequence in R.\n"
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

//------------------------------------------------------------------------------
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
//------------------------------------------------------------------------------

#include "delta.hh"
#include "tigrinc.hh"
#include "translate.hh"
#include "sw_alignscore.hh"
#include <vector>
#include <algorithm>
#include <string>
#include <sstream>
#include <cstring>
#include <map>
#include <set>
#include <cstdio>
using namespace std;




//=============================================================== Options ====//
string  OPT_AlignName;                  // delta file name
string  OPT_ReferenceName;              // reference sequence file name
string  OPT_QueryName;                  // query sequence file name

bool    OPT_SortReference = false;      // -r option
bool    OPT_SortQuery     = false;      // -q option
bool    OPT_ShowLength    = false;      // -l option
bool    OPT_ShowConflict  = false;      // -c option
bool    OPT_ShowIndels    = true;       // -I option
bool    OPT_PrintTabular  = false;      // -T option
bool    OPT_PrintHeader   = true;       // -H option
bool    OPT_SelectAligns  = false;      // -S option

int     OPT_Context       = 0;          // -x option

set<string> OPT_Aligns;                 // -S option




//============================================================= Constants ====//
typedef unsigned char DataType_t;
const DataType_t   NULL_DATA = 0;
const DataType_t NUCMER_DATA = 1;
const DataType_t PROMER_DATA = 2;

typedef unsigned char Dir_t;
const Dir_t FORWARD_DIR = 0;
const Dir_t REVERSE_DIR = 1;

const char  INDEL_CHAR = '.';
const char SEQEND_CHAR = '-';



struct DeltaEdgelet_t;
struct DeltaEdge_t;
struct DeltaNode_t;

struct SNP_t
     //!< A single nuc/aa poly
{
  char cQ, cR;
  long int pQ, pR;
  int conQ, conR;
  string ctxQ, ctxR;
  DeltaEdgelet_t * lp;
  DeltaEdge_t * ep;

  SNP_t ( )
  {
    cQ = cR = 0;
    pQ = pR = 0;
    conQ = conR = 0;
  };
};


//-- Data Structure Explanation
//
//   A bipartite graph with two partite sets, R and Q, where R is the set of
//   reference sequences and Q is the set of query sequences. These nodes are
//   named "DeltaNode_t". We connect a node in R to a node in Q if an alignment
//   is present between the two sequences. The group of all alignments between
//   the two is named "DeltaEdge_t" and a single alignment between the two is
//   named a "DeltaEdgelet_t". Alignment coordinates reference the forward
//   strand and are stored lo before hi.
//
struct DeltaEdgelet_t
     //!< A piece of a delta graph edge, a single alignment
{
  unsigned char dirR   : 1;      // reference match direction
  unsigned char dirQ   : 1;      // query match direction

  long int loQ, hiQ, loR, hiR;   // alignment bounds
  int frmQ, frmR;

  vector<long int> delta;        // delta information
  vector<SNP_t> snps;            // snps for this edgelet

  DeltaEdgelet_t ( )
  {
    dirR = dirQ = FORWARD_DIR;
    loQ = hiQ = loR = hiR = 0;
    frmQ = frmR = 1;
  }
};


struct DeltaEdge_t
     //!< A delta graph edge, alignments between a single reference and query
{
  DeltaNode_t * refnode;      // the adjacent reference node
  DeltaNode_t * qrynode;      // the adjacent query node
  vector<DeltaEdgelet_t *> edgelets;     // the set of individual alignments

  DeltaEdge_t ( )
  {
    refnode = qrynode = NULL;
  }

  ~DeltaEdge_t ( )
  {
    vector<DeltaEdgelet_t *>::iterator i;
    for ( i = edgelets . begin( ); i != edgelets . end( ); ++ i )
      delete (*i);
  }
};


struct DeltaNode_t
     //!< A delta graph node, contains the sequence information
{
  const string * id;               // the id of the sequence
  char * seq;                      // the DNA sequence
  long int len;                    // the length of the sequence
  vector<DeltaEdge_t *> edges;     // the set of related edges

  DeltaNode_t ( )
  {
    id = NULL;
    seq = NULL;
    len = 0;
  }

  ~DeltaNode_t ( )
  {
    free (seq);
    // DeltaGraph_t will take care of destructing the edges
  }
};


struct DeltaGraph_t
    //!< A delta graph of sequences and their alignments
{
  //-- The reference and query delta graph nodes (1 node per sequence)
  map<string, DeltaNode_t> refnodes;
  map<string, DeltaNode_t> qrynodes;

  string refpath;
  string qrypath;
  DataType_t datatype;

  DeltaGraph_t ( )
  {
    datatype = NULL_DATA;
  }

  ~DeltaGraph_t ( )
  {
    //-- Make sure the edges only get destructed once
    map<string, DeltaNode_t>::iterator i;
    vector<DeltaEdge_t *>::iterator j;
    for ( i = refnodes . begin( ); i != refnodes . end( ); ++ i )
      for ( j  = i -> second . edges . begin( );
            j != i -> second . edges . end( ); ++ j )
        delete (*j);
  }
};


struct SNP_R_Sort
{
  bool operator() (const SNP_t * a, const SNP_t * b)
  {
    int i = a->ep->refnode->id->compare (*(b->ep->refnode->id));

    if ( i < 0 )
      return true;
    else if ( i > 0 )
      return false;
    else
      {
        if ( a -> pR < b -> pR )
          return true;
        else if ( a -> pR > b -> pR )
          return false;
        else
          {
            int j = a->ep->qrynode->id->compare (*(b->ep->qrynode->id));

            if ( j < 0 )
              return true;
            else if ( j > 0 )
              return false;
            else
              {
                if ( a -> pQ < b -> pQ )
                  return true;
                else
                  return false;
              }
          }
      }
  }
};


struct SNP_Q_Sort
{
  bool operator() (const SNP_t * a, const SNP_t * b)
  {
    int i = a->ep->qrynode->id->compare (*(b->ep->qrynode->id));

    if ( i < 0 )
      return true;
    else if ( i > 0 )
      return false;
    else
      {
        if ( a -> pQ < b -> pQ )
          return true;
        else if ( a -> pQ > b -> pQ )
          return false;
        else
          {
            int j = a->ep->refnode->id->compare (*(b->ep->refnode->id));

            if ( j < 0 )
              return true;
            else if ( j > 0 )
              return false;
            else
              {
                if ( a -> pR < b -> pR )
                  return true;
                else
                  return false;
              }
          }
      }
  }
};




//========================================================== Fuction Decs ====//
//------------------------------------------------------------------ RevC ----//
inline long int RevC (long int coord, long int len)
{
  return len - coord + 1;
}


//------------------------------------------------------------------ Norm ----//
inline long int Norm (long int c, long int l, int f, DataType_t d)
{
  long int retval = (d == PROMER_DATA ? c * 3 - (3 - abs(f)) : c);
  if ( f < 0 ) retval = RevC (retval, l);
  return retval;
}


//------------------------------------------------------------------ Swap ----//
inline void Swap (long int & a, long int & b)
{
  long int t = a; a = b; b = t;
}


//------------------------------------------------------------- BuildEdge ----//
void BuildEdge (DeltaEdge_t & edge, const DeltaRecord_t & rec);


//------------------------------------------------------------ BuildGraph ----//
void BuildGraph (DeltaGraph_t & graph);


//------------------------------------------------------------- CheckSNPs ----//
void CheckSNPs (DeltaGraph_t & graph);


//-------------------------------------------------------------- FindSNPs ----//
void FindSNPs (DeltaGraph_t & graph);


//------------------------------------------------------------ PrintHuman ----//
void PrintHuman (const vector<const SNP_t *> & snps,
                 const DeltaGraph_t & graph);


//---------------------------------------------------------- PrintTabular ----//
void PrintTabular (const vector<const SNP_t *> & snps,
                   const DeltaGraph_t & graph);


//--------------------------------------------------------- ReadSequences ----//
void ReadSequences (DeltaGraph_t & graph);


//---------------------------------------------------------- SelectAligns ----//
void SelectAligns ( );


//------------------------------------------------------------- ParseArgs ----//
void ParseArgs (int argc, char ** argv);


//------------------------------------------------------------- PrintHelp ----//
void PrintHelp (const char * s);


//------------------------------------------------------------ PrintUsage ----//
void PrintUsage (const char * s);




//========================================================= Function Defs ====//
int main (int argc, char ** argv)
{
  vector<const SNP_t *> snps;
  DeltaGraph_t graph;


  //-- Command line parsing
  ParseArgs (argc, argv);

  //-- Select alignments
  if ( OPT_SelectAligns )
    SelectAligns ( );

  //-- Build the alignment graph from the delta file
  BuildGraph (graph);

  //-- Read sequences
  ReadSequences (graph);

  //-- Locate the SNPs
  FindSNPs (graph);

  //-- Check for ambiguous alignment regions
  CheckSNPs (graph);


  //-- Collect and sort the SNPs
  map<string, DeltaNode_t>::iterator mi;
  vector<DeltaEdge_t *>::iterator ei;
  vector<DeltaEdgelet_t *>::iterator li;
  vector<SNP_t>::iterator si;
  for ( mi = graph.refnodes.begin( ); mi != graph.refnodes.end( ); ++ mi )
    for ( ei = mi->second.edges.begin( ); ei != mi->second.edges.end( ); ++ ei )
      for ( li = (*ei)->edgelets.begin( ); li != (*ei)->edgelets.end( ); ++ li )
        for ( si = (*li)->snps.begin( ); si != (*li)->snps.end( ); ++ si )
          snps . push_back (&(*si));

  if ( OPT_SortReference )
    sort (snps . begin( ), snps . end( ), SNP_R_Sort( ));
  else
    sort (snps . begin( ), snps . end( ), SNP_Q_Sort( ));


  //-- Output data to stdout
  if ( OPT_PrintTabular )
    PrintTabular (snps, graph);
  else
    PrintHuman (snps, graph);


  return EXIT_SUCCESS;
}




//------------------------------------------------------------- BuildEdge ----//
void BuildEdge (DeltaEdge_t & edge, const DeltaRecord_t & rec)
{
  DeltaEdgelet_t * p;

  vector<DeltaAlignment_t>::const_iterator i;
  for ( i = rec . aligns . begin( ); i != rec . aligns . end( ); ++ i )
    {
      //-- Set the edgelet
      p = new DeltaEdgelet_t( );

      p -> loR = i -> sR;
      p -> hiR = i -> eR;
      p -> loQ = i -> sQ;
      p -> hiQ = i -> eQ;

      p -> dirR = p -> hiR < p -> loR ? REVERSE_DIR : FORWARD_DIR;
      p -> dirQ = p -> hiQ < p -> loQ ? REVERSE_DIR : FORWARD_DIR;

      //-- Store the delta information
      p -> delta = i -> deltas;

      //-- Force loR < hiR && loQ < hiQ
      if ( p -> dirR == REVERSE_DIR )
        Swap (p -> loR, p -> hiR);
      if ( p -> dirQ == REVERSE_DIR )
        Swap (p -> loQ, p -> hiQ);

      edge . edgelets . push_back (p);
    }
}




//------------------------------------------------------------ BuildGraph ----//
void BuildGraph (DeltaGraph_t & graph)
{
  DeltaReader_t dr;
  DeltaEdge_t * dep;
  pair<map<string, DeltaNode_t>::iterator, bool> insret;


  //-- Open the delta file and read in the alignment information
  dr . open (OPT_AlignName);

  graph . refpath = dr . getReferencePath( );
  graph . qrypath = dr . getQueryPath( );

  if ( dr . getDataType( ) == NUCMER_STRING )
    graph . datatype = NUCMER_DATA;
  else if ( dr . getDataType( ) == PROMER_STRING )
    graph . datatype = PROMER_DATA;
  else
    graph . datatype = NULL_DATA;

  //-- Read in the next graph edge, i.e. a new delta record
  while ( dr . readNext( ) )
    {
      dep = new DeltaEdge_t( );

      //-- Build the edge
      BuildEdge (*dep, dr . getRecord( ));

      if ( dep -> edgelets . empty( ) )
        {
          delete dep;
          continue;
        }


      //-- Find the reference node in the graph, add a new one if necessary
      insret = graph . refnodes . insert
        (map<string, DeltaNode_t>::value_type
         (dr . getRecord( ) . idR, DeltaNode_t( )));
      dep -> refnode = &((insret . first) -> second);

      //-- If a new reference node
      if ( insret . second )
        {
          dep -> refnode -> id  = &((insret . first) -> first);
          dep -> refnode -> len = dr . getRecord( ) . lenR;
        }


      //-- Find the query node in the graph, add a new one if necessary
      insret = graph . qrynodes . insert
        (map<string, DeltaNode_t>::value_type
         (dr . getRecord( ) . idQ, DeltaNode_t( )));
      dep -> qrynode = &((insret . first) -> second);

      //-- If a new query node
      if ( insret . second )
        {
          dep -> qrynode -> id  = &((insret . first) -> first);
          dep -> qrynode -> len = dr . getRecord( ) . lenQ;
        }


      //-- Link the nodes
      dep -> refnode -> edges . push_back (dep);
      dep -> qrynode -> edges . push_back (dep);
    }
  dr . close ( );
}




//------------------------------------------------------------- CheckSNPs ----//
void CheckSNPs (DeltaGraph_t & graph)
{
  map<string, DeltaNode_t>::const_iterator mi;
  vector<DeltaEdge_t *>::const_iterator ei;
  vector<DeltaEdgelet_t *>::iterator eli;
  vector<SNP_t>::iterator si;
  long int i;

  //-- For each reference sequence
  long int ref_size = 0;
  long int ref_len = 0;
  unsigned char * ref_cov = NULL;
  for ( mi = graph.refnodes.begin( ); mi != graph.refnodes.end( ); ++ mi )
    {
      //-- Reset the reference coverage array
      ref_len = (mi -> second) . len;
      if ( ref_len > ref_size )
        {
          ref_cov = (unsigned char *) Safe_realloc (ref_cov, ref_len + 1);
          ref_size = ref_len;
        }
      for ( i = 1; i <= ref_len; ++ i )
        ref_cov[i] = 0;

      //-- Add to the reference coverage
      for ( ei  = (mi -> second) . edges . begin( );
            ei != (mi -> second) . edges . end( ); ++ ei )
        for ( eli  = (*ei) -> edgelets . begin( );
              eli != (*ei) -> edgelets . end( ); ++ eli )
          for ( i = (*eli) -> loR; i <= (*eli) -> hiR; i ++ )
            if ( ref_cov [i] < UCHAR_MAX )
              ref_cov [i] ++;

      //-- Set the SNP conflict counter
      for ( ei  = (mi -> second) . edges . begin( );
            ei != (mi -> second) . edges . end( ); ++ ei )
        for ( eli  = (*ei) -> edgelets . begin( );
              eli != (*ei) -> edgelets . end( ); ++ eli )
          for ( si = (*eli)->snps.begin( ); si != (*eli)->snps.end( ); ++ si )
            si -> conR = ref_cov [si->pR] - 1;
    }
  free (ref_cov);


  //-- For each query sequence
  long int qry_size = 0;
  long int qry_len = 0;
  unsigned char * qry_cov = NULL;
  for ( mi = graph.qrynodes.begin( ); mi != graph.qrynodes.end( ); ++ mi )
    {
      //-- Reset the query coverage array
      qry_len = (mi -> second) . len;
      if ( qry_len > qry_size )
        {
          qry_cov = (unsigned char *) Safe_realloc (qry_cov, qry_len + 1);
          qry_size = qry_len;
        }
      for ( i = 1; i <= qry_len; ++ i )
        qry_cov[i] = 0;

      //-- Add to the query coverage
      for ( ei  = (mi -> second) . edges . begin( );
            ei != (mi -> second) . edges . end( ); ++ ei )
        for ( eli  = (*ei) -> edgelets . begin( );
              eli != (*ei) -> edgelets . end( ); ++ eli )
          for ( i = (*eli) -> loQ; i <= (*eli) -> hiQ; i ++ )
            if ( qry_cov [i] < UCHAR_MAX )
              qry_cov [i] ++;

      //-- Set the SNP conflict counter
      for ( ei  = (mi -> second) . edges . begin( );
            ei != (mi -> second) . edges . end( ); ++ ei )
        for ( eli  = (*ei) -> edgelets . begin( );
              eli != (*ei) -> edgelets . end( ); ++ eli )
          for ( si = (*eli)->snps.begin( ); si != (*eli)->snps.end( ); ++ si )
            si -> conQ = qry_cov [si->pQ] - 1;
    }
  free (qry_cov);
}




//-------------------------------------------------------------- FindSNPs ----//
void FindSNPs (DeltaGraph_t & graph)
{
  map<string, DeltaNode_t>::iterator mi;
  vector<DeltaEdge_t *>::iterator ei;
  vector<DeltaEdgelet_t *>::iterator li;

  //-- For each alignment, identify the SNPs
  for ( mi = graph.refnodes.begin( ); mi != graph.refnodes.end( ); ++ mi )
    for ( ei = mi->second.edges.begin( ); ei != mi->second.edges.end( ); ++ ei )
      {
        SNP_t snp;
        int ri, qi;
        char * R[] = {(*ei)->refnode->seq, NULL, NULL, NULL, NULL, NULL, NULL};
        char * Q[] = {(*ei)->qrynode->seq, NULL, NULL, NULL, NULL, NULL, NULL};

        long int i;
        long int lenR = (*ei) -> refnode -> len;
        long int lenQ = (*ei) -> qrynode -> len;

        snp . ep = *ei;

        for (li = (*ei)->edgelets.begin( ); li != (*ei)->edgelets.end( ); ++ li)
          {
            long int delta;
            int frameR, frameQ, sign;
            long int sR, eR, sQ, eQ;
            long long int rpos, qpos, remain;
            long long int rctx, qctx;
            vector<long int>::iterator dp;
            long int alenR = lenR;
            long int alenQ = lenQ;

            snp . lp = *li;

            //-- Only do the ones requested by user
            if ( OPT_SelectAligns )
              {
                ostringstream ss;
                set<string>::iterator si;

                if ( (*li) -> dirR == FORWARD_DIR )
                  ss << (*li) -> loR << ' ' << (*li) -> hiR << ' ';
                else
                  ss << (*li) -> hiR << ' ' << (*li) -> loR << ' ';

                if ( (*li) -> dirQ == FORWARD_DIR )
                  ss << (*li) -> loQ << ' ' << (*li) -> hiQ << ' ';
                else
                  ss << (*li) -> hiQ << ' ' << (*li) -> loQ << ' ';

                ss << *((*ei)->refnode->id) << ' ' << *((*ei)->qrynode->id);

                si = OPT_Aligns . find (ss .str( ));
                if ( si == OPT_Aligns . end( ) )
                  continue;
                else
                  OPT_Aligns . erase (si);
              }

            //-- Point the coords the right direction
            frameR = 1;
            if ( (*li) -> dirR == REVERSE_DIR )
              {
                sR = RevC ((*li) -> hiR, lenR);
                eR = RevC ((*li) -> loR, lenR);
                frameR += 3;
              }
            else
              {
                sR = (*li) -> loR;
                eR = (*li) -> hiR;
              }

            frameQ = 1;
            if ( (*li) -> dirQ == REVERSE_DIR )
              {
                sQ = RevC ((*li) -> hiQ, lenQ);
                eQ = RevC ((*li) -> loQ, lenQ);
                frameQ += 3;
              }
            else
              {
                sQ = (*li) -> loQ;
                eQ = (*li) -> hiQ;
              }

            //-- Translate coords to AA if necessary
            if ( graph . datatype == PROMER_DATA )
              {
                alenR /= 3;
                alenQ /= 3;

                frameR += (sR + 2) % 3;
                frameQ += (sQ + 2) % 3;

                // remeber that eR and eQ point to the last base in the codon
                sR = (sR + 2) / 3;
                eR = eR / 3;
                sQ = (sQ + 2) / 3;
                eQ = eQ / 3;
              }

            ri = frameR;
            qi = frameQ;

            if ( frameR > 3 )
              frameR = -(frameR - 3);
            if ( frameQ > 3 )
              frameQ = -(frameQ - 3);

            //-- Generate the sequences if needed
            if ( R [ri] == NULL )
              {
                if ( graph . datatype == PROMER_DATA )
                  {
                    R [ri] = (char *) Safe_malloc (alenR + 2);
                    R [ri][0] = '\0';
                    Translate_DNA (R [0], R [ri], ri);
                  }
                else
                  {
                    R [ri] = (char *) Safe_malloc (alenR + 2);
                    R [ri][0] = '\0';
                    strcpy (R [ri] + 1, R [0] + 1);
                    if ( (*li) -> dirR == REVERSE_DIR )
                      Reverse_Complement (R [ri], 1, lenR);
                  }
              }
            if ( Q [qi] == NULL )
              {
                if ( graph . datatype == PROMER_DATA )
                  {
                    Q [qi] = (char *) Safe_malloc (alenQ + 2);
                    Q [qi][0] = '\0';
                    Translate_DNA (Q [0], Q [qi], qi);
                  }
                else
                  {
                    Q [qi] = (char *) Safe_malloc (alenQ + 2);
                    Q [qi][0] = '\0';
                    strcpy (Q [qi] + 1, Q [0] + 1);
                    if ( (*li) -> dirQ == REVERSE_DIR )
                      Reverse_Complement (Q [qi], 1, lenQ);
                  }
              }

            //-- Locate the SNPs
            rpos = sR;
            qpos = sQ;
            remain = eR - sR + 1;

            (*li) -> frmR = frameR;
            (*li) -> frmQ = frameQ;

            for ( dp  = (*li) -> delta . begin( );
                  dp != (*li) -> delta . end( )  &&  *dp != 0; ++ dp )
              {
                delta = *dp;
                sign = delta > 0 ? 1 : -1;
                delta = labs (delta);

                //-- For all SNPs before the next indel
                for ( i = 1; i < delta; i ++ )
                  if ( R [ri] [rpos ++] != Q [qi] [qpos ++] )
                    {
                      snp . pR = Norm (rpos-1, lenR, frameR, graph.datatype);
                      snp . pQ = Norm (qpos-1, lenQ, frameQ, graph.datatype);
                      snp . cR = toupper (R [ri] [rpos-1]);
                      snp . cQ = toupper (Q [qi] [qpos-1]);

                      snp . ctxR . erase( );
                      for ( rctx = rpos - OPT_Context - 1;
                            rctx < rpos + OPT_Context; rctx ++ )
                        if ( rctx < 1  ||  rctx > alenR )
                          snp . ctxR . push_back (SEQEND_CHAR);
                        else if ( rctx == rpos - 1 )
                          snp . ctxR . push_back (snp . cR);
                        else
                          snp . ctxR . push_back (toupper (R [ri] [rctx]));

                      snp . ctxQ . erase( );
                      for ( qctx = qpos - OPT_Context - 1;
                            qctx < qpos + OPT_Context; qctx ++ )
                        if ( qctx < 1  ||  qctx > alenQ )
                          snp . ctxQ . push_back (SEQEND_CHAR);
                        else if ( qctx == qpos - 1 )
                          snp . ctxQ . push_back (snp . cQ);
                        else
                          snp . ctxQ . push_back (toupper (Q [qi] [qctx]));

                      (*li) -> snps . push_back (snp);
                    }

                //-- For the indel
                snp . ctxR . erase( );
                for ( rctx = rpos - OPT_Context; rctx < rpos; rctx ++ )
                  if ( rctx < 1 )
                    snp . ctxR . push_back (SEQEND_CHAR);
                  else
                    snp . ctxR . push_back (toupper (R [ri] [rctx]));

                snp . ctxQ . erase( );
                for ( qctx = qpos - OPT_Context; qctx < qpos; qctx ++ )
                  if ( qctx < 1 )
                    snp . ctxQ . push_back (SEQEND_CHAR);
                  else
                    snp . ctxQ . push_back (toupper (Q [qi] [qctx]));

                if ( sign > 0 )
                  {
                    snp . pR = Norm (rpos, lenR, frameR, graph.datatype);
                    if ( frameQ > 0 )
                      snp . pQ = Norm (qpos - 1, lenQ, frameQ, graph.datatype);
                    else
                      snp . pQ = Norm (qpos, lenQ, frameQ, graph.datatype);

                    snp . cR = toupper (R [ri] [rpos ++]);
                    snp . cQ = INDEL_CHAR;
                    snp . ctxR . push_back (snp . cR);
                    snp . ctxQ . push_back (snp . cQ);
                    remain -= i;
                    rctx ++;
                  }
                else
                  {
                    snp . pQ = Norm (qpos, lenQ, frameQ, graph.datatype);
                    if ( frameR > 0 )
                      snp . pR = Norm (rpos - 1, lenR, frameR, graph.datatype);
                    else
                      snp . pR = Norm (rpos, lenR, frameR, graph.datatype);

                    snp . cR = INDEL_CHAR;
                    snp . cQ = toupper (Q [qi] [qpos ++]);
                    snp . ctxR . push_back (snp . cR);
                    snp . ctxQ . push_back (snp . cQ);
                    remain -= i - 1;
                    qctx ++;
                  }

                for ( ; rctx < rpos + OPT_Context; rctx ++ )
                  if ( rctx > alenR )
                    snp . ctxR . push_back (SEQEND_CHAR);
                  else
                    snp . ctxR . push_back (toupper (R [ri] [rctx]));

                for ( ; qctx < qpos + OPT_Context; qctx ++ )
                  if ( qctx > alenQ )
                    snp . ctxQ . push_back (SEQEND_CHAR);
                  else
                    snp . ctxQ . push_back (toupper (Q [qi] [qctx]));


                (*li) -> snps . push_back (snp);
              }

            //-- For all SNPs after the final indel
            for ( i = 0; i < remain; i ++ )
              if ( R [ri] [rpos ++] != Q [qi] [qpos ++] )
                {
                  snp . pR = Norm (rpos-1, lenR, frameR, graph.datatype);
                  snp . pQ = Norm (qpos-1, lenQ, frameQ, graph.datatype);
                  snp . cR = toupper (R [ri] [rpos-1]);
                  snp . cQ = toupper (Q [qi] [qpos-1]);

                  snp . ctxR . erase( );
                  for ( rctx = rpos - OPT_Context - 1;
                        rctx < rpos + OPT_Context; rctx ++ )
                    if ( rctx < 1  ||  rctx > alenR )
                      snp . ctxR . push_back (SEQEND_CHAR);
                    else if ( rctx == rpos - 1 )
                      snp . ctxR . push_back (snp . cR);
                    else
                      snp . ctxR . push_back (toupper (R [ri] [rctx]));
                      
                  snp . ctxQ . erase( );
                  for ( qctx = qpos - OPT_Context - 1;
                        qctx < qpos + OPT_Context; qctx ++ )
                    if ( qctx < 1  ||  qctx > alenQ )
                      snp . ctxQ . push_back (SEQEND_CHAR);
                    else if ( qctx == qpos - 1 )
                      snp . ctxQ . push_back (snp . cQ);
                    else
                      snp . ctxQ . push_back (toupper (Q [qi] [qctx]));
                
                  (*li) -> snps . push_back (snp);
                }
          }

        //-- Clear up the seq
        for ( i = 1; i <= 6; i ++ )
          {
            free (R[i]);
            free (Q[i]);
          }
      }

  if ( OPT_SelectAligns  &&  ! OPT_Aligns . empty( ) )
    {
      cerr << "ERROR: One or more alignments from stdin could not be found\n";
      exit (EXIT_FAILURE);
    }
}




//------------------------------------------------------------ PrintHuman ----//
void PrintHuman (const vector<const SNP_t *> & snps,
                 const DeltaGraph_t & graph)
{
  vector<const SNP_t *>::const_iterator si, psi, nsi;
  long int diff, dist, distR, distQ;
  int ctxw = 2 * OPT_Context + 1;
  int ctxc = ctxw < 7 ? 7 : ctxw;

  if ( OPT_PrintHeader )
    {
      printf ("%s %s\n%s\n\n",
              graph . refpath . c_str( ), graph . qrypath . c_str( ),
              graph . datatype == NUCMER_DATA ? "NUCMER" : "PROMER");
      printf ("%8s  %5s  %-8s  | ", "[P1]", "[SUB]", "[P2]");
      printf ("%8s %8s  | ", "[DIFF]", "[DIST]");
      if ( OPT_ShowConflict )
        printf ("%4s %4s  | ", "[R]", "[Q]");
      if ( OPT_ShowLength )
        printf ("%8s %8s  | ", "[LEN R]", "[LEN Q]");
      if ( OPT_Context != 0 )
        {
          for ( int i = 0; i < ctxc - 7; i ++ )
            putchar (' ');
          printf (" [CTX R]  ");
          for ( int i = 0; i < ctxc - 7; i ++ )
            putchar (' ');
          printf ("[CTX Q]  | ");
        }
      printf ("%5s  ", "[FRM]");
      printf ("%s", "[TAGS]");
      printf ("\n");

      if ( OPT_ShowConflict )
        printf ("=============");
      if ( OPT_ShowLength )
        printf ("=====================");
      if ( OPT_Context != 0 )
        for ( int i = 0; i < 2 * ctxc + 7; i ++ )
          putchar ('=');
      printf("================================="
             "==================================\n");
    }

  for ( si = snps . begin( ); si != snps . end( ); ++ si )
    {
      psi = si - 1;
      nsi = si + 1;

      //-- Distance to the nearest alignment bound
      distR = (*si)->pR - (*si)->lp->loR < (*si)->lp->hiR - (*si)->pR ?
        (*si)->pR - (*si)->lp->loR : (*si)->lp->hiR - (*si)->pR;
      distQ = (*si)->pQ - (*si)->lp->loQ < (*si)->lp->hiQ - (*si)->pQ ?
        (*si)->pQ - (*si)->lp->loQ : (*si)->lp->hiQ - (*si)->pQ;
      dist = distR < distQ ? distR : distQ;

      //-- Distance to the nearest mismatch
      if ( OPT_SortReference )
        {
          diff = distR;

          if ( psi >= snps . begin( )  &&
               (*psi)->ep->refnode->id == (*si)->ep->refnode->id  &&
               (*si)->pR - (*psi)->pR < diff )
            diff = (*si)->pR - (*psi)->pR;

          if ( nsi < snps . end( )  &&
               (*nsi)->ep->refnode->id == (*si)->ep->refnode->id  &&
               (*nsi)->pR - (*si)->pR < diff )
            diff = (*nsi)->pR - (*si)->pR;
        }
      else
        {
          diff = distQ;

          if ( psi >= snps . begin( )  &&
               (*psi)->ep->refnode->id == (*si)->ep->refnode->id  &&
               (*si)->pQ - (*psi)->pQ < diff )
            diff = (*si)->pQ - (*psi)->pQ;

          if ( nsi < snps . end( )  &&
               (*nsi)->ep->refnode->id == (*si)->ep->refnode->id  &&
               (*nsi)->pQ - (*si)->pQ < diff )
            diff = (*nsi)->pQ - (*si)->pQ;
        }

      //-- Skip this SNP if we are hiding it
      if ( ! OPT_ShowConflict  &&  ((*si)->conR != 0  ||  (*si)->conQ != 0) )
        continue;

      printf ("%8ld   %c %c   %-8ld  | ",
              (*si)->pR, (*si)->cR, (*si)->cQ, (*si)->pQ);
      printf ("%8ld %8ld  | ", diff, dist);
      if ( OPT_ShowConflict )
        printf ("%4d %4d  | ", (*si)->conR, (*si)->conQ);
      if ( OPT_ShowLength )
        printf ("%8ld %8ld  | ",
                (*si)->ep->refnode->len, (*si)->ep->qrynode->len);
      if ( OPT_Context != 0 )
        {
          for ( int i = 0; i < ctxc - ctxw; i ++ )
            putchar (' ');
          printf (" %s  ", (*si)->ctxR.c_str( ));
          for ( int i = 0; i < ctxc - ctxw; i ++ )
            putchar (' ');
          printf ("%s  | ", (*si)->ctxQ.c_str( ));
        }
      printf ("%2d %2d  ", (*si)->lp->frmR, (*si)->lp->frmQ);
      printf ("%s\t%s",
              (*si)->ep->refnode->id->c_str( ),
              (*si)->ep->qrynode->id->c_str( ));
      printf ("\n");
    }
}




//---------------------------------------------------------- PrintTabular ----//
void PrintTabular (const vector<const SNP_t *> & snps,
                   const DeltaGraph_t & graph)
{
  vector<const SNP_t *>::const_iterator si, psi, nsi;
  long int diff, dist, distR, distQ;

  if ( OPT_PrintHeader )
    {
      printf ("%s %s\n%s\n\n",
              graph . refpath . c_str( ), graph . qrypath . c_str( ),
              graph . datatype == NUCMER_DATA ? "NUCMER" : "PROMER");
      printf ("%s\t%s\t%s\t%s\t", "[P1]", "[SUB]", "[SUB]", "[P2]");
      printf ("%s\t%s\t", "[DIFF]", "[DIST]");
      if ( OPT_ShowConflict )
        printf ("%s\t%s\t", "[R]", "[Q]");
      if ( OPT_ShowLength )
        printf ("%s\t%s\t", "[LEN R]", "[LEN Q]");
      if ( OPT_Context != 0 )
        printf ("%s\t%s\t", "[CTX R]", "[CTX Q]");
      printf ("%s\t", "[FRM]");
      printf ("%s\n", "[TAGS]");
    }

  for ( si = snps . begin( ); si != snps . end( ); ++ si )
    {
      psi = si - 1;
      nsi = si + 1;

      //-- Distance to the nearest alignment bound
      distR = (*si)->pR - (*si)->lp->loR < (*si)->lp->hiR - (*si)->pR ?
        (*si)->pR - (*si)->lp->loR : (*si)->lp->hiR - (*si)->pR;
      distQ = (*si)->pQ - (*si)->lp->loQ < (*si)->lp->hiQ - (*si)->pQ ?
        (*si)->pQ - (*si)->lp->loQ : (*si)->lp->hiQ - (*si)->pQ;
      dist = distR < distQ ? distR : distQ;

      //-- Distance to the nearest mismatch
      if ( OPT_SortReference )
        {
          diff = distR;

          if ( psi >= snps . begin( )  &&
               (*psi)->ep->refnode->id == (*si)->ep->refnode->id  &&
               (*si)->pR - (*psi)->pR < diff )
            diff = (*si)->pR - (*psi)->pR;

          if ( nsi < snps . end( )  &&
               (*nsi)->ep->refnode->id == (*si)->ep->refnode->id  &&
               (*nsi)->pR - (*si)->pR < diff )
            diff = (*nsi)->pR - (*si)->pR;
        }
      else
        {
          diff = distQ;

          if ( psi >= snps . begin( )  &&
               (*psi)->ep->refnode->id == (*si)->ep->refnode->id  &&
               (*si)->pQ - (*psi)->pQ < diff )
            diff = (*si)->pQ - (*psi)->pQ;

          if ( nsi < snps . end( )  &&
               (*nsi)->ep->refnode->id == (*si)->ep->refnode->id  &&
               (*nsi)->pQ - (*si)->pQ < diff )
            diff = (*nsi)->pQ - (*si)->pQ;
        }

      //-- Skip this SNP if we are hiding it
      if ( (! OPT_ShowConflict  &&
            ((*si)->conR != 0  ||  (*si)->conQ != 0))
           ||
           (! OPT_ShowIndels  &&
            ((*si)->cR == INDEL_CHAR  ||  (*si)->cQ == INDEL_CHAR)) )
        continue;

      printf ("%ld\t%c\t%c\t%ld\t",
              (*si)->pR, (*si)->cR, (*si)->cQ, (*si)->pQ);
      printf ("%ld\t%ld\t", diff, dist);
      if ( OPT_ShowConflict )
        printf ("%d\t%d\t", (*si)->conR, (*si)->conQ);
      if ( OPT_ShowLength )
        printf ("%ld\t%ld\t",
                (*si)->ep->refnode->len, (*si)->ep->qrynode->len);
      if ( OPT_Context != 0 )
        printf ("%s\t%s\t", (*si)->ctxR.c_str( ), (*si)->ctxQ.c_str( ));
      printf ("%d\t%d\t", (*si)->lp->frmR, (*si)->lp->frmQ);
      printf ("%s\t%s",
              (*si)->ep->refnode->id->c_str( ),
              (*si)->ep->qrynode->id->c_str( ));
      printf ("\n");
    }
}




//--------------------------------------------------------- ReadSequences ----//
void ReadSequences (DeltaGraph_t & graph)
{
  map<string, DeltaNode_t>::iterator mi;

  FILE * qryfile, * reffile;
  char * R = NULL;
  char * Q = NULL;
  char id [MAX_LINE];
  long int initsize;
  long int len;

  //-- Read in the reference sequences
  reffile = File_Open (graph . refpath . c_str( ), "r");
  initsize = INIT_SIZE;
  R = (char *) Safe_malloc (initsize);
  while ( Read_String (reffile, R, initsize, id, FALSE) )
    if ( (mi = graph . refnodes . find (id)) != graph . refnodes . end( ) )
      {
        len = strlen (R + 1);
        mi -> second . seq = (char *) Safe_malloc (len + 2);
        mi -> second . seq[0] = '\0';
        strcpy (mi -> second . seq + 1, R + 1);
        if ( len != mi -> second . len )
          {
            cerr << "ERROR: Reference input does not match delta file\n";
            exit (EXIT_FAILURE);
          }
      }
  fclose (reffile);
  free (R);

  //-- Read in the query sequences
  qryfile = File_Open (graph . qrypath . c_str( ), "r");
  initsize = INIT_SIZE;
  Q = (char *) Safe_malloc (initsize);
  while ( Read_String (qryfile, Q, initsize, id, FALSE) )
    if ( (mi = graph . qrynodes . find (id)) != graph . qrynodes . end( ) )
      {
        len = strlen (Q + 1);
        mi -> second . seq = (char *) Safe_malloc (len + 2);
        mi -> second . seq[0] = '\0';
        strcpy (mi -> second . seq + 1, Q + 1);
        if ( len != mi -> second . len )
          {
            cerr << "ERROR: Query input does not match delta file\n";
            exit (EXIT_FAILURE);
          }
      }
  fclose (qryfile);
  free (Q);


  //-- Check that we found all the sequences
  for ( mi = graph.refnodes.begin( ); mi != graph.refnodes.end( ); ++ mi )
    if ( mi -> second . seq == NULL )
      {
        cerr << "ERROR: '" << mi -> first << "' not found in reference file\n";
        exit (EXIT_FAILURE);
      }

  for ( mi = graph.qrynodes.begin( ); mi != graph.qrynodes.end( ); ++ mi )
    if ( mi -> second . seq == NULL )
      {
        cerr << "ERROR: '" << mi -> first << "' not found in query file\n";
        exit (EXIT_FAILURE);
      }
}




//---------------------------------------------------------- SelectAligns ----//
void SelectAligns ( )
{
  string line, part;
  string s1, e1, s2, e2, tR, tQ;
  istringstream sin;
  ostringstream sout;

  while ( cin )
    {
      getline (cin, line);
      if ( line . empty( ) )
        continue;

      sin . clear( );
      sin . str (line);
      sin >> s1 >> e1 >> s2;
      if ( s2 == "|" ) sin >> s2;
      sin >> e2 >> tR >> tQ;

      if ( sin . fail( ) )
        {
          cerr << "ERROR: Could not parse input from stdin\n";
          exit (EXIT_FAILURE);
        }

      while ( sin >> part )
        {
          tR = tQ;
          tQ = part;
        }

      sout . clear( );
      sout . str ("");
      sout << s1 << ' ' << e1 << ' '
           << s2 << ' ' << e2 << ' '
           << tR << ' ' << tQ;

      OPT_Aligns . insert (sout . str( ));
    }
}




//------------------------------------------------------------- ParseArgs ----//
void ParseArgs (int argc, char ** argv)
{
  int ch, errflg = 0;
  optarg = NULL;
  
  while ( !errflg  &&
          ((ch = getopt (argc, argv, "chHIlqrSTx:")) != EOF) )
    switch (ch)
      {
      case 'c':
        OPT_ShowConflict = true;
        break;

      case 'h':
        PrintHelp (argv[0]);
        exit (EXIT_SUCCESS);
        break;

      case 'H':
        OPT_PrintHeader = false;
        break;

      case 'I':
        OPT_ShowIndels = false;
        break;

      case 'l':
        OPT_ShowLength = true;
        break;

      case 'q':
        OPT_SortQuery = true;
        break;

      case 'r':
        OPT_SortReference = true;
        break;

      case 'S':
        OPT_SelectAligns = true;
        break;

      case 'T':
        OPT_PrintTabular = true;
        break;

      case 'x':
        OPT_Context = atoi (optarg);
        break;

      default:
        errflg ++;
      }

  if ( OPT_Context < 0 )
    {
      cerr << "ERROR: SNP context must be a positive int\n";
      errflg ++;
    }

  if ( OPT_SortReference  &&  OPT_SortQuery )
    cerr << "WARNING: both -r and -q were passed, -q ignored\n";

  if ( !OPT_SortReference  &&  !OPT_SortQuery )
    OPT_SortReference = true;

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
    << "-c            Report SNPs from conflicting regions, i.e. regions of\n"
    << "              the reference or query with more than 1 representative\n"
    << "              alignment. This will also force two additional output\n"
    << "              columns reporting the number of conflicts.\n"
    << "-h            Display help information\n"
    << "-H            Do not print the output header\n"
    << "-I            Do not output indels\n"            
    << "-l            Include sequence length information in the output\n"
    << "-q            Sort output lines by query IDs and SNP positions\n"
    << "-r            Sort output lines by reference IDs and SNP positions\n"
    << "-S            Specify which alignments to report by passing\n"
    << "              'show-coords' lines to stdin\n"
    << "-T            Switch to tab-delimited format\n"
    << "-x int        Include x characters of surrounding SNP context in the\n"
    << "              output, default "
    << OPT_Context << endl
    << endl;

  cerr
    << "  Input is the .delta output of either the nucmer or promer program\n"
    << "passed on the command line.\n"
    << "  Output is to stdout, and consists of a list of SNPs (or amino acid\n"
    << "substitutions for promer) with positions and other useful\n"
    << "information. SNPs will only be reported from regions with an\n"
    << "unambiguous mapping unless the -c option is specified. Output will be\n"
    << "sorted with -r by default and the difference column will always\n"
    << "reference the sequence whose positions have been sorted.\n"
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

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

bool           OPT_QLIS         = false;      // do query based LIS
bool           OPT_RLIS         = false;      // do reference based LIS
bool           OPT_GLIS         = false;      // do global LIS
long int       OPT_MinLength    = 0;          // minimum alignment length
float          OPT_MinIdentity  = 0.0;        // minimum %identity
float          OPT_MinUnique    = 0.0;        // minimum %unique
float          OPT_MaxOverlap   = 75.0;       // maximum olap as % of align len




//============================================================= Constants ====//
typedef unsigned char DataType_t;
const DataType_t NULL_DATA = 0;
const DataType_t NUCMER_DATA = 1;
const DataType_t PROMER_DATA = 2;

typedef unsigned char Dir_t;
const Dir_t FORWARD_DIR = 0;
const Dir_t REVERSE_DIR = 1;


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
struct DeltaNode_t;
struct DeltaEdgelet_t
     //!< A piece of a delta graph edge, a single alignment
{
  unsigned char isGOOD : 1;   // meets the requirements
  unsigned char isQLIS : 1;   // is part of the query's LIS
  unsigned char isRLIS : 1;   // is part of the reference's LIS
  unsigned char isGLIS : 1;   // is part of the reference/query LIS
  unsigned char dirR   : 1;   // reference match direction
  unsigned char dirQ   : 1;   // query match direction

  float idy, sim, stp;                    // percent identity [0 - 1]
  unsigned long int idyc, simc, stpc;     // idy, sim, stp counts
  unsigned long int loQ, hiQ, loR, hiR;   // alignment bounds

  string delta;                           // delta information

  DeltaEdgelet_t ( )
  {
    isGOOD = true;
    isQLIS = isRLIS = isGLIS = false;
    dirR = dirQ = FORWARD_DIR;
    idy = sim = stp = 0;
    idyc = simc = stpc = 0;
    loQ = hiQ = loR = hiR = 0;
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
  unsigned long int len;           // the length of the sequence
  vector<DeltaEdge_t *> edges;     // the set of related edges

  DeltaNode_t ( )
  {
    id = NULL;
    len = 0;
  }

  // DeltaGraph_t will take care of destructing the edges
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


struct EdgeletQCmp_t
    //!< Compares query lo coord
{
  bool operator( ) (const DeltaEdgelet_t * i, const DeltaEdgelet_t * j) const
  { return ( i -> loQ < j -> loQ ); }
};


struct EdgeletRCmp_t
    //!< Compares reference lo coord
{
  bool operator( ) (const DeltaEdgelet_t * i, const DeltaEdgelet_t * j) const
  { return ( i -> loR < j -> loR ); }
};


struct LIS_t
    //!< LIS score
{
  DeltaEdgelet_t * a;
  long int score;
  long int from;
};




//========================================================== Fuction Decs ====//
//------------------------------------------------------------------ RevC ----//
inline unsigned long int RevC (const unsigned long int & coord,
                               const unsigned long int & len)
{
  return len - coord + 1;
}


//------------------------------------------------------------ ScoreLocal ----//
inline long int ScoreLocal
(long int scorej, long int leni, long int lenj, long int olap, float idyi)
{
  if ( olap > 0  &&
       ((float)olap / (float)leni * 100.0 > OPT_MaxOverlap  ||
	(float)olap / (float)lenj * 100.0 > OPT_MaxOverlap) )
    return -1;
  else
    return (scorej + (long int)((leni - olap) * pow (idyi, 2)));
}


//----------------------------------------------------------- ScoreGlobal ----//
inline long int ScoreGlobal
(long int scorej, long int leni, long int olap, float idyi)
{
  return (scorej + (long int)((leni - olap) * pow (idyi, 2)));
}


//------------------------------------------------------------------ Swap ----//
inline void Swap (unsigned long int & a, unsigned long int & b)
{
  unsigned long int t = a; a = b; b = t;
}


//------------------------------------------------------------- BuildEdge ----//
void BuildEdge (DeltaEdge_t & edge, const DeltaRecord_t & rec);


//------------------------------------------------------------ BuildGraph ----//
void BuildGraph (DeltaGraph_t & graph);


//-------------------------------------------------------------- FlagGLIS ----//
void FlagGLIS (DeltaGraph_t & graph);


//------------------------------------------------------------- FlagScore ----//
void FlagScore (DeltaGraph_t & graph);


//-------------------------------------------------------------- FlagQLIS ----//
void FlagQLIS (DeltaGraph_t & graph);


//-------------------------------------------------------------- FlagRLIS ----//
void FlagRLIS (DeltaGraph_t & graph);


//-------------------------------------------------------------- FlagUNIQ ----//
void FlagUNIQ (DeltaGraph_t & graph);


//------------------------------------------------------------ PrintDelta ----//
void PrintDelta (const DeltaGraph_t & graph);


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


  //-- Command line parsing
  ParseArgs (argc, argv);

  //-- Build the alignment graph from the delta file
  BuildGraph (graph);


  //-- Identity requirements
  if ( OPT_MinIdentity > 0  ||  OPT_MinLength > 0 )
    FlagScore (graph);

  //-- Uniqueness requirements
  if ( OPT_MinUnique > 0 )
    FlagUNIQ (graph);

  //-- Query-based LIS
  if ( OPT_QLIS )
    FlagQLIS (graph);

  //-- Reference-based LIS
  if ( OPT_RLIS )
    FlagRLIS (graph);

  //-- Global LIS
  if ( OPT_GLIS )
    FlagGLIS (graph);


  //-- Output the filtered delta file
  PrintDelta (graph);

  return EXIT_SUCCESS;
}




//------------------------------------------------------------- BuildEdge ----//
void BuildEdge (DeltaEdge_t & edge, const DeltaRecord_t & rec)
{
  stringstream ss;
  vector<long int>::const_iterator di;
  DeltaEdgelet_t * p;

  vector<DeltaAlignment_t>::const_iterator i;
  for ( i = rec . aligns . begin( ); i != rec . aligns . end( ); ++ i )
    {
      //-- Set the edgelet
      p = new DeltaEdgelet_t( );

      p -> idy = i -> idy / 100.0;
      p -> sim = i -> sim / 100.0;
      p -> stp = i -> stp / 100.0;

      p -> idyc = i -> idyc;
      p -> simc = i -> simc;
      p -> stpc = i -> stpc;

      p -> loR = i -> sR;
      p -> hiR = i -> eR;
      p -> loQ = i -> sQ;
      p -> hiQ = i -> eQ;

      p -> dirR = p -> hiR < p -> loR ? REVERSE_DIR : FORWARD_DIR;
      p -> dirQ = p -> hiQ < p -> loQ ? REVERSE_DIR : FORWARD_DIR;

      //-- Get the delta information
      for ( di = i -> deltas . begin( ); di != i -> deltas . end( ); ++ di )
        {
          ss << *di << '\n';
          p -> delta . append (ss . str( ));
          ss . str ("");
        }

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


      //-- Build the edge
      BuildEdge (*dep, dr . getRecord( ));
      dep -> refnode -> edges . push_back (dep);
      dep -> qrynode -> edges . push_back (dep);
    }
  dr . close ( );
}




//-------------------------------------------------------------- FlagGLIS ----//
void FlagGLIS (DeltaGraph_t & graph)
{
  LIS_t * lis = NULL;
  long int lis_size = 0;
  long int i, j, n, best;
  long int olap, olapQ, olapR, len, lenQ, lenR, score;

  vector<DeltaEdgelet_t *> edgelets;

  map<string, DeltaNode_t>::const_iterator mi;
  vector<DeltaEdge_t *>::const_iterator ei;
  vector<DeltaEdgelet_t *>::iterator eli;


  //-- For each reference sequence
  for ( mi = graph.refnodes.begin( ); mi != graph.refnodes.end( ); ++ mi )
    {
      //-- For each query aligning to this reference
      for ( ei  = (mi -> second) . edges . begin( );
            ei != (mi -> second) . edges . end( ); ++ ei )
        {
          //-- Collect all the good edgelets
          edgelets . clear( );
          for ( eli  = (*ei) -> edgelets . begin( );
                eli != (*ei) -> edgelets . end( ); ++ eli )
            if ( (*eli) -> isGOOD )
              {
                edgelets . push_back (*eli);

                //-- Fix the coordinates to make global LIS work
                if ( (*eli) -> dirR == (*eli) -> dirQ )
                  {
                    (*eli) -> dirQ = FORWARD_DIR;
                  }
                else
                  {
                    if ( (*eli) -> dirQ == REVERSE_DIR )
                      Swap ((*eli) -> loQ, (*eli) -> hiQ);
                    (*eli) -> loQ = RevC ((*eli) -> loQ, (*ei)->qrynode->len);
                    (*eli) -> hiQ = RevC ((*eli) -> hiQ, (*ei)->qrynode->len);
                    (*eli) -> dirQ = REVERSE_DIR;
                  }
              }

          //-- Resize
          n = edgelets . size( );
          if ( n > lis_size )
            {
              lis = (LIS_t *) Safe_realloc (lis, sizeof (LIS_t) * n);
              lis_size = n;
            }

          //-- Sort by lo query coord
          sort (edgelets . begin( ), edgelets . end( ), EdgeletQCmp_t( ));

          //-- Dynamic
          for ( i = 0; i < n; ++ i )
            {
              lis[i] . a = edgelets[i];
              lis[i] . score = lis[i] . a -> hiQ - lis[i] . a -> loQ + 1;
              lis[i] . from = -1;

              for ( j = 0; j < i; ++ j )
                {
                  if ( lis[i] . a -> dirQ != lis[j] . a -> dirQ )
                    continue;
                  
                  lenR = lis[i] . a -> hiR - lis[i] . a -> loR + 1;
                  lenQ = lis[i] . a -> hiQ - lis[i] . a -> loQ + 1;
                  len = lenR > lenQ ? lenQ : lenR;
                  
                  olapR = lis[j] . a -> hiR - lis[i] . a -> loR + 1;
                  olapQ = lis[j] . a -> hiQ - lis[i] . a -> loQ + 1;
                  olap = olapR > olapQ ? olapR : olapQ;
                  if ( olap < 0 )
                    olap = 0;

                  score =
                    ScoreGlobal
                    (lis[j] . score, len, olap, lis[i] . a -> idy);
                  if ( score > lis[i] . score )
                    {
                      lis[i] . from = j;
                      lis[i] . score = score;
                    }
                }
            }

          //-- Backtrack and flag the GLIS edgelets
          best = 0;
          for ( i = 1; i < n; ++ i )
            if ( lis[i] . score > lis[best] . score )
              best = i;

          for ( i = best; i >= 0  &&  i < n; i = lis[i] . from )
            lis[i] . a -> isGLIS = true;

          //-- Flag the edgelets not in the GLIS
          for ( eli = edgelets . begin( ); eli != edgelets . end( ); ++ eli )
            {
              if ( ! (*eli) -> isGLIS )
                (*eli) -> isGOOD = false;

              //-- Bring the coordinates back to normal
              if ( (*eli) -> dirQ == FORWARD_DIR )
                {
                  (*eli) -> dirQ = (*eli) -> dirR;
                }
              else
                {
                  if ( (*eli) -> dirR == FORWARD_DIR )
                    Swap ((*eli) -> loQ, (*eli) -> hiQ);
                  (*eli) -> loQ = RevC ((*eli) -> loQ, (*ei)->qrynode->len);
                  (*eli) -> hiQ = RevC ((*eli) -> hiQ, (*ei)->qrynode->len);
                  (*eli) -> dirQ =
                    (*eli) -> dirR == FORWARD_DIR ? REVERSE_DIR : FORWARD_DIR;
                }
            }
        }
    }

  free (lis);
}





//------------------------------------------------------------- FlagScore ----//
void FlagScore (DeltaGraph_t & graph)
{
  map<string, DeltaNode_t>::const_iterator mi;
  vector<DeltaEdge_t *>::const_iterator ei;
  vector<DeltaEdgelet_t *>::iterator eli;

  for ( mi = graph.refnodes.begin( ); mi != graph.refnodes.end( ); ++ mi )
    for ( ei  = (mi -> second) . edges . begin( );
          ei != (mi -> second) . edges . end( ); ++ ei )
      for ( eli  = (*ei) -> edgelets . begin( );
            eli != (*ei) -> edgelets . end( ); ++ eli )
        if ( (*eli) -> isGOOD )
          {
            //-- Flag low identities
            if ( (*eli)->idy * 100.0 < OPT_MinIdentity )
              (*eli) -> isGOOD = false;

            //-- Flag small lengths
            if ( (*eli)->hiR - (*eli)->loR + 1 < (unsigned long)OPT_MinLength ||
                 (*eli)->hiQ - (*eli)->loQ + 1 < (unsigned long)OPT_MinLength )
              (*eli) -> isGOOD = false;
          }
}




//-------------------------------------------------------------- FlagQLIS ----//
void FlagQLIS (DeltaGraph_t & graph)
{
  LIS_t * lis = NULL;
  long int lis_size = 0;
  long int i, j, n, best;
  long int olap, leni, lenj, score;

  vector<DeltaEdgelet_t *> edgelets;

  map<string, DeltaNode_t>::const_iterator mi;
  vector<DeltaEdge_t *>::const_iterator ei;
  vector<DeltaEdgelet_t *>::iterator eli;


  //-- For each query sequence
  for ( mi = graph.qrynodes.begin( ); mi != graph.qrynodes.end( ); ++ mi )
    {
      //-- Collect all the good edgelets
      edgelets . clear( );
      for ( ei  = (mi -> second) . edges . begin( );
            ei != (mi -> second) . edges . end( ); ++ ei )
        for ( eli  = (*ei) -> edgelets . begin( );
              eli != (*ei) -> edgelets . end( ); ++ eli )
          if ( (*eli) -> isGOOD )
            edgelets . push_back (*eli);

      //-- Resize
      n = edgelets . size( );
      if ( n > lis_size )
        {
          lis = (LIS_t *) Safe_realloc (lis, sizeof (LIS_t) * n);
          lis_size = n;
        }

      //-- Sort by lo query coord
      sort (edgelets . begin( ), edgelets . end( ), EdgeletQCmp_t( ));

      //-- Dynamic
      for ( i = 0; i < n; ++ i )
        {
          lis[i] . a = edgelets[i];
          lis[i] . score = lis[i] . a -> hiQ - lis[i] . a -> loQ + 1;
          lis[i] . from = -1;

          for ( j = 0; j < i; ++ j )
            {
              leni = lis[i] . a -> hiQ - lis[i] . a -> loQ + 1;
              lenj = lis[j] . a -> hiQ - lis[j] . a -> loQ + 1;
              olap = lis[j] . a -> hiQ - lis[i] . a -> loQ + 1;
              if ( olap < 0 )
                olap = 0;
              
              score = ScoreLocal
                (lis[j] . score, leni, lenj, olap, lis[i] . a -> idy);
              if ( score > lis[i] . score )
                {
                  lis[i] . from = j;
                  lis[i] . score = score;
                }
            }
        }

      //-- Backtrack and flag the QLIS edgelets
      best = 0;
      for ( i = 1; i < n; ++ i )
        if ( lis[i] . score > lis[best] . score )
          best = i;

      for ( i = best; i >= 0  &&  i < n; i = lis[i] . from )
        lis[i] . a -> isQLIS = true;

      //-- Flag the edgelets not in the QLIS
      for ( eli = edgelets . begin( ); eli != edgelets . end( ); ++ eli )
        if ( ! (*eli) -> isQLIS )
          (*eli) -> isGOOD = false;
    }

  free (lis);
}




//-------------------------------------------------------------- FlagRLIS ----//
void FlagRLIS (DeltaGraph_t & graph)
{
  LIS_t * lis = NULL;
  long int lis_size = 0;
  long int i, j, n, best;
  long int olap, leni, lenj, score;

  vector<DeltaEdgelet_t *> edgelets;

  map<string, DeltaNode_t>::const_iterator mi;
  vector<DeltaEdge_t *>::const_iterator ei;
  vector<DeltaEdgelet_t *>::iterator eli;


  //-- For each reference sequence
  for ( mi = graph.refnodes.begin( ); mi != graph.refnodes.end( ); ++ mi )
    {
      //-- Collect all the good edgelets
      edgelets . clear( );
      for ( ei  = (mi -> second) . edges . begin( );
            ei != (mi -> second) . edges . end( ); ++ ei )
        for ( eli  = (*ei) -> edgelets . begin( );
              eli != (*ei) -> edgelets . end( ); ++ eli )
          if ( (*eli) -> isGOOD )
            edgelets . push_back (*eli);

      //-- Resize
      n = edgelets . size( );
      if ( n > lis_size )
        {
          lis = (LIS_t *) Safe_realloc (lis, sizeof (LIS_t) * n);
          lis_size = n;
        }

      //-- Sort by lo reference coord
      sort (edgelets . begin( ), edgelets . end( ), EdgeletRCmp_t( ));

      //-- Dynamic
      for ( i = 0; i < n; ++ i )
        {
          lis[i] . a = edgelets[i];
          lis[i] . score = lis[i] . a -> hiR - lis[i] . a -> loR + 1;
          lis[i] . from = -1;

          for ( j = 0; j < i; ++ j )
            {
              leni = lis[i] . a -> hiR - lis[i] . a -> loR + 1;
              lenj = lis[j] . a -> hiR - lis[j] . a -> loR + 1;
              olap = lis[j] . a -> hiR - lis[i] . a -> loR + 1;
              if ( olap < 0 )
                olap = 0;

              score = ScoreLocal
                (lis[j] . score, leni, lenj, olap, lis[i] . a -> idy);
              if ( score > lis[i] . score )
                {
                  lis[i] . from = j;
                  lis[i] . score = score;
                }
            }
        }

      //-- Backtrack and flag the RLIS edgelets
      best = 0;
      for ( i = 1; i < n; ++ i )
        if ( lis[i] . score > lis[best] . score )
          best = i;

      for ( i = best; i >= 0  &&  i < n; i = lis[i] . from )
        lis[i] . a -> isRLIS = true;

      //-- Flag the edgelets not in the RLIS
      for ( eli = edgelets . begin( ); eli != edgelets . end( ); ++ eli )
        if ( ! (*eli) -> isRLIS )
          (*eli) -> isGOOD = false;
    }

  free (lis);
}




//-------------------------------------------------------------- FlagUNIQ ----//
void FlagUNIQ (DeltaGraph_t & graph)
{
  unsigned long int i, uniq, len;

  vector<DeltaEdgelet_t *> edgelets;

  map<string, DeltaNode_t>::const_iterator mi;
  vector<DeltaEdge_t *>::const_iterator ei;
  vector<DeltaEdgelet_t *>::iterator eli;


  //-- For each reference sequence
  unsigned long int ref_size = 0;
  unsigned long int ref_len = 0;
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

      //-- Collect all the good edgelets
      edgelets . clear( );
      for ( ei  = (mi -> second) . edges . begin( );
            ei != (mi -> second) . edges . end( ); ++ ei )
        for ( eli  = (*ei) -> edgelets . begin( );
              eli != (*ei) -> edgelets . end( ); ++ eli )
          if ( (*eli) -> isGOOD )
            {
              edgelets . push_back (*eli);

              //-- Add to the reference coverage
              for ( i = (*eli) -> loR; i <= (*eli) -> hiR; i ++ )
                if ( ref_cov[i] < UCHAR_MAX )
                  ref_cov[i] ++;
            }

      //-- Calculate the uniqueness of each edgelet
      for ( eli = edgelets . begin( ); eli != edgelets . end( ); ++ eli )
        {
          uniq = 0;
          len = (*eli) -> hiR - (*eli) -> loR + 1;
          for ( i = (*eli) -> loR; i <= (*eli) -> hiR; i ++ )
            if ( ref_cov[i] == 1 )
              uniq ++;

          //-- Flag low reference uniqueness
          if ( (float)uniq / (float)len * 100.0 < OPT_MinUnique )
            (*eli) -> isGOOD = false;
        }
    }
  free (ref_cov);


  //-- For each query sequence
  unsigned long int qry_size = 0;
  unsigned long int qry_len = 0;
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

      //-- Collect all the good edgelets
      edgelets . clear( );
      for ( ei  = (mi -> second) . edges . begin( );
            ei != (mi -> second) . edges . end( ); ++ ei )
        for ( eli  = (*ei) -> edgelets . begin( );
              eli != (*ei) -> edgelets . end( ); ++ eli )
          if ( (*eli) -> isGOOD )
            {
              edgelets . push_back (*eli);

              //-- Add to the query coverage
              for ( i = (*eli) -> loQ; i <= (*eli) -> hiQ; i ++ )
                if ( qry_cov[i] < UCHAR_MAX )
                  qry_cov[i] ++;
            }

      //-- Calculate the uniqueness of each edgelet
      for ( eli = edgelets . begin( ); eli != edgelets . end( ); ++ eli )
        {
          uniq = 0;
          len = (*eli) -> hiQ - (*eli) -> loQ + 1;
          for ( i = (*eli) -> loQ; i <= (*eli) -> hiQ; i ++ )
            if ( qry_cov[i] == 1 )
              uniq ++;
          
          //-- Flag low query uniqueness
          if ( (float)uniq / (float)len * 100.0 < OPT_MinUnique )
            (*eli) -> isGOOD = false;
        }
    }
  free (qry_cov);
}




//------------------------------------------------------------ PrintDelta ----//
void PrintDelta (const DeltaGraph_t & graph)
{
  bool header;
  unsigned long int s1, e1, s2, e2;

  map<string, DeltaNode_t>::const_iterator mi;
  vector<DeltaEdge_t *>::const_iterator ei;
  vector<DeltaEdgelet_t *>::const_iterator eli;

  //-- Print the file header
  cout
    << graph . refpath << ' ' << graph . qrypath << endl
    << (graph.datatype == PROMER_DATA ? PROMER_STRING : NUCMER_STRING) << endl;

  for ( mi = graph.qrynodes.begin( ); mi != graph.qrynodes.end( ); ++ mi )
    {
      for ( ei  = (mi -> second) . edges . begin( );
            ei != (mi -> second) . edges . end( ); ++ ei )
        {
          header = false;

          for ( eli  = (*ei) -> edgelets . begin( );
                eli != (*ei) -> edgelets . end( ); ++ eli )
            {
              if ( ! (*eli) -> isGOOD )
                continue;

              //-- Print the sequence header
              if ( ! header )
                {
                  cout
                    << '>'
                    << *((*ei) -> refnode -> id) << ' '
                    << *((*ei) -> qrynode -> id) << ' '
                    << (*ei) -> refnode -> len << ' '
                    << (*ei) -> qrynode -> len << endl;
                  header = true;
                }

              //-- Print the alignment
              s1 = (*eli) -> loR;
              e1 = (*eli) -> hiR;
              s2 = (*eli) -> loQ;
              e2 = (*eli) -> hiQ;
              if ( (*eli) -> dirR == REVERSE_DIR )
                Swap (s1, e1);
              if ( (*eli) -> dirQ == REVERSE_DIR )
                Swap (s2, e2);

              cout
                << s1 << ' ' << e1 << ' ' << s2 << ' ' << e2 << ' '
                << (*eli) -> idyc << ' '
                << (*eli) -> simc << ' '
                << (*eli) -> stpc << endl
                << (*eli) -> delta;
            }
        }
    }
}




//------------------------------------------------------------- ParseArgs ----//
void ParseArgs (int argc, char ** argv)
{
  int ch, errflg = 0;
  optarg = NULL;
  
  while ( !errflg  &&
          ((ch = getopt (argc, argv, "ghi:l:o:qru:")) != EOF) )
    switch (ch)
      {
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

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
#include <iostream>
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
bool    OPT_ShowConflict  = true;       // -C option
bool    OPT_ShowIndels    = true;       // -I option
bool    OPT_PrintTabular  = false;      // -T option
bool    OPT_PrintHeader   = true;       // -H option
bool    OPT_SelectAligns  = false;      // -S option
bool    OPT_SNPs          = true;       // -V option

bool    OPT_filtNs        = false;       // -N option
unsigned int     OPT_minNlen	      = 10;         // -L option
float   OPT_minNper       = 0.70;        // -% option

int     OPT_Context       = 0;          // -x option

set<string> OPT_Aligns;                 // -S option

//============================================================= Constants ====//
const char  INDEL_CHAR = '.';
const char SEQEND_CHAR = '-';



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


struct VCF_R_Sort
{
  bool operator() (const VCFe_t * a, const VCFe_t * b)
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


struct VCF_Q_Sort
{
  bool operator() (const VCFe_t * a, const VCFe_t * b)
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
inline long RevC (long coord, long len)
{
  return len - coord + 1;
}


//------------------------------------------------------------------ Norm ----//
inline long Norm (long c, long l, int f, AlignmentType_t d)
{
  long retval = (d == PROMER_DATA ? c * 3 - (3 - abs(f)) : c);
  if ( f < 0 ) retval = RevC (retval, l);
  return retval;
}


//------------------------------------------------------------------ Swap ----//
inline void Swap (long & a, long & b)
{
  long t = a; a = b; b = t;
}

//----------------------------------------------------------- toStdNucleo ----//
char toStdNucleo (char x);

//----------------------------------------------------------- isStdNucleo ----//
bool isStdNucleo(char x);

//------------------------------------------------------------- CheckSNPs ----//
void CheckSNPs (DeltaGraph_t & graph);

//------------------------------------------------------------- CheckVCFs ----//
void CheckVCFs (DeltaGraph_t & graph);

//-------------------------------------------------------------- FindSNPs ----//
void FindSNPs (DeltaGraph_t & graph);

//-------------------------------------------------------- getIndelSNPvcf ----//

VCFe_t getIndelSNPvcf(DeltaGraph_t & graph,
                      const vector<long>& deltas,
                      vector<long>::iterator& dit, 
                      long delta,
                      long& rpos, long& qpos,
                      vector< DeltaEdge_t * >::iterator ei,
                      vector< DeltaEdgelet_t * >::iterator li,
                      char ** R, char ** Q,
                      int ri, int qi,  
                      int frameR, int frameQ,
                      long lenR, long lenQ,
                      bool add2VCFs
                      );

//-------------------------------------------------------------- FindVCFs ----//
void FindVCFs (DeltaGraph_t & graph);

//------------------------------------------------------------ PrintHuman ----//
void PrintHuman (const vector<const SNP_t *> & snps,
                 const DeltaGraph_t & graph);


//------------------------------------------------------------- PrintVCFs ----//
void PrintVCFs (const vector<const VCFe_t *> & vcfs, int argc, char ** argv, DeltaGraph_t & graph);

//---------------------------------------------------------- PrintTabular ----//
void PrintTabular (const vector<const SNP_t *> & snps,
                   const DeltaGraph_t & graph);


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
  DeltaGraph_t graph;
  
  //-- Command line parsing
  ParseArgs (argc, argv);

  //-- Select alignments
  if ( OPT_SelectAligns )
    SelectAligns ( );

  //-- Build the alignment graph from the delta file
  graph . build (OPT_AlignName, true);

  //-- Read sequences
  graph . loadSequences ( );

  if(OPT_SNPs){
    vector<const SNP_t *> snps;
    
    //-- Locate the SNPs
    FindSNPs (graph);

    //-- Check for ambiguous alignment regions
    CheckSNPs (graph);


    //-- Collect and sort the SNPs
    map<string, DeltaNode_t>::iterator mi;
    vector<DeltaEdge_t *>::iterator ei;
    vector<DeltaEdgelet_t *>::iterator li;
    vector<SNP_t *>::iterator si;
    for ( mi = graph.refnodes.begin( ); mi != graph.refnodes.end( ); ++ mi )
      for ( ei = mi->second.edges.begin( ); ei != mi->second.edges.end( ); ++ ei )
        for ( li = (*ei)->edgelets.begin( ); li != (*ei)->edgelets.end( ); ++ li )
          for ( si = (*li)->snps.begin( ); si != (*li)->snps.end( ); ++ si )
            if ( (OPT_ShowConflict ||
                  ((*si)->conR == 0 && (*si)->conQ == 0))
                &&
                (OPT_ShowIndels ||
                  ((*si)->cR != INDEL_CHAR && (*si)->cQ != INDEL_CHAR)) )
              snps . push_back (*si);

    if ( OPT_SortReference )
      sort (snps . begin( ), snps . end( ), SNP_R_Sort( ));
    else
      sort (snps . begin( ), snps . end( ), SNP_Q_Sort( ));
      
      //-- Output data to stdout
    if ( OPT_PrintTabular )
      PrintTabular (snps, graph);
    else
      PrintHuman (snps, graph); 
      
  }else{
    vector<const VCFe_t *> vcfs;
    
    if(graph.datatype == PROMER_DATA){
      cerr<<"The VCF format is only designed for\n"
          <<"DNA variants. Please provide a delta\n"
          <<"file generated with nucmer\n";
      exit(1);
    }
    
    //-- Locate the VCFevents
    FindVCFs (graph);
    //-- Check for ambiguous alignment regions
    CheckVCFs (graph);
    //-- Collect and sort the SNPs
    map<string, DeltaNode_t>::iterator mi;
    vector<DeltaEdge_t *>::iterator ei;
    vector<DeltaEdgelet_t *>::iterator li;
    vector<VCFe_t *>::iterator vi;
    for ( mi = graph.refnodes.begin( ); mi != graph.refnodes.end( ); ++ mi )
      for ( ei = mi->second.edges.begin( ); ei != mi->second.edges.end( ); ++ ei )
        for ( li = (*ei)->edgelets.begin( ); li != (*ei)->edgelets.end( ); ++ li )
          for ( vi = (*li)->vcfs.begin( ); vi != (*li)->vcfs.end( ); ++ vi )
            if ( (OPT_ShowConflict ||
                  ((*vi)->conR == 0 && (*vi)->conQ == 0))
                &&
                (OPT_ShowIndels ||
                  ((*vi)->type != "INDEL")) )
              vcfs . push_back (*vi);

    if ( OPT_SortReference )
      sort (vcfs . begin( ), vcfs . end( ), VCF_R_Sort( ));
    else
      sort (vcfs . begin( ), vcfs . end( ), VCF_Q_Sort( ));

    PrintVCFs(vcfs,argc,argv, graph);
  
	}
  
  return EXIT_SUCCESS;
}

//----------------------------------------------------------- toStdNucleo ----//

char toStdNucleo (char x){
	switch (toupper(x)){
    case 'R':
      return('A');
      break;
    
    case 'Y':
      return('C');
      break;
      
    case 'M':
      return('A');
      break;
      
    case 'K':
      return('G');
      break;
      
    case 'S':
      return('C');
      break;
      
    case 'W':
      return('A');
      break;
      
    case 'B':
      return('C');
      break;
      
    case 'D':
      return('A');
      break;
      
    case 'H':
      return('A');
      break;
      
    case 'V':
      return('A');
      break;
      
    default:
      return(x);
      
  }
      
}

//----------------------------------------------------------- isStdNucleo ----//
bool isStdNucleo(char x){
  switch(toupper(x)){
    case 'A':
      return true;
      
    case 'T':
      return true;
      
    case 'C':
      return true;
      
    case 'G':
      return true;
      
    case 'N':
      return true;
      
    default:
      return false;
    
  }
}

//------------------------------------------------------------- CheckSNPs ----//
void CheckSNPs (DeltaGraph_t & graph)
{
  map<string, DeltaNode_t>::const_iterator mi;
  vector<DeltaEdge_t *>::const_iterator ei;
  vector<DeltaEdgelet_t *>::iterator eli;
  vector<SNP_t *>::iterator si;
  long i;

  //-- For each reference sequence
  long ref_size = 0;
  long ref_len = 0;
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
            (*si) -> conR = ref_cov [(*si)->pR] - 1;
    }
  free (ref_cov);


  //-- For each query sequence
  long qry_size = 0;
  long qry_len = 0;
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
            (*si) -> conQ = qry_cov [(*si)->pQ] - 1;
    }
  free (qry_cov);
}

//------------------------------------------------------------- CheckVCFs ----//
void CheckVCFs (DeltaGraph_t & graph)
{
  map<string, DeltaNode_t>::const_iterator mi;
  vector<DeltaEdge_t *>::const_iterator ei;
  vector<DeltaEdgelet_t *>::iterator eli;
  vector<VCFe_t *>::iterator si;
  long i;

  //-- For each reference sequence
  long ref_size = 0;
  long ref_len = 0;
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
          for ( si = (*eli)->vcfs.begin( ); si != (*eli)->vcfs.end( ); ++ si )
            if((*si)->pR >=1 && (*si)->pR <= ref_len)
				(*si) -> conR = ref_cov [(*si)->pR] - 1;
			else
				cerr<<"Out of range:\n"
				    <<"\t"<<*((*si)->ep->refnode->id)<<"\n"				
				    <<"\t"<<(*si)->pR<<"\n";
    }
  free (ref_cov);


  //-- For each query sequence
  long qry_size = 0;
  long qry_len = 0;
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
          for ( si = (*eli)->vcfs.begin( ); si != (*eli)->vcfs.end( ); ++ si )
            if((*si)->pQ >=1 && (*si)->pQ <= qry_len)
				(*si) -> conQ = qry_cov [(*si)->pQ] - 1;
			else
				cerr<<"Out of range:\n"
				    <<"\t"<<*((*si)->ep->qrynode->id)<<"\n"				
				    <<"\t"<<(*si)->pQ<<"\n";
    }
  free (qry_cov);
}


//-------------------------------------------------------------- FindSNPs ----//
void FindSNPs (DeltaGraph_t & graph)
{
  map<string, DeltaNode_t>::iterator mi;
  vector<DeltaEdge_t *>::iterator ei;
  vector<DeltaEdgelet_t *>::iterator li;
  vector<SNP_t *>::iterator si, psi, nsi;

  //-- For each alignment, identify the SNPs
  for ( mi = graph.refnodes.begin( ); mi != graph.refnodes.end( ); ++ mi )
    for ( ei = mi->second.edges.begin( ); ei != mi->second.edges.end( ); ++ ei )
      {
        SNP_t * snp;
        int ri, qi;
        char * R[] = {(*ei)->refnode->seq, NULL, NULL, NULL, NULL, NULL, NULL};
        char * Q[] = {(*ei)->qrynode->seq, NULL, NULL, NULL, NULL, NULL, NULL};

        long i;
        long lenR = (*ei) -> refnode -> len;
        long lenQ = (*ei) -> qrynode -> len;

        for (li = (*ei)->edgelets.begin( ); li != (*ei)->edgelets.end( ); ++ li)
          {
            long delta;
            int frameR, frameQ, sign;
            long sR, eR, sQ, eQ;
            long rpos, qpos, remain;
            long rctx, qctx;
            long alenR = lenR;
            long alenQ = lenQ;

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

            istringstream ss;
            ss . str ((*li)->delta);

            while ( ss >> delta && delta != 0 )
              {
                sign = delta > 0 ? 1 : -1;
                delta = labs (delta);

                //-- For all SNPs before the next indel
                for ( i = 1; i < delta; i ++ )
                  if ( R [ri] [rpos ++] != Q [qi] [qpos ++] )
                    {
                      if ( graph.datatype == NUCMER_DATA &&
                           CompareIUPAC (R [ri][rpos-1], Q [qi][qpos-1]) )
                        continue;

                      snp = new SNP_t;
                      snp -> ep = *ei;
                      snp -> lp = *li;
                      snp -> pR = Norm (rpos-1, lenR, frameR, graph.datatype);
                      snp -> pQ = Norm (qpos-1, lenQ, frameQ, graph.datatype);
                      snp -> cR = toupper (R [ri] [rpos-1]);
                      snp -> cQ = toupper (Q [qi] [qpos-1]);

                      for ( rctx = rpos - OPT_Context - 1;
                            rctx < rpos + OPT_Context; rctx ++ )
                        if ( rctx < 1  ||  rctx > alenR )
                          snp -> ctxR . push_back (SEQEND_CHAR);
                        else if ( rctx == rpos - 1 )
                          snp -> ctxR . push_back (snp -> cR);
                        else
                          snp -> ctxR . push_back (toupper (R [ri] [rctx]));
                          
                      for ( qctx = qpos - OPT_Context - 1;
                            qctx < qpos + OPT_Context; qctx ++ )
                        if ( qctx < 1  ||  qctx > alenQ )
                          snp -> ctxQ . push_back (SEQEND_CHAR);
                        else if ( qctx == qpos - 1 )
                          snp -> ctxQ . push_back (snp -> cQ);
                        else
                          snp -> ctxQ . push_back (toupper (Q [qi] [qctx]));

                      (*li) -> snps . push_back (snp);
                    }

                //-- For the indel
                snp = new SNP_t;
                snp -> ep = *ei;
                snp -> lp = *li;

                for ( rctx = rpos - OPT_Context; rctx < rpos; rctx ++ )
                  if ( rctx < 1 )
                    snp -> ctxR . push_back (SEQEND_CHAR);
                  else
                    snp -> ctxR . push_back (toupper (R [ri] [rctx]));
                
                for ( qctx = qpos - OPT_Context; qctx < qpos; qctx ++ )
                  if ( qctx < 1 )
                    snp -> ctxQ . push_back (SEQEND_CHAR);
                  else
                    snp -> ctxQ . push_back (toupper (Q [qi] [qctx]));

                if ( sign > 0 )
                  {
                    snp -> pR = Norm (rpos, lenR, frameR, graph.datatype);
                    if ( frameQ > 0 )
                      snp -> pQ = Norm (qpos - 1, lenQ, frameQ, graph.datatype);
                    else
                      snp -> pQ = Norm (qpos, lenQ, frameQ, graph.datatype);

                    snp -> cR = toupper (R [ri] [rpos ++]);
                    snp -> cQ = INDEL_CHAR;

                    remain -= i;
                    rctx ++;
                  }
                else
                  {
                    snp -> pQ = Norm (qpos, lenQ, frameQ, graph.datatype);
                    if ( frameR > 0 )
                      snp -> pR = Norm (rpos - 1, lenR, frameR, graph.datatype);
                    else
                      snp -> pR = Norm (rpos, lenR, frameR, graph.datatype);

                    snp -> cR = INDEL_CHAR;
                    snp -> cQ = toupper (Q [qi] [qpos ++]);

                    remain -= i - 1;
                    qctx ++;
                  }

                snp -> ctxR . push_back (snp -> cR);
                for ( ; rctx < rpos + OPT_Context; rctx ++ )
                  if ( rctx > alenR )
                    snp -> ctxR . push_back (SEQEND_CHAR);
                  else
                    snp -> ctxR . push_back (toupper (R [ri] [rctx]));
                
                snp -> ctxQ . push_back (snp -> cQ);
                for ( ; qctx < qpos + OPT_Context; qctx ++ )
                  if ( qctx > alenQ )
                    snp -> ctxQ . push_back (SEQEND_CHAR);
                  else
                    snp -> ctxQ . push_back (toupper (Q [qi] [qctx]));
                
                (*li) -> snps . push_back (snp);
              }

            //-- For all SNPs after the final indel
            for ( i = 0; i < remain; i ++ )
              if ( R [ri] [rpos ++] != Q [qi] [qpos ++] )
                {
                  if ( graph.datatype == NUCMER_DATA &&
                       CompareIUPAC (R [ri][rpos-1], Q [qi][qpos-1]) )
                    continue;

                  snp = new SNP_t;
                  snp -> ep = *ei;
                  snp -> lp = *li;
                  snp -> pR = Norm (rpos-1, lenR, frameR, graph.datatype);
                  snp -> pQ = Norm (qpos-1, lenQ, frameQ, graph.datatype);
                  snp -> cR = toupper (R [ri] [rpos-1]);
                  snp -> cQ = toupper (Q [qi] [qpos-1]);

                  for ( rctx = rpos - OPT_Context - 1;
                        rctx < rpos + OPT_Context; rctx ++ )
                    if ( rctx < 1  ||  rctx > alenR )
                      snp -> ctxR . push_back (SEQEND_CHAR);
                    else if ( rctx == rpos - 1 )
                      snp -> ctxR . push_back (snp -> cR);
                    else
                      snp -> ctxR . push_back (toupper (R [ri] [rctx]));
                  
                  for ( qctx = qpos - OPT_Context - 1;
                        qctx < qpos + OPT_Context; qctx ++ )
                    if ( qctx < 1  ||  qctx > alenQ )
                      snp -> ctxQ . push_back (SEQEND_CHAR);
                    else if ( qctx == qpos - 1 )
                      snp -> ctxQ . push_back (snp -> cQ);
                    else
                      snp -> ctxQ . push_back (toupper (Q [qi] [qctx]));

                  (*li) -> snps . push_back (snp);
                }


            //-- Sort SNPs and calculate distances
            if ( OPT_SortReference )
              {
                sort ((*li)->snps.begin( ), (*li)->snps.end( ), SNP_R_Sort( ));

                for ( si = (*li)->snps.begin(); si != (*li)->snps.end(); ++ si )
                  {
                    psi = si - 1;
                    nsi = si + 1;

                    (*si) -> buff = 1 +
                      ((*si)->pR - (*li)->loR < (*li)->hiR - (*si)->pR ?
                       (*si)->pR - (*li)->loR : (*li)->hiR - (*si)->pR);

                    if ( psi >= (*li) -> snps . begin( )  &&
                         (*si)->pR - (*psi)->pR < (*si)->buff )
                      (*si) -> buff = (*si)->pR - (*psi)->pR;
                    
                    if ( nsi < (*li) -> snps . end( )  &&
                         (*nsi)->pR - (*si)->pR < (*si)->buff )
                      (*si) -> buff = (*nsi)->pR - (*si)->pR;
                  }
              }
            else
              {
                sort ((*li)->snps.begin( ), (*li)->snps.end( ), SNP_Q_Sort( ));

                for ( si = (*li)->snps.begin(); si != (*li)->snps.end(); ++ si )
                  {
                    psi = si - 1;
                    nsi = si + 1;
 
                    (*si) -> buff = 1 +
                      ((*si)->pQ - (*li)->loQ < (*li)->hiQ - (*si)->pQ ?
                       (*si)->pQ - (*li)->loQ : (*li)->hiQ - (*si)->pQ);

                    if ( psi >= (*li) -> snps . begin( )  &&
                         (*si)->pQ - (*psi)->pQ < (*si)->buff )
                      (*si) -> buff = (*si)->pQ - (*psi)->pQ;
                    
                    if ( nsi < (*li) -> snps . end( )  &&
                         (*nsi)->pQ - (*si)->pQ < (*si)->buff )
                      (*si) -> buff = (*nsi)->pQ - (*si)->pQ;
                  }
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


//-------------------------------------------------------------- FindVCFs ----//
void FindVCFs (DeltaGraph_t & graph)
{
  map< string, DeltaNode_t >::iterator mi;
  vector< VCFe_t * >::iterator vi, pvi, nvi;
  vector< DeltaEdge_t * >::iterator ei;
  vector< DeltaEdgelet_t * >::iterator li;
  
  //-- For each alignment, identify the SNPs
  for ( mi = graph.refnodes.begin( ); mi != graph.refnodes.end( ); ++ mi ){
    for ( ei = mi->second.edges.begin( ); ei != mi->second.edges.end( ); ei ++ )
      {
        VCFe_t * vcf;
        char * R[] = {(*ei)->refnode->seq, NULL, NULL, NULL, NULL, NULL, NULL};
        char * Q[] = {(*ei)->qrynode->seq, NULL, NULL, NULL, NULL, NULL, NULL};
        int ri, qi;
        
        long i;
        long lenR = (*ei) -> refnode -> len;
        long lenQ = (*ei) -> qrynode -> len;

        for (li = (*ei)->edgelets.begin( ); li != (*ei)->edgelets.end( ); ++ li)
          {
            long delta;
            int frameR, frameQ;
            long sR, eR, sQ, eQ;
            long rpos, qpos, remain;
            long alenR = lenR;
            long alenQ = lenQ;

            //-- Only do the ones requested by user
            if ( OPT_SelectAligns )
              {
                ostringstream ss;
                set<string>::iterator vi;

                if ( (*li) -> dirR == FORWARD_DIR )
                  ss << (*li) -> loR << ' ' << (*li) -> hiR << ' ';
                else
                  ss << (*li) -> hiR << ' ' << (*li) -> loR << ' ';

                if ( (*li) -> dirQ == FORWARD_DIR )
                  ss << (*li) -> loQ << ' ' << (*li) -> hiQ << ' ';
                else
                  ss << (*li) -> hiQ << ' ' << (*li) -> loQ << ' ';

                ss << *((*ei)->refnode->id) << ' ' << *((*ei)->qrynode->id);

                vi = OPT_Aligns . find (ss .str( ));
                if ( vi == OPT_Aligns . end( ) )
                  continue;
                else
                  OPT_Aligns . erase (vi);
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

            istringstream ss;
            ss . str ((*li)->delta);
            
            vector<long> deltas;
            while(ss >> delta){
              deltas.push_back(delta);
            }
            
            vector<long>::iterator dit;
            for ( dit=deltas.begin();dit<deltas.end() && *dit != 0; dit++ )
                getIndelSNPvcf(graph,deltas, dit, *dit,
                                rpos, qpos, 
                                ei, li,
                                R, Q, 
                                ri, qi,
                                frameR, frameQ,
                                lenR, lenQ,
                                true);  
            
            //-- For all vcfs after the final indel
            remain = eR - rpos + 1;
            for ( i = 0; i < remain; i ++ )
              if ( R [ri] [rpos ++] != Q [qi] [qpos ++] )
                {
                  if ( graph.datatype == NUCMER_DATA &&
                           CompareIUPAC (R [ri][rpos-1], Q [qi][qpos-1]) )
                        continue;

                      vcf = new VCFe_t;
                      vcf -> ep = *ei;
                      vcf -> lp = *li;
                      vcf -> pR = Norm (rpos-1, lenR, frameR, graph.datatype);
                      vcf -> stdN = vcf->stdN && isStdNucleo(R [0] [vcf -> pR]);
                      vcf -> cR = toStdNucleo(toupper (R [0] [vcf -> pR])); //use letter in forward reference
                      if(frameR > 0){
                        vcf -> stdN = vcf->stdN && isStdNucleo(toupper( Q [qi] [qpos-1]));
                        vcf -> cQ = toStdNucleo(toupper (Q [qi] [qpos-1]));
                        vcf -> pQ = qpos-1;
                      }else{
                        vcf -> stdN = vcf->stdN && isStdNucleo(Complement( Q [qi] [qpos-1]));
                        vcf -> cQ = toStdNucleo(toupper (Complement(Q [qi] [qpos-1])));
                        vcf -> pQ = RevC(qpos-1,lenQ);
											}
                      vcf -> ctxR.push_back(vcf->cR);
                      vcf -> ctxQ.push_back(vcf->cQ);                      
                      vcf -> type = "SNP";
                      vcf -> infType = "SNP";
                      vcf -> inf_sQ = vcf -> pQ;
                      vcf -> inf_fQ = frameQ * frameR;
                      (*li) -> vcfs . push_back (vcf);
                }

            //-- Sort vcfs
            if ( OPT_SortReference )
              {
                sort ((*li)->vcfs.begin( ), (*li)->vcfs.end( ), VCF_R_Sort( ));
              }
            else
              {
                sort ((*li)->vcfs.begin( ), (*li)->vcfs.end( ), VCF_Q_Sort( ));
              }
          }

        //-- Clear up the seq
        for ( i = 1; i <= 6; i ++ )
          {
            free (R[i]);
            free (Q[i]);
          }
      }
  }

  if ( OPT_SelectAligns  &&  ! OPT_Aligns . empty( ) )
    {
      cerr << "ERROR: One or more alignments from stdin could not be found\n";
      exit (EXIT_FAILURE);
    }
}

//--------------------------------------------------------- getIndelSNPvcf ----//
VCFe_t getIndelSNPvcf(DeltaGraph_t & graph,
                      const vector<long>& deltas,
                      vector<long>::iterator& dit, 
                      long delta,
                      long& rpos, long& qpos,
                      vector< DeltaEdge_t * >::iterator ei,
                      vector< DeltaEdgelet_t * >::iterator li,
                      char ** R, char ** Q,
                      int ri, int qi,  
                      int frameR, int frameQ,
                      long lenR, long lenQ,
                      bool add2VCFs
                      ){
                VCFe_t* vcf;
                
                long rctx;
                int sign = *dit > 0 ? 1 : -1;
                delta = labs (delta);

                long i;
                //-- For all SNPs before the next indel
                for ( i = 1; i < delta; i ++ )
                  if ( R [ri] [rpos ++] != Q [qi] [qpos ++] )
                    {
                      if ( graph.datatype == NUCMER_DATA &&
                           CompareIUPAC (R [ri][rpos-1], Q [qi][qpos-1]) )
                        continue;

                      vcf = new VCFe_t;
                      vcf -> ep = *ei;
                      vcf -> lp = *li;
                      vcf -> pR = Norm (rpos-1, lenR, frameR, graph.datatype);
                      vcf -> stdN = vcf -> stdN && isStdNucleo(R [0] [vcf -> pR]);
                      vcf -> cR = toStdNucleo(toupper (R [0] [vcf -> pR])); //use letter in forward reference
                      if(frameR > 0){
                        vcf -> stdN = vcf -> stdN && isStdNucleo(Q [qi] [qpos-1]);
                        vcf -> cQ = toStdNucleo(toupper (Q [qi] [qpos-1]));
                        vcf -> pQ = qpos-1;
                      }else{
                        vcf -> stdN = vcf -> stdN && isStdNucleo(Complement(Q [qi] [qpos-1]));
                        vcf -> cQ = toStdNucleo(toupper (Complement(Q [qi] [qpos-1])));
                        vcf -> pQ = RevC(qpos-1,lenQ);
                      }
                      vcf -> ctxR.push_back(vcf->cR);
                      vcf -> ctxQ.push_back(vcf->cQ);                      
                      vcf -> type = "SNP";
                      vcf -> infType = "SNP";
                      vcf -> inf_sQ = vcf -> pQ;
                      vcf -> inf_fQ = frameQ * frameR;
                      (*li) -> vcfs . push_back (vcf);
                    }
                //-- For the indel
                
                delta = *dit;
                delta = labs(delta);
                vcf = new VCFe_t;
                vcf -> ep = *ei;
                vcf -> lp = *li;
                vcf -> type = "INDEL";

                long inlen =1;
                int l;
                int nsign;
                for(l=1; dit+l < deltas.end(); l++){
                  nsign= *(dit+l) > 0 ? 1 : -1;
                  if(labs(*(dit+l))==1 && sign==nsign)
                    inlen++;
                  else 
                    break;
                }
                dit=dit+inlen-1;
              
                //Whether insert is N-insert
                unsigned int Nlength=0;
                unsigned int Nlentmp=0;
                unsigned int Ncount=0;
                long ni;
                if(sign>0){
                  for(ni=rpos; ni<=rpos+inlen-1; ni++){
                      if(toupper(R [ri] [ni]) == 'N'){
                        Ncount++;
                        Nlentmp++;
                      }else{
                        Nlength=max(Nlength,Nlentmp);
                        Nlentmp=0;
                      }
                  }
                }else{
                  for(ni=qpos; ni<=qpos+inlen-1; ni++){
                      if(toupper(Q [qi] [ni]) == 'N'){
                        Ncount++;
                        Nlentmp++;
                      }else{
                        Nlength=max(Nlength,Nlentmp);
                        Nlentmp=0;
                      }
                  }
                }
                Nlength=max(Nlength,Nlentmp);
                vcf -> Nlength=Nlength;
                vcf -> Nper = (float) Ncount / (float) inlen;
                                
                long k;
                if(sign > 0){
                  if(frameR>0){
                    //Find the first coincident nucleotide
                    rctx=0;
                    do{
                      rctx++;
                    }while( rpos - rctx >= 1 &&
                            qpos - rctx >= 1 &&
                            R [0] [rpos-rctx] != Q [qi][qpos -rctx] && 
                            rctx < delta);
                    
                    if(R [0] [rpos-rctx] != Q [qi][qpos -rctx] &&
                       dit > deltas.begin() && rctx == delta){
                      //IF bump into previous indel before finding 
                      //the first coincident nucleotide.
                      vector<VCFe_t *>::reverse_iterator rit;
                      for(rit= (*li) -> vcfs . rbegin(); 
                          rit < (*li) -> vcfs . rend() &&  (*rit) -> type == "SNPS";
                          rit++){
                          }
                      if(rit < (*li) -> vcfs . rend() && 
                         (*rit) -> type == "INDEL" ){
                        vcf -> ctxR = (*rit) -> ctxR;
                        vcf ->   pR = (*rit) -> pR;
                        vcf ->   cR = (*rit) -> cR;
                        vcf -> ctxQ = (*rit) -> ctxQ;
                        vcf ->   pQ = (*rit) -> pQ;
                        vcf ->   cQ = (*rit) -> cQ;
                        vcf -> stdN = vcf -> stdN && (*rit) -> stdN;
                      }
                    } 
                    
                    if(vcf -> ctxR.size() == 0) {
                      k=rpos-rctx;
                      if ( k < 1 )
                        vcf -> ctxR . push_back (SEQEND_CHAR);
                      else{
                        vcf -> ctxR . push_back (toStdNucleo(toupper (R [0] [k])));
                        vcf -> stdN = vcf -> stdN && isStdNucleo(R [0] [k]);
                      }
                        
                      vcf ->   cR = vcf -> ctxR.c_str()[0];
                    }
                    
                    if(vcf -> ctxQ.size() == 0){    
                      k=qpos-rctx;
                      if ( k < 1 )
                        vcf -> ctxQ . push_back (SEQEND_CHAR);
                      else{
                        vcf -> ctxQ . push_back (toStdNucleo(toupper (Q [qi] [k])));
                        vcf -> stdN = vcf -> stdN && isStdNucleo(toupper( Q [qi] [k]));
                      }
						
                      vcf ->   cQ = vcf -> ctxQ.c_str()[0];
                    }
                    
                    if(!vcf -> pR)
                      vcf -> pR = max((long) 1,rpos-rctx);
                    
                    if(!vcf -> pQ)
                      vcf -> pQ = max((long)1,qpos-rctx);
                    
                    for(k=rpos-rctx+1; k<=rpos+inlen-1; k++){
                      if ( k < 1  || k > lenR)
                        vcf -> ctxR . push_back (SEQEND_CHAR);
                      else{
                        vcf -> ctxR . push_back (toStdNucleo(toupper (R [0] [k]))); 
                        vcf -> stdN = vcf -> stdN && isStdNucleo(toupper (R [0] [k]));
                      }
                    }
                    
                    for(k=qpos-rctx+1; k<qpos; k++){
                      if ( k < 1 || k > lenQ)
                        vcf -> ctxQ . push_back (SEQEND_CHAR);
                      else{
                        vcf -> ctxQ . push_back (toStdNucleo( toupper (Q [qi] [k])));
                        vcf -> stdN = vcf -> stdN && isStdNucleo(toupper(Q [qi] [k]));
                      }
                    }
                    
                    vcf -> infType = "INS";
                    vcf -> inf_sQ   = qpos-1;
                    vcf -> inf_fQ   = frameQ;
                    
                    rpos=rpos+inlen;
                    qpos=qpos;
                  }else{
                    //Find first concident nucleotide
                    long indend=rpos+inlen-1;
                    rctx=0;
                    do{
                      rctx++;
                    }while( indend + rctx <= lenR &&
                            qpos -1 + rctx <= lenQ &&
                            R [ri] [indend+rctx] != Q [qi][qpos -1 + rctx] && 
                            ( dit+1==deltas.end() || rctx < labs(*(dit+1))));

                    if(R [ri] [indend+rctx] != Q [qi][qpos -1 + rctx] &&
                       dit+1 < deltas.end() && rctx == labs(*(dit+1))){
                      //IF bump into next indel before finding 
                      //the first coincident nucleotide.
                      VCFe_t nInd;
                      long tmprpos=indend+rctx;
                      long tmpqpos=qpos-1+rctx;
                      vector<long>::iterator ndit=dit+1;
                      nInd=getIndelSNPvcf(graph,deltas, ndit, 0, tmprpos, tmpqpos, 
                                          ei, li, 
                                          R, Q, 
                                          ri, qi, 
                                          frameR, frameQ, 
                                          lenR, lenQ, 
                                          false);
                      vcf -> ctxR = nInd . ctxR;
                      vcf -> pR   = nInd . pR;
                      vcf -> cR   = nInd . cR;
                      vcf -> ctxQ = nInd . ctxQ;
                      vcf -> pQ   = nInd . pQ;
                      vcf -> cQ   = nInd . cQ;
                      vcf -> stdN = vcf -> stdN && nInd . stdN;
                    } 
                    
                    if(vcf -> ctxR.size() == 0) {
                      k=Norm(min(lenR,indend+rctx), lenR, frameR, graph.datatype);
                      if ( k < 1 || k > lenR )
                        vcf -> ctxR . push_back (SEQEND_CHAR);
                      else{
                        vcf -> ctxR . push_back (toStdNucleo( toupper (R [0] [k])));
                        vcf -> stdN = vcf -> stdN && isStdNucleo(toupper(R [0] [k]));
                      }
                        
                      vcf ->   cR = vcf -> ctxR.c_str()[0];
                    }
                    
                    if(vcf -> ctxQ.size() == 0){    
                      k=min(lenQ,qpos-1+rctx);
                      if ( k<1 || k > lenQ )
                        vcf -> ctxQ . push_back (SEQEND_CHAR);
                      else{
                        vcf -> ctxQ . push_back (toStdNucleo(toupper (Complement ( Q [qi] [k]))));
                        vcf -> stdN = vcf -> stdN && isStdNucleo(Complement( Q [qi] [k]));
                      }
						
                      vcf ->   cQ = vcf -> ctxQ.c_str()[0];
                    }
                    
                    if(!vcf -> pR)
                      vcf -> pR = Norm(min(lenR,indend+rctx), lenR, frameR, graph.datatype);
                    
                    if(!vcf -> pQ)
                      vcf -> pQ = RevC(min(lenQ,qpos-1+rctx), lenQ);
                    
                    for(k=Norm(min(lenR,indend+rctx-1), lenR, frameR, graph.datatype); 
                        k<=Norm(rpos, lenR, frameR, graph.datatype); 
                        k++)
                      if ( k < 1 || k > lenR)
                        vcf -> ctxR . push_back (SEQEND_CHAR);
                      else{
                        vcf -> ctxR . push_back (toStdNucleo (toupper (R [0] [k]))); 
                        vcf -> stdN = vcf -> stdN && isStdNucleo(toupper(R [0] [k]));
                      }
                    
                    
                    for(k=min(lenQ,qpos-1+rctx-1); k>=qpos; k--)
                      if ( k< 1 || k > lenQ )
                        vcf -> ctxQ . push_back (SEQEND_CHAR);
                      else{
                        vcf -> ctxQ . push_back (toStdNucleo( toupper (Complement ( Q [qi] [k]))));
                        vcf -> stdN = vcf -> stdN && isStdNucleo(Complement(Q [qi] [k]));
                      }
                        
                    vcf -> infType = "INS";
                    vcf -> inf_sQ  = RevC(qpos-1,lenQ);
                    vcf -> inf_fQ  = -1*frameQ;
                    
                        
                    rpos=indend+1;
                    qpos=qpos;
                  }
                }else{
                  if(frameR > 0){
                    rctx=0;
                    do{
                      rctx++;
                    }while( rpos - rctx >= 1 &&
                            qpos - rctx >= 1 &&
                            R [0] [rpos-rctx] != Q [qi][qpos -rctx] && 
                            rctx < delta);
                    
                    if(R [0] [rpos-rctx] != Q [qi][qpos -rctx] && 
                      dit>deltas.begin() && rctx == delta){
                      //IF bump into previous indel before finding 
                      //the first coincident nucleotide.
                      vector<VCFe_t *>::reverse_iterator rit;
                      for(rit= (*li) -> vcfs . rbegin(); 
                          rit < (*li) -> vcfs . rend() &&  (*rit) -> type == "SNPS";
                          rit++){
                          }
                      if(rit < (*li) -> vcfs . rend() && 
                         (*rit) -> type == "INDEL"){
                        vcf -> ctxR = (*rit) -> ctxR;
                        vcf ->   pR = (*rit) -> pR;
                        vcf ->   cR = (*rit) -> cR;
                        vcf -> ctxQ = (*rit) -> ctxQ;
                        vcf ->   pQ = (*rit) -> pQ;
                        vcf ->   cQ = (*rit) -> cQ;
                        vcf -> stdN = vcf-> stdN && (*rit) -> stdN;
                      }
                    } 
                    
                    if(vcf -> ctxR.size() == 0) {
                      k=rpos-rctx;
                      if ( k < 1)
                        vcf -> ctxR . push_back (SEQEND_CHAR);
                      else{
                        vcf -> ctxR . push_back (toStdNucleo( toupper (R [0] [k])));
                        vcf -> stdN = vcf -> stdN && isStdNucleo(toupper(R [0] [k]));
                      }
                        
                      vcf ->   cR = vcf -> ctxR.c_str()[0];
                    }
                    
                    if(vcf -> ctxQ.size() == 0){    
                      k=qpos-rctx;
                      if ( k < 1 )
                        vcf -> ctxQ . push_back (SEQEND_CHAR);
                      else{
                        vcf -> ctxQ . push_back (toStdNucleo( toupper (Q [qi] [k])));
                        vcf -> stdN = vcf -> stdN && isStdNucleo(toupper(Q [qi] [k]));
                      }
						
                      vcf ->   cQ = vcf -> ctxQ.c_str()[0];
                    }
                    
                    if(!vcf -> pR)
                      vcf -> pR = max((long)1,rpos-rctx);
                    
                    if(!vcf -> pQ)
                      vcf -> pQ = max((long)1,qpos-rctx);
                    
                    for(k=rpos-rctx+1; k<rpos; k++){
                      if ( k < 1 )
                        vcf -> ctxR . push_back (SEQEND_CHAR);
                      else{
                        vcf -> ctxR . push_back (toStdNucleo( toupper (R [0] [k]))); 
                        vcf -> stdN = vcf -> stdN && isStdNucleo(toupper(R [0] [k]));
                      }
                    }
                    
                    for(k=qpos-rctx+1; k<=qpos+inlen-1; k++){
                      if ( k < 1 )
                        vcf -> ctxQ . push_back (SEQEND_CHAR);
                      else{
                        vcf -> ctxQ . push_back (toStdNucleo( toupper (Q [qi] [k])));
                        vcf -> stdN = vcf -> stdN && isStdNucleo(toupper(Q [qi] [k]));
                      }
                    }
                    
                    vcf -> infType="DEL";
                    vcf -> inf_sQ=qpos-1;
                    vcf -> inf_eQ=qpos-1+inlen;
                    vcf -> inf_fQ=frameQ;                    
                    
                    rpos=rpos;
                    qpos=qpos+inlen;
                  }else{
                    long indend=qpos+inlen-1;
                    rctx=0;
                    do{
                      rctx++;
                    }while( rpos - 1 + rctx <= lenR &&
                            indend + rctx <= lenQ &&
                            R [ri] [rpos-1+rctx] != Q [qi][indend + rctx] && 
                          ( dit+1==deltas.end() || rctx < labs(*(dit+1))));
                    
                    if(R [ri] [rpos-1+rctx] != Q [qi][indend + rctx] && 
                       dit+1 < deltas.end() && rctx == labs(*(dit+1))){
                      //IF bump into next indel before finding 
                      //the first coincident nucleotide.
                      VCFe_t nInd;
                      long tmprpos=rpos-1+rctx;
                      long tmpqpos=indend+rctx;
                      vector<long>::iterator ndit=dit+1;
                      nInd=getIndelSNPvcf(graph,deltas, ndit, 0, tmprpos, tmpqpos, 
                                          ei, li, 
                                          R, Q, 
                                          ri, qi, 
                                          frameR, frameQ, 
                                          lenR, lenQ,
                                          false);
                      vcf -> ctxR = nInd . ctxR;
                      vcf -> pR   = nInd . pR;
                      vcf -> cR   = nInd . cR;
                      vcf -> ctxQ = nInd . ctxQ;
                      vcf -> pQ   = nInd . pQ;
                      vcf -> cQ   = nInd . cQ;
                      vcf -> stdN = vcf -> stdN && nInd.stdN; 
                    }
                    
                    if(!vcf -> pR)
                      vcf -> pR = Norm(min(lenR,rpos-1+rctx), lenR, frameR, graph.datatype);
                    
                    if(!vcf -> pQ)
                      vcf -> pQ = RevC(min(lenQ,indend+rctx), lenQ);
                    
                    for(k=vcf->pR; k<=Norm(rpos, lenR, frameR, graph.datatype); k++)
                      if ( k < 1 )
                        vcf -> ctxR . push_back (SEQEND_CHAR);
                      else{
                        vcf -> ctxR . push_back (toStdNucleo( toupper (R [0] [k]))); 
                        vcf -> stdN = vcf -> stdN && isStdNucleo(toupper(R [0] [k]));
                      }
                    
                    
                    for(k=min(lenQ,indend+rctx); k>=qpos; k--)
                      if ( k > lenQ )
                        vcf -> ctxQ . push_back (SEQEND_CHAR);
                      else{
                        vcf -> ctxQ . push_back (toStdNucleo( toupper ( Complement(Q [qi] [k]))));
                        vcf -> stdN = vcf -> stdN && isStdNucleo(Complement(Q [qi] [k]));
                      }
                        
                    vcf -> infType="DEL";
                    vcf -> inf_sQ=RevC(qpos-1,lenQ);
                    vcf -> inf_eQ=RevC(qpos-1+inlen,lenQ);
                    vcf -> inf_fQ=-1*frameQ;  
                        
                    rpos=rpos;
                    qpos=indend+1;
                  }
                }
                
                if(add2VCFs)
                  (*li) -> vcfs . push_back (vcf);

                return(*vcf);
}

//------------------------------------------------------------ PrintHuman ----//
void PrintHuman (const vector<const SNP_t *> & snps,
                 const DeltaGraph_t & graph)
{
  vector<const SNP_t *>::const_iterator si;
  long dist, distR, distQ;
  int ctxw = 2 * OPT_Context + 1;
  int ctxc = ctxw < 7 ? 7 : ctxw;

  if ( OPT_PrintHeader )
    {
      printf ("%s %s\n%s\n\n",
              graph . refpath . c_str( ), graph . qrypath . c_str( ),
              graph . datatype == NUCMER_DATA ? "NUCMER" : "PROMER");
      printf ("%8s  %5s  %-8s  | ", "[P1]", "[SUB]", "[P2]");
      printf ("%8s %8s  | ", "[BUFF]", "[DIST]");
      if ( OPT_ShowConflict )
        printf ("%4s %4s  | ", "[R]", "[Q]");
      if ( OPT_ShowLength )
        printf ("%8s %8s  | ", "[LEN R]", "[LEN Q]");
      if ( OPT_Context )
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
      if ( OPT_Context )
        for ( int i = 0; i < 2 * ctxc + 7; i ++ )
          putchar ('=');
      printf("================================="
             "==================================\n");
    }

  for ( si = snps . begin( ); si != snps . end( ); ++ si )
    {
      distR = (*si)->pR < (signed long)(*si)->ep->refnode->len - (*si)->pR + 1 ?
        (*si)->pR : (*si)->ep->refnode->len - (*si)->pR + 1;
      distQ = (*si)->pQ < (signed long)(*si)->ep->qrynode->len - (*si)->pQ + 1 ?
        (*si)->pQ : (*si)->ep->qrynode->len - (*si)->pQ + 1;
      dist = distR < distQ ? distR : distQ;

      printf ("%8ld   %c %c   %-8ld  | ",
              (*si)->pR, (*si)->cR, (*si)->cQ, (*si)->pQ);
      printf ("%8ld %8ld  | ", (*si)->buff, dist);
      if ( OPT_ShowConflict )
        printf ("%4d %4d  | ", (*si)->conR, (*si)->conQ);
      if ( OPT_ShowLength )
        printf ("%8ld %8ld  | ",
                (*si)->ep->refnode->len, (*si)->ep->qrynode->len);
      if ( OPT_Context )
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
  vector<const SNP_t *>::const_iterator si;
  long dist, distR, distQ;

  if ( OPT_PrintHeader )
    {
      printf ("%s %s\n%s\n\n",
              graph . refpath . c_str( ), graph . qrypath . c_str( ),
              graph . datatype == NUCMER_DATA ? "NUCMER" : "PROMER");
      printf ("%s\t%s\t%s\t%s\t", "[P1]", "[SUB]", "[SUB]", "[P2]");
      printf ("%s\t%s\t", "[BUFF]", "[DIST]");
      if ( OPT_ShowConflict )
        printf ("%s\t%s\t", "[R]", "[Q]");
      if ( OPT_ShowLength )
        printf ("%s\t%s\t", "[LEN R]", "[LEN Q]");
      if ( OPT_Context )
        printf ("%s\t%s\t", "[CTX R]", "[CTX Q]");
      printf ("%s\t", "[FRM]");
      printf ("%s\n", "[TAGS]");
    }

  for ( si = snps . begin( ); si != snps . end( ); ++ si )
    {
      distR = (*si)->pR < (signed long)(*si)->ep->refnode->len - (*si)->pR + 1 ?
        (*si)->pR : (*si)->ep->refnode->len - (*si)->pR + 1;
      distQ = (*si)->pQ < (signed long)(*si)->ep->qrynode->len - (*si)->pQ + 1 ?
        (*si)->pQ : (*si)->ep->qrynode->len - (*si)->pQ + 1;
      dist = distR < distQ ? distR : distQ;

      printf ("%ld\t%c\t%c\t%ld\t",
              (*si)->pR, (*si)->cR, (*si)->cQ, (*si)->pQ);
      printf ("%ld\t%ld\t", (*si)->buff, dist);
      if ( OPT_ShowConflict )
        printf ("%d\t%d\t", (*si)->conR, (*si)->conQ);
      if ( OPT_ShowLength )
        printf ("%ld\t%ld\t",
                (*si)->ep->refnode->len, (*si)->ep->qrynode->len);
      if ( OPT_Context )
        printf ("%s\t%s\t", (*si)->ctxR.c_str( ), (*si)->ctxQ.c_str( ));
      printf ("%d\t%d\t", (*si)->lp->frmR, (*si)->lp->frmQ);
      printf ("%s\t%s",
              (*si)->ep->refnode->id->c_str( ),
              (*si)->ep->qrynode->id->c_str( ));
      printf ("\n");
    }
}

//------------------------------------------------------------- PrintVCFs ----//
void PrintVCFs (const vector<const VCFe_t *> & vcfs,int argc, char ** argv, DeltaGraph_t & graph){
  vector<const VCFe_t *>::const_iterator vi;

  //Print VCF Headers:
  printf("##fileformat=VCFv4.1\n");

  time_t now;
  now=time(0);
  time(&now);
  struct tm *current=localtime(&now);
  printf("##fileDate=%d/%d/%d\n",1900+current->tm_year,current->tm_mon,current->tm_mday);
  
  printf("##source= ");
  int i;
  for(i=0; i < argc ; i++ )
    printf("%s ",argv[i]);
  printf("\n");
  printf("##reference=%s\n",graph.refpath.c_str());
  std::map<std::string, DeltaNode_t>::iterator nit;
  for(nit=graph.refnodes.begin(); nit != graph.refnodes.end(); nit++)
    printf("##contig=<ID=%s,length=%ld,assembly=%s>\n",nit->second.id->c_str(),nit->second.len,graph.refpath.c_str());
  printf("##INFO=<ID=SNP,Number=0,Type=Flag,Description=\"Whether is SNP\">\n");
  printf("##INFO=<ID=N,Number=1,Type=String,Description=\"Type of N-insert\">\n");
  printf("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n");
  printf("##INFO=<ID=QID,Number=1,Type=String,Description=\"Query ID\">\n");
  printf("##INFO=<ID=QS,Number=.,Type=Integer,Description=\"Query 1-based start position\">\n");
  printf("##INFO=<ID=QE,Number=.,Type=Integer,Description=\"Query 1-based end position\">\n");
  printf("##INFO=<ID=QF,Number=.,Type=Integer,Description=\"Query orientation\">\n");
  printf("##INFO=<ID=AMB,Number=0,Type=Flag,Description=\"Whether ambiguous bases where found and changed\">\n");
  printf("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
  
  
  for ( vi = vcfs . begin( ); vi != vcfs . end( ); ++ vi ){
    if((*vi) -> infType != "SNP" && 
       OPT_filtNs && 
       (*vi) -> Nlength >= OPT_minNlen &&
       (*vi) -> Nper    >= OPT_minNper)
       continue;
    
	  //Parse Ref name
      size_t f;
      string refID((*vi)->ep->refnode->id->c_str());
      do{
        f=refID.find(':',0);
        if(f!=std::string::npos){
          refID.replace(f,1,1,'-');  
        }
      }while(f!=std::string::npos);
      do{
        f=refID.find(' ',0);
        if(f!=std::string::npos){
          refID.replace(f,1,1,'_');  
        }
      }while(f!=std::string::npos);
      //
      
      //Parse qry name
      string qryID((*vi)->ep->qrynode->id->c_str());
      do{
        f=qryID.find(';',0);
        if(f != std::string::npos){
          qryID.replace(f,1,1,'-');
        }
      }while(f!= std::string::npos);
      do{
        f=qryID.find(' ',0);
        if(f != std::string::npos){
          qryID.replace(f,1,1,'_');
        }
      }while(f!= std::string::npos);
      do{
        f=qryID.find(',',0);
        if(f != std::string::npos){
          qryID.replace(f,1,1,'_');
        }
      }while(f!= std::string::npos);
      //
    
      if((*vi)->infType == "SNP"){
				     //CHR  POS  ID
				printf("%s\t%ld\t.\t",
                    refID.c_str( ), 
                    (*vi)->pR);
             // REF ALT QUA FILTER
        printf("%s\t%s\t.\tPASS\t",(*vi)->ctxR.c_str(),(*vi)->ctxQ.c_str());
        printf("%s;QID=%s;QS=%ld",(*vi)->infType.c_str(),qryID.c_str(),(*vi)->inf_sQ);
        if(!(*vi)->stdN)
          printf(";AMB");
        printf("\n");
			}else if((*vi)->infType == "INS"){
				     //CHR  POS  ID
				printf("%s\t%ld\t.\t",
                    refID.c_str( ), 
                    (*vi)->pR);
             // REF ALT QUA FILTER
        printf("%s\t%s\t.\tPASS\t",(*vi)->ctxR.c_str(),(*vi)->ctxQ.c_str());
        printf("SVTYPE=%s;QID=%s;QS=%ld;QF=%d",(*vi)->infType.c_str(),qryID.c_str(),(*vi)->inf_sQ,(*vi)->inf_fQ);
        if(!(*vi)->stdN)
          printf(";AMB");
        printf("\n");
		  
			}else if((*vi)->infType == "DEL"){
				     //CHR  POS  ID
				printf("%s\t%ld\t.\t",
                    refID.c_str( ), 
                    (*vi)->pR);
             // REF ALT QUA FILTER
        printf("%s\t%s\t.\tPASS\t",(*vi)->ctxR.c_str(),(*vi)->ctxQ.c_str());
        printf("SVTYPE=%s;QID=%s;QS=%ld;QE=%ld;QF=%d",(*vi)->infType.c_str(),qryID.c_str(),(*vi)->inf_sQ,(*vi)->inf_eQ,(*vi)->inf_fQ);
        if(!(*vi)->stdN)
          printf(";AMB");
        printf("\n");
		  
			}
      
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
          ((ch = getopt (argc, argv, "ChHIlqrSTVx:NL:%:")) != EOF) )
    switch (ch)
      {
      case 'C':
        OPT_ShowConflict = false;
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
        
      case 'V':
				OPT_SNPs=false;
				break;
        
      case 'N':
        OPT_filtNs=true;
        break;
        
      case 'L':
        OPT_minNlen=atoi(optarg);
        break;
        
      case '%':
        OPT_minNper=atof(optarg);
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
    << "-C            Do not report SNPs from alignments with an ambiguous\n"
    << "              mapping, i.e. only report SNPs where the [R] and [Q]\n"
    << "              columns equal 0 and do not output these columns\n"
    << "-h            Display help information\n"
    << "-H            Do not print the output header\n"
    << "-I            Do not report indels\n"            
    << "-l            Include sequence length information in the output\n"
    << "-q            Sort output lines by query IDs and SNP positions\n"
    << "-r            Sort output lines by reference IDs and SNP positions\n"
    << "-S            Specify which alignments to report by passing\n"
    << "              'show-coords' lines to stdin\n"
    << "-T            Switch to tab-delimited format\n"
    << "-x int        Include x characters of surrounding SNP context in the\n"
    << "              output, default "<< OPT_Context << endl
    << "-V            Output in VCF format\n"
    << "-N -L -%      A sequence range containing at least one contiguous row of Ns\n"
    << "              longer or equal to -L and such that -% of the sequence\n"
    << "              is composed of Ns is considered a N-insert.\n" 
    << "              Setting -N will filter N-inserts from the VCF file.\n"
    << "              (Defaults: -N:false, -L:"<<OPT_minNlen<<", -%:"<<OPT_minNper<<")\n"
    << endl;

  cerr
    << "  Input is the .delta output of either the nucmer or promer program\n"
    << "passed on the command line.\n"
    << "  Output is to stdout, and consists of a list of SNPs (or amino acid\n"
    << "substitutions for promer) with positions and other useful info.\n"
    << "Output will be sorted with -r by default and the [BUFF] column will\n"
    << "always refer to the sequence whose positions have been sorted. This\n"
    << "value specifies the distance from this SNP to the nearest mismatch\n"
    << "(end of alignment, indel, SNP, etc) in the same alignment, while the\n"
    << "[DIST] column specifies the distance from this SNP to the nearest\n"
    << "sequence end. SNPs for which the [R] and [Q] columns are greater than\n"
    << "0 should be evaluated with caution, as these columns specify the\n"
    << "number of other alignments which overlap this position. Use -C to\n"
    << "assure SNPs are only reported from unique alignment regions.\n"
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

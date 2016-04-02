#ifndef  GEAUXDOCK_H
#define  GEAUXDOCK_H


#include <string>

#include "size.h"
#include "toggle.h"



#if defined(__CUDACC__) // NVCC
   #define ALIGN(n) __align__(n)
#elif defined(__GNUC__) // GCC
  #define ALIGN(n) __attribute__((aligned(n)))
#elif defined(__INTEL_COMPILER) // INTEL
  #define ALIGN(n) __attribute__((aligned(n)))
#elif defined(_MSC_VER) // MSVC
  #define ALIGN(n) __declspec(align(n))
#else
  #error "Please provide a definition for MY_ALIGN macro for your host compiler!"
#endif

// 64 bit memory alignment for CPU/Xeon Phi
#define MYALIGN ALIGN(64)

#define RESTRICT __restrict__









struct TraceFile
{
  std::string path;  // ligand trace path
};

struct LigandFile
{
  std::string path; // lig file path
  std::string conf_path;
  std::string id;	// ligand name
  std::string molid; // MOLID;

  int conf_total;	// number of conformation for one lig
  int raw_conf;		// raw total conf in the .sdf without small rmsd excluded

  int lig_natom;	// total effective pts number
  int lnb;	// total bonds number

  //int ensemble_total;	// ENSEMBLE_TOTAL in the .sdf file
};


struct ProteinFile
{
  std::string path;  // prt file path
  std::string id;    // prt name

  int conf_total;	// number of conformations for one lig
  int prt_npoint;	// total effective pts number
  int pnr;	// total bonds number
};


struct LhmFile
{
  std::string path;
  std::string ligand_id;
  int mcs_nrow; // number of mcs positions
};


struct EneParaFile
{
  std::string path;
};

struct NorParaFile
{
  std::string path_a;
  std::string path_b;
};

struct WeightFile
{
  std::string path;
};

struct InputFiles
{
  std::string lig_list;
  std::string prt_list;

  int lig_files_num;
  LigandFile * lig_files;
  LhmFile * lhm_files;

  ProteinFile prt_file;
  EneParaFile enepara_file;
  WeightFile weight_file;
  NorParaFile norpara_file;
  TraceFile trace_file;
};













struct MYALIGN Ligand
{
  // coord_xyz "orig" is under the ligand_ceter system


  // xyz might be rearranged, for PRT computation
  float x1[MAXLIG];		// Ligand x coords, optimized
  float y1[MAXLIG];		// Ligand y coords, optimized
  float z1[MAXLIG];		// Ligand z coords, optimized
  int   t1[MAXLIG];		// atom type, index, might sort ligand point by t
  float c1[MAXLIG];		// atom charge, optimized


  // xyz no rearrange, for MCS computation
  float x2[MAXLIG];
  float y2[MAXLIG];
  float z2[MAXLIG];
  int   t2[MAXLIG];
  float c2[MAXLIG];

  // xyz sort by "t", for KDE computation
  float x3[MAXLIG];
  float y3[MAXLIG];
  float z3[MAXLIG];
  int   t3[MAXLIG];
  float c3[MAXLIG];
  int kde_begin_idx[MAXLIG];
  int kde_end_idx[MAXLIG];




  // coord_center is under the lab system
  float center[3];		// ligand geometric center

  int lig_natom;			// number of ligand atoms
};


struct MYALIGN Protein
{
  // the order reflect the reference order from calc energy function

  int t[MAXPRO];		// effective point type, index, might sort prt points by t
  float x[MAXPRO];		// residue x coord, optimized
  float y[MAXPRO];		// residue y coord, optimized
  float z[MAXPRO];		// residue z coord, optimized

  int cdc[MAXPRO];              // might sort
  float hpp[MAXPRO];            // enepara->hpp[prt->d[i]], optimized
  float ele[MAXPRO];            // dt = prt->t[i] == 0 ? prt->d[i] + 30 : prt->t[i];
                                // enepara->ele[dt], optimized
  int c[MAXPRO];		// effective point class, might sort prt point by c
  int seq3r[MAXPRO];            // index, might sort


  float pocket_center[3];

  int prt_npoint;			// number of protein effective points
};


struct MYALIGN Psp
{
  float psp[MAXLIG][MAXPRO]; // not benifical to implement sparse matrix
};


struct MYALIGN Kde
{
  int t[MAXKDE];		// KDE atom type, sort kde point by t
  float x[MAXKDE];		// KDE x coord, optimized
  float y[MAXKDE];		// KDE y coord, optimized
  float z[MAXKDE];		// KDE z coord, optimized

  int kde_npoint;		// number of kde points
};


// sparse matrix, ELLPACK format
struct MYALIGN Mcs_ELL
{
  int   i[MAX_MCS_ROW][MAX_MCS_COL];         // index
  float x[MAX_MCS_ROW][MAX_MCS_COL];
  float y[MAX_MCS_ROW][MAX_MCS_COL];
  float z[MAX_MCS_ROW][MAX_MCS_COL];
  int ncol[MAX_MCS_ROW];
  float tcc[MAX_MCS_ROW];
};


struct MYALIGN EnePara
{
  // L-J
  float p1a[MAXTP2][MAXTP1];
  float p2a[MAXTP2][MAXTP1];
  float lj0, lj1;

  // electrostatic
  float el1;
  float el0;
  float a1; // 4.0f - 3.0f * el0;
  float b1; // 2.0f * el0 - 3.0f;

  // contact
  float pmf0[MAXTP2][MAXTP1];  // interchange to [lig][prt]
  float pmf1[MAXTP2][MAXTP1];  // interchange to [lig][prt]
  float hdb0[MAXTP2][MAXTP1];  // interchange to [lig][prt]
  float hdb1[MAXTP2][MAXTP1];  // interchange to [lig][prt]

  // hydrophobic
  float hpp[MAXTP4];
  float hpl0[MAXTP2];
  float hpl1[MAXTP2];
  float hpl2[MAXTP2];

  // kde
  float kde2; // -0.5f / (kde * kde)
  float kde3; // powf (kde * sqrtf (2.0f * PI), 3.0f)

  // weights for energy terms
  float a_para[MAXWEI];         // the a parameter for normalization
  float b_para[MAXWEI];         // the b parameter for normalization
  float w[MAXWEI];
};












struct LigCoord
{
  float x[MAXLIG];		// Ligand x coords
  float y[MAXLIG];		// Ligand y coords
  float z[MAXLIG];		// Ligand z coords
  float center[3];		// ligand geometric center
};



struct Ligand0
{
  LigCoord coord_orig;           //                                      used

  int t[MAXLIG];		// atom type                            used
  float c[MAXLIG];		// atom charge                          used
  int n[MAXLIG];		// atom number                          after loading, lig0.n[i] == i + 1

  int lig_natom;			// number of ligand atoms               used
  int lnb;			// number of ligand bonds               NOT USED

  float pocket_center[3];	// pocket center                        used
                                // should belong to "Protein" structure


  float lens_rmsd;		// ensemble rmsd                        NOT USED

  float mw;			// ligand molecular weight              NOT USED
  float logp;			// ligand water/octanol partition coeff NOT USED
  float psa;			// ligand polar surface area            NOT USED
  float mr;			// ligand molar refractivity            NOT USED

  int hbd;			// ligand hydrogen bond donors          NOT USED
  int hba;			// ligand hydrogen bond acceptors       NOT USED

  std::string a[MAXLIG];	// atom name                            NOT USED
  std::string id;		// ligand id                            NOT USED
  std::string smiles;		// ligand smiles                        NOT USED
};



struct Protein0
{
  float x[MAXPRO];		// residue x coord                      used
  float y[MAXPRO];		// residue y coord                      used
  float z[MAXPRO];		// residue z coord                      used
  int n[MAXPRO];		// effective point number               NOT USED
  int t[MAXPRO];		// effective point type                 used
  int c[MAXPRO];		// effective point class                used
  int d[MAXPRO];		// redidue code                         used

  int prt_npoint;			// number of protein effective points   used
  int pnr;			// number of protein residues           NOT USED

  int r[MAXPRO];		// residue number                       replaced
  int seq3[MAXPRO];		// aa sequence numbering                replaced


  int cdc[MAXPRO]; // cdc = (prt_c == 0 && prt_d == 12) || (prt_c == 2)
  int seq3r[MAXPRO]; // seq3r[i] == prt->seq3[prt->r[i]]


  //std::string protein_seq1;  // aa sequence
  //char protein_seq2[MAXPRO]; // aa sequence
};




struct Psp0
{
  float psp[MAXPRO][MAXLIG];                                           // replaced
  int n;			// total number of PSP point           NOT USED
};



struct Kde0
{
  float x[MAXKDE];		// KDE x coord                          used
  float y[MAXKDE];		// KDE y coord                          used
  float z[MAXKDE];		// KDE z coord                          used
  int n[MAXKDE];		// KDE point number                     NOT USED
  int t[MAXKDE];		// KDE atom type                        used

  int kde_npoint;			// number of kde points                 used
  int pns[MAXTP2];		// number of specific kde points        NOT USED
};


struct Mcs0
{
  float x[MAX_MCS_COL];              //                                      used
  float y[MAX_MCS_COL];              //                                      used
  float z[MAX_MCS_COL];              //                                      used


  int   idx_col[MAX_MCS_COL];         // index for sparse matrix              used
  float x2[MAX_MCS_COL];              //                                      used
  float y2[MAX_MCS_COL];              //                                      used
  float z2[MAX_MCS_COL];              //                                      used
  int ncol;			// column number

  float tcc;                    //                                      used
};



struct EnePara0
{
  float vdw[MAXTP1][MAXTP2][2];	// L-J potential                        vdw[prt_t][lig_t][]
  float ele[MAXTP3];		// electrostatic potential              ele[prt_d + 30], ele[prt_t]
  float pmf[MAXTP1][MAXTP2][2];	// contact potential                    pmf[prt_t][lig_t][]
  float hpp[MAXTP4];		// protein hydrophobicity               hpp[prt_d]
  float hpl[MAXTP2][2];		// ligand hydrophobicity                hpl[lig_t][]
  float hdb[MAXTP1][MAXTP2][2];	// ligand hydrophobicity                hdb[prt_t][lig_t][]

  float lj[3];			// L-J params
  float el[2];			// electrostatic params
  float kde;			// kde bandwidth

  float w[MAXWEI];		// weights for energy terms
  float a_para[MAXWEI];         // the a parameter in normalization
  float b_para[MAXWEI];         // the b parameter in normalization
};













struct Temp
{
  float t;
  float minus_beta;
  int order;
};





struct MYALIGN ReplicaMC
{
  int idx_rep; // n_rep, replica

  int idx_lig; // n_lig, ligand
  int idx_prt; // n_prt, protein
  int idx_tmp; // n_tmp, temperature

  float movematrix[6];       // translation x y z, rotation x y z
  float energy[MAXWEI];

  int step; // step counter
  int is_accept; // was the last purturb accepted? record data if is_accpet == 1
};



struct ExchgPara
{
  float floor_temp;   // lowest temperature in all replicas
  float ceiling_temp;   // highest temperature in all replicas
  int num_temp;      // number of temperatures for the same ligand and protein conformations
};


struct McPara
{
  int steps_total;
  int steps_per_dump;
  int steps_per_exchange;

  float move_scale[6]; // translation x y z, rotation x y z

  char outputdir[MAXSTRINGLENG];
  char outputfile[MAXSTRINGLENG];
};

struct McLog
{
  int ac_mc;  // MC acceptance counter
  int acs_mc[MAX_REP];   // MC acceptance counter for all replicas
  int ac_temp_exchg;
  int acs_temp_exchg[MAX_REP]; 
  
  // int ac_lig_exchg;
  // int acs_lig_exchg[MAX_REP];
};




struct Record
{
  ReplicaMC replica[STEPS_PER_DUMP];
  int next_entry;
};




struct MoveVector
{
  float ele[6]; // translation xyz, rotation xyz
};



struct ComplexSize
{
  // replica numbers
  int n_lig; // number of ligand conf
  int n_prt; // number of protein conf
  int n_tmp; // number of temperature
  int n_rep; // n_rep = n_lig * n_prt * n_tmp;

  // residue numbers (of per replica)
  int lig_natom; // number of ligand points
  int prt_npoint; // number of protein points
  int kde_npoint; // number of kde points
  int mcs_nrow; // number of mcs positions
};








struct Complex
{
  // GPU read only arrays
  Ligand lig[MAX_CONF_LIG];
  Protein prt[MAX_CONF_PRT];
  Psp psp;
  Kde kde;
  Mcs_ELL mcs_ell;
  EnePara enepara;
  Temp temp[MAX_TMP];

  // GPU writable arrays that duplicate on multiple GPUs
  ReplicaMC replica[MAX_REP];
  //float *etotal; // used by replica exchange
  //MoveVector *movevector; // used by replcia exchange

  // extra parameters
  ComplexSize size;
  McPara mcpara;


  // file names
  char lig_file[MAX_STR_LENG];
  char prt_file[MAX_STR_LENG];
  char lhm_file[MAX_STR_LENG];
  char enepara_file[MAX_STR_LENG];
  char weight_file[MAX_STR_LENG];


  // GPU read only scalars
  // sizes for multi-device decomposition
  int rep_begin;
  int rep_end;
  //int n_rep;
  //int record_sz;


  // MPI message signal
  int signal;
};


#endif // GEAUXDOCK_H


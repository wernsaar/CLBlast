
// =================================================================================================
// This file is part of the CLBlast project. The project is licensed under Apache Version 2.0. This
// project loosely follows the Google C++ styleguide and uses a tab-size of two spaces and a max-
// width of 100 characters per line.
//
// Author(s):
//   Cedric Nugteren <www.cedricnugteren.nl>
//
// This file contains an optimized matrix-multiplication kernel inspired by the paper by Matsumoto
// et al. and the tutorial on http://www.cedricnugteren.nl/tutorial.php. It is fully configurable
// (and tunable!) using more or less the same parameters/naming conventions as in the paper. It
// supports different data-types (SGEMM/DGEMM/CGEMM/ZGEMM/HGEMM) through a pre-processor define.
//
// Matrices are accessed as follows:
// A: [k*M + m], with 'k' ranging from 0:K and 'm' from 0:M (m,k,m)
// B: [k*N + n], with 'k' ranging from 0:K and 'n' from 0:N (n,k,n)
// C: [n*M + m], with 'n' ranging from 0:N and 'm' from 0:M (m,n,m)
//
// Or as an image (assuming column-major)
//       K                      
//    o-------o                 
//    |       |                 
//  N | [B^T] |                 
//    |       |                 
//    o-------o                 
//        K               N     
//    o-------o        o-----o  
//  M |  [A]  |      M | [C] |  
//    |       |        |     |  
//    o-------o        o-----o  
//                              
//
// This kernel is seperated into two files. This is part 1 out of 2.
//
// =================================================================================================

// Enables loading of this file using the C++ pre-processor's #include (C++11 standard raw string
// literal). Comment-out this line for syntax-highlighting when developing.
R"(

// =================================================================================================

// Parameters set by the tuner or by the database. Here they are given a basic default value in case
// this kernel file is used outside of the CLBlast library.
#ifndef MWG
  #define MWG 8      // Tile-size in dimension M (e.g. 64, 128)
#endif
#ifndef NWG
  #define NWG 8      // Tile-size in dimension N (e.g. 64, 128)
#endif
#ifndef KWG
  #define KWG 8      // Tile-size in dimension K (e.g. 8, 16)
#endif
#ifndef MDIMC
  #define MDIMC 8    // Threads per workgroup in M-dimension (e.g. 8, 16, 32)
#endif
#ifndef NDIMC
  #define NDIMC 8    // Threads per workgroup in N-dimension (e.g. 8, 16, 32)
#endif
#ifndef MDIMA
  #define MDIMA 8    // Re-shaped tile dimension of matrix A: KDIMA * MDIMA
#endif
#ifndef NDIMB
  #define NDIMB 8    // Re-shaped tile dimension of matrix B: KDIMB * NDIMB
#endif
#ifndef KWI
  #define KWI 1      // Unroll factor of the KWG loop (smaller or equal than KWG)
#endif
#ifndef VWM
  #define VWM 1      // Vector width of matrices A and C
#endif
#ifndef VWN
  #define VWN 1      // Vector width of matrix B
#endif
#ifndef STRM
  #define STRM 0     // Use strided access within a thread in the M-dimension (1) or not (0)
#endif
#ifndef STRN
  #define STRN 0     // Use strided access within a thread in the N-dimension (1) or not (0)
#endif
#ifndef SA
  #define SA 0       // Use local/shared memory to cache matrix A (1) or not (0)
#endif
#ifndef SB
  #define SB 0       // Use local/shared memory to cache matrix B (1) or not (0)
#endif


// Settings
#ifndef USE_VECTOR_MAD
  #define USE_VECTOR_MAD 0      // Unroll (0) or don't (1) unroll the vector MAD manually
#endif
#ifndef GLOBAL_MEM_FENCE
  #define GLOBAL_MEM_FENCE 0    // Global synchronisation barrier for potential better performance
#endif

// =================================================================================================

// Data-widths in dimension M
#if VWM == 1
    #define VWM_SHIFT 0
    typedef real realM;
#elif VWM == 2
    #define VWM_SHIFT 1
    typedef real2 realM;
#elif VWM == 4
    #define VWM_SHIFT 2
    typedef real4 realM;
#elif VWM == 8
    #define VWM_SHIFT 3
    typedef real8 realM;
#elif VWM == 16
    #define VWM_SHIFT 4
    typedef real16 realM;
#endif

// Data-widths in dimension N
#if VWN == 1
    #define VWN_SHIFT 0
    typedef real realN;
#elif VWN == 2
    #define VWN_SHIFT 1
    typedef real2 realN;
#elif VWN == 4
    #define VWN_SHIFT 2
    typedef real4 realN;
#elif VWN == 8
    #define VWN_SHIFT 3
    typedef real8 realN;
#elif VWN == 16
    #define VWN_SHIFT 4
    typedef real16 realN;
#endif

#if MDIMA == 1
    #define MDIMA_SHIFT 0
#elif MDIMA == 2
    #define MDIMA_SHIFT 1
#elif MDIMA == 4
    #define MDIMA_SHIFT 2
#elif MDIMA == 8
    #define MDIMA_SHIFT 3
#elif MDIMA == 16
    #define MDIMA_SHIFT 4
#elif MDIMA == 32
    #define MDIMA_SHIFT 5
#elif MDIMA == 64
    #define MDIMA_SHIFT 6
#endif

#if NDIMB == 1
    #define NDIMB_SHIFT 0
#elif NDIMB == 2
    #define NDIMB_SHIFT 1
#elif NDIMB == 4
    #define NDIMB_SHIFT 2
#elif NDIMB == 8
    #define NDIMB_SHIFT 3
#elif NDIMB == 16
    #define NDIMB_SHIFT 4
#elif NDIMB == 32
    #define NDIMB_SHIFT 5
#elif NDIMB == 64
    #define NDIMB_SHIFT 6
#endif

#if MDIMC == 1
    #define MDIMC_SHIFT 0
#elif MDIMC == 2
    #define MDIMC_SHIFT 1
#elif MDIMC == 4
    #define MDIMC_SHIFT 2
#elif MDIMC == 8
    #define MDIMC_SHIFT 3
#elif MDIMC == 16
    #define MDIMC_SHIFT 4
#elif MDIMC == 32
    #define MDIMC_SHIFT 5
#elif MDIMC == 64
    #define MDIMC_SHIFT 6
#endif

#if NDIMC == 1
    #define NDIMC_SHIFT 0
#elif NDIMC == 2
    #define NDIMC_SHIFT 1
#elif NDIMC == 4
    #define NDIMC_SHIFT 2
#elif NDIMC == 8
    #define NDIMC_SHIFT 3
#elif NDIMC == 16
    #define NDIMC_SHIFT 4
#elif NDIMC == 32
    #define NDIMC_SHIFT 5
#elif NDIMC == 64
    #define NDIMC_SHIFT 6
#endif


#if NWG == 1
    #define NWG_SHIFT 0
#elif NWG == 2
    #define NWG_SHIFT 1
#elif NWG == 4
    #define NWG_SHIFT 2
#elif NWG == 8
    #define NWG_SHIFT 3
#elif NWG == 16
    #define NWG_SHIFT 4
#elif NWG == 32
    #define NWG_SHIFT 5
#elif NWG == 64
    #define NWG_SHIFT 6
#elif NWG == 128
    #define NWG_SHIFT 7
#endif

#if MWG == 1
    #define MWG_SHIFT 0
#elif MWG == 2
    #define MWG_SHIFT 1
#elif MWG == 4
    #define MWG_SHIFT 2
#elif MWG == 8
    #define MWG_SHIFT 3
#elif MWG == 16
    #define MWG_SHIFT 4
#elif MWG == 32
    #define MWG_SHIFT 5
#elif MWG == 64
    #define MWG_SHIFT 6
#elif MWG == 128
    #define MWG_SHIFT 7
#endif

#if KWG == 1
    #define KWG_SHIFT 0
#elif KWG == 2
    #define KWG_SHIFT 1
#elif KWG == 4
    #define KWG_SHIFT 2
#elif KWG == 8
    #define KWG_SHIFT 3
#elif KWG == 16
    #define KWG_SHIFT 4
#elif KWG == 32
    #define KWG_SHIFT 5
#elif KWG == 64
    #define KWG_SHIFT 6
#elif KWG == 128
    #define KWG_SHIFT 7
#endif

// Helper parameters based on the above tuning parameters
#define MWI (MWG/MDIMC)               // Work per work-item (M-dimension)
#define NWI (NWG/NDIMC)               // Work per work-item (N-dimension)
#define KDIMA ((MDIMC*NDIMC)/(MDIMA)) // Re-shaped tile dimension of matrix A: KDIMA * MDIMA
#define KDIMB ((MDIMC*NDIMC)/(NDIMB)) // Re-shaped tile dimension of matrix B: KDIMB * NDIMB
#define MWA (MWG/MDIMA)               // Amount of loads-per-thread for matrix A (M-dimension)
#define KWA (KWG/KDIMA)               // Amount of loads-per-thread for matrix A (K-dimension)
#define KWB (KWG/KDIMB)               // Amount of loads-per-thread for matrix B (K-dimension)
#define NWB (NWG/NDIMB)               // Amount of loads-per-thread for matrix B (N-dimension)

#define MWI_SHIFT   ((MWG_SHIFT - MDIMC_SHIFT) > 0 ? (MWG_SHIFT - MDIMC_SHIFT) : 0 )
#define NWI_SHIFT   ((NWG_SHIFT - NDIMC_SHIFT) > 0 ? (NWG_SHIFT - NDIMC_SHIFT) : 0 )
#define KDIMA_SHIFT ((MDIMC_SHIFT + NDIMC_SHIFT - MDIMA_SHIFT) > 0 ? (MDIMC_SHIFT + NDIMC_SHIFT - MDIMA_SHIFT) : 0 )
#define KDIMB_SHIFT ((MDIMC_SHIFT + NDIMC_SHIFT - NDIMB_SHIFT) > 0 ? (MDIMC_SHIFT + NDIMC_SHIFT - NDIMB_SHIFT) : 0 )
#define MWA_SHIFT   ((MWG_SHIFT - MDIMA_SHIFT) > 0 ? (MWG_SHIFT - MDIMA_SHIFT) : 0 )
#define KWA_SHIFT   ((KWG_SHIFT - KDIMA_SHIFT) > 0 ? (KWG_SHIFT - KDIMA_SHIFT) : 0 )
#define KWB_SHIFT   ((KWG_SHIFT - KDIMB_SHIFT) > 0 ? (KWG_SHIFT - KDIMB_SHIFT) : 0 )
#define NWB_SHIFT   ((NWG_SHIFT - NDIMB_SHIFT) > 0 ? (NWG_SHIFT - NDIMB_SHIFT) : 0 )







// =================================================================================================

// Initializes the accumulation registers to zero
inline void InitAccRegisters(realM cpm[NWI][MWI >> VWM_SHIFT]) {
  #pragma unroll
  for (int mi=0; mi<(MWI>>VWM_SHIFT); ++mi) {
    #pragma unroll
    for (int ni=0; ni<NWI; ++ni) {
      #if VWM == 1
        SetToZero(cpm[ni][mi]);
      #elif VWM == 2
        SetToZero(cpm[ni][mi].x);
        SetToZero(cpm[ni][mi].y);
      #elif VWM == 4
        SetToZero(cpm[ni][mi].x);
        SetToZero(cpm[ni][mi].y);
        SetToZero(cpm[ni][mi].z);
        SetToZero(cpm[ni][mi].w);
      #elif VWM == 8
        SetToZero(cpm[ni][mi].s0);
        SetToZero(cpm[ni][mi].s1);
        SetToZero(cpm[ni][mi].s2);
        SetToZero(cpm[ni][mi].s3);
        SetToZero(cpm[ni][mi].s4);
        SetToZero(cpm[ni][mi].s5);
        SetToZero(cpm[ni][mi].s6);
        SetToZero(cpm[ni][mi].s7);
      #elif VWM == 16
        SetToZero(cpm[ni][mi].s0);
        SetToZero(cpm[ni][mi].s1);
        SetToZero(cpm[ni][mi].s2);
        SetToZero(cpm[ni][mi].s3);
        SetToZero(cpm[ni][mi].s4);
        SetToZero(cpm[ni][mi].s5);
        SetToZero(cpm[ni][mi].s6);
        SetToZero(cpm[ni][mi].s7);
        SetToZero(cpm[ni][mi].s8);
        SetToZero(cpm[ni][mi].s9);
        SetToZero(cpm[ni][mi].sA);
        SetToZero(cpm[ni][mi].sB);
        SetToZero(cpm[ni][mi].sC);
        SetToZero(cpm[ni][mi].sD);
        SetToZero(cpm[ni][mi].sE);
        SetToZero(cpm[ni][mi].sF);
      #endif
    }
  }
}

// =================================================================================================

// Caches global off-chip memory into local (shared) memory on-chip. This function is specific for
// caching the A input matrix.
#if SA == 1
inline void GlobalToLocalA(const __global realM* restrict agm, __local realM* alm,
                           const int kSizeM, const int tid, const int kwg) {

  #if STRM == 0
    const int la0_M1 = (tid - (tid & -MDIMA)) << ((MWA_SHIFT - VWM_SHIFT) > 0 ? (MWA_SHIFT - VWM_SHIFT) : 0);
  #else
    const int la0 = tid - (tid & -MDIMA);
  #endif

  const int GroupID0_M1 = GetGroupID0() << ((MWG_SHIFT - VWM_SHIFT) > 0 ? (MWG_SHIFT - VWM_SHIFT) : 0);
  const int  la1_M1 = (tid >> MDIMA_SHIFT) << KWA_SHIFT;
  const int  kSizeMxVWM = kSizeM >> VWM_SHIFT;
 

  #pragma unroll
  for (int mia=0; mia<(MWA >> VWM_SHIFT); ++mia) {

    #if STRM == 0
      const int mg = mia + la0_M1;
    #elif STRM == 1
      const int mg = la0 + (mia << MDIMA_SHIFT);
    #endif

    const int idm = mg + GroupID0_M1;

    #pragma unroll
    for (int kia=0; kia<KWA; ++kia) {
        // Computes the indices for the global memory
        int kg = kia + la1_M1;

        // Loads the data from global memory (not transposed) into the local memory
        #if USE_MAD24 == 1
          alm[(kg << ((MWG_SHIFT - VWM_SHIFT) >0 ? (MWG_SHIFT - VWM_SHIFT) : 0 )) + mg] = agm[mad24(kg+kwg,kSizeMxVWM , idm)];
        #else
          alm[(kg << ((MWG_SHIFT - VWM_SHIFT) >0 ? (MWG_SHIFT - VWM_SHIFT) : 0 )) + mg] = agm[(kg+kwg)*kSizeMxVWM + idm];
        #endif 

    }
  }
}
#endif

// Same as above, but now for the B input matrix
#if SB == 1
inline void GlobalToLocalB(const __global realN* restrict bgm, __local realN* blm,
                           const int kSizeN, const int tid, const int kwg) {

  const int lb1_M1 = (tid >> NDIMB_SHIFT) << KWB_SHIFT;

  #if STRN == 0
    const int lb0_M1 = (tid - (tid & -NDIMB)) << ((NWB_SHIFT - VWN_SHIFT) >0 ? (NWB_SHIFT - VWN_SHIFT) : 0);
  #else
    const int lb0 = tid - (tid & -NDIMB);
  #endif

  const int GroupID1_M1 = GetGroupID1() << ((NWG_SHIFT - VWN_SHIFT) >0 ? (NWG_SHIFT - VWN_SHIFT) : 0);
  const int kSizeNxVWN = kSizeN >> VWN_SHIFT;

  #pragma unroll
  for (int kib=0; kib<KWB; ++kib) {

    const int kg = kib + lb1_M1;
    const int kg_M1 = kg << ((NWG_SHIFT - VWN_SHIFT) >0 ? (NWG_SHIFT - VWN_SHIFT) : 0);
    #if USE_MAD24 == 1
      const int idk_M1 = mad24(kg + kwg , kSizeNxVWN , GroupID1_M1);
    #else
      const int idk_M1 = (kg + kwg) * kSizeNxVWN + GroupID1_M1;
    #endif

    #pragma unroll
    for (int nib=0; nib<(NWB >> VWN_SHIFT); ++nib) {

        // Computes the indices based on strided/non-strided access
        #if STRN == 0
          int ng = nib + lb0_M1;
        #elif STRN == 1
          int ng = lb0 + (nib << NDIMB_SHIFT);
        #endif

        // Loads the data from global memory (transposed) into the local memory
        blm[kg_M1 + ng] = bgm[idk_M1 + ng];

    }
  }
}
#endif

// =================================================================================================

// Caches global off-chip memory directly into per-thread private memory (registers). This function
// is specific for caching the A input matrix.
#if SA == 0
inline void GlobalToPrivateA(const __global realM* restrict agm, realM apm[MWI>>VWM_SHIFT],
                             const int kSizeM, const int idk, const int kwg) {


  #if STRM == 0
    int LocalID0_M1 = get_local_id(0) << (MWI_SHIFT - VWM_SHIFT);
  #else
    int local_id0 = get_local_id(0);
  #endif 

  #if USE_MAD24 == 1
    const int idk_M1 = mad24(idk,(kSizeM >> VWM_SHIFT) , (int) GetGroupID0() << ((MWG_SHIFT - VWM_SHIFT) >0 ? (MWG_SHIFT - VWM_SHIFT) : 0));
  #else
    const int idk_M1 = idk*(kSizeM >> VWM_SHIFT) + ((int) GetGroupID0() << ((MWG_SHIFT - VWM_SHIFT) >0 ? (MWG_SHIFT - VWM_SHIFT) : 0));
  #endif

  #pragma unroll
  for (int mi=0; mi<(MWI >> VWM_SHIFT); ++mi) {

      // Computes the indices based on strided/non-strided access
      #if STRM == 0
        int mg = mi + LocalID0_M1;
      #elif STRM == 1
        int mg = local_id0 + (mi << MDIMC_SHIFT);
      #endif

      // Loads the data from global memory (not transposed) and stores into registers
      apm[mi] = agm[idk_M1 + mg];


  }
}
#endif

// Same as above, but now for the B input matrix
#if SB == 0
inline void GlobalToPrivateB(const __global realN* restrict bgm, realN bpm[NWI >> VWN_SHIFT],
                             const int kSizeN, const int idk) {

  #if STRN == 0
    int LocalID1_M1 = get_local_id(1) << ((NWI_SHIFT - VWN_SHIFT) >0 ? (NWI_SHIFT - VWN_SHIFT) : 0);
  #else
    int local_id1 = get_local_id(1); 
  #endif

  #if USE_MAD24 == 1
    const int idk_M1 = mad24(idk,(kSizeN >> VWN_SHIFT) , ((int) GetGroupID1() << ((NWG_SHIFT - VWN_SHIFT) >0 ? (NWG_SHIFT - VWN_SHIFT) : 0)));
  #else
    const int idk_M1 = idk*(kSizeN >> VWN_SHIFT) + ((int) GetGroupID1() << ((NWG_SHIFT - VWN_SHIFT) >0 ? (NWG_SHIFT - VWN_SHIFT) : 0));
  #endif

  #pragma unroll
  for (int ni=0; ni<(NWI >> VWN_SHIFT); ++ni) {

      // Computes the indices based on strided/non-strided access
      #if STRN == 0
        int ng = ni + LocalID1_M1;
      #elif STRN == 1
        int ng = local_id1 + (ni << NDIMC_SHIFT);
      #endif

      // Loads the data from global memory (transposed) and stores into registers
      bpm[ni] = bgm[idk_M1 + ng];

  }
}
#endif

// =================================================================================================

// Caches on-chip local memory into per-thread private memory (registers). This function is specific
// for caching the A input matrix.
#if SA == 1
inline void LocalToPrivateA(__local realM* alm, realM apm[MWI >> VWM_SHIFT], const int kg) {

  #if STRM == 0
    int LocalID0_M1 = get_local_id(0) << ((MWI_SHIFT - VWM_SHIFT) >0 ? (MWI_SHIFT - VWM_SHIFT) : 0);
  #else
    int local_id0 = get_local_id(0);
  #endif

  int kg_M1 = kg << ((MWG_SHIFT - VWM_SHIFT) >0 ? (MWG_SHIFT - VWM_SHIFT) : 0);

  #pragma unroll
  for (int mi=0; mi<(MWI >> VWM_SHIFT); ++mi) {

      #if STRM == 0
        int mg = mi + LocalID0_M1;
      #elif STRM == 1
        int mg = local_id0 + (mi << MDIMC_SHIFT);
      #endif

      apm[mi] = alm[kg_M1 + mg];

  }
}
#endif

// Same as above, but now for the B input matrix
#if SB == 1
inline void LocalToPrivateB(__local realN* blm, realN bpm[NWI >> VWN_SHIFT], const int kg) {

  #if STRN == 0
    int local_id1_M1 = get_local_id(1) << ((NWI_SHIFT - VWN_SHIFT) >0 ? (NWI_SHIFT - VWN_SHIFT) : 0);
  #else
    int local_id1 = get_local_id(1);
  #endif

  const int kg_M1 = kg << ((NWG_SHIFT - VWN_SHIFT) >0 ? (NWG_SHIFT - VWN_SHIFT) : 0);

  #pragma unroll
  for (int ni=0; ni<(NWI >> VWN_SHIFT); ++ni) {

    #if STRN == 0
      int ng = ni + local_id1_M1;
    #elif STRN == 1
      int ng = local_id1 + (ni << NDIMC_SHIFT);
    #endif

    bpm[ni] = blm[kg_M1 + ng];
  }
}
#endif

// =================================================================================================

// End of the C++11 raw string literal
)"

// =================================================================================================

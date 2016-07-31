
// =================================================================================================
// This file is part of the CLBlast project. The project is licensed under Apache Version 2.0. This
// project loosely follows the Google C++ styleguide and uses a tab-size of two spaces and a max-
// width of 100 characters per line.
//
// Author(s):
//   Cedric Nugteren <www.cedricnugteren.nl>
//
// This is part 2 of 2 of the GEMM kernel. See part 1 for more information.
//
// =================================================================================================

// Enables loading of this file using the C++ pre-processor's #include (C++11 standard raw string
// literal). Comment-out this line for syntax-highlighting when developing.
R"(

// =================================================================================================

// The vectorised multiply-add function
inline realM MultiplyAddVector(realM cvec, const realM avec, const real bval) {
  #if (USE_VECTOR_MAD == 1) && ((PRECISION == 32) || (PRECISION == 64)) 
    #if USE_CL_MAD == 1
      cvec = mad(avec,bval,cvec);
    #elif USE_CL_FMA == 1
      cvec = fma(avec,bval,cvec);
    #else
      cvec += avec * bval;
    #endif 

  #else
    #if VWM == 1
      MultiplyAdd(cvec,    avec,    bval);
    #elif VWM == 2
      MultiplyAdd(cvec.x , avec.x,  bval);
      MultiplyAdd(cvec.y , avec.y,  bval);
    #elif VWM == 4
      MultiplyAdd(cvec.x , avec.x,  bval);
      MultiplyAdd(cvec.y , avec.y,  bval);
      MultiplyAdd(cvec.z , avec.z,  bval);
      MultiplyAdd(cvec.w , avec.w,  bval);
    #elif VWM == 8
      MultiplyAdd(cvec.s0, avec.s0, bval);
      MultiplyAdd(cvec.s1, avec.s1, bval);
      MultiplyAdd(cvec.s2, avec.s2, bval);
      MultiplyAdd(cvec.s3, avec.s3, bval);
      MultiplyAdd(cvec.s4, avec.s4, bval);
      MultiplyAdd(cvec.s5, avec.s5, bval);
      MultiplyAdd(cvec.s6, avec.s6, bval);
      MultiplyAdd(cvec.s7, avec.s7, bval);
    #elif VWM == 16
      MultiplyAdd(cvec.s0, avec.s0, bval);
      MultiplyAdd(cvec.s1, avec.s1, bval);
      MultiplyAdd(cvec.s2, avec.s2, bval);
      MultiplyAdd(cvec.s3, avec.s3, bval);
      MultiplyAdd(cvec.s4, avec.s4, bval);
      MultiplyAdd(cvec.s5, avec.s5, bval);
      MultiplyAdd(cvec.s6, avec.s6, bval);
      MultiplyAdd(cvec.s7, avec.s7, bval);
      MultiplyAdd(cvec.s8, avec.s8, bval);
      MultiplyAdd(cvec.s9, avec.s9, bval);
      MultiplyAdd(cvec.sA, avec.sA, bval);
      MultiplyAdd(cvec.sB, avec.sB, bval);
      MultiplyAdd(cvec.sC, avec.sC, bval);
      MultiplyAdd(cvec.sD, avec.sD, bval);
      MultiplyAdd(cvec.sE, avec.sE, bval);
      MultiplyAdd(cvec.sF, avec.sF, bval);
    #endif
  #endif
  return cvec;
}

// Performs the actual computation: Cpm += Apm * Bpm
inline void MultiplyAccumulate(realM cpm[NWI][MWI >> VWM_SHIFT], realM apm[MWI >> VWM_SHIFT], realN bpm[NWI >> VWN_SHIFT]) {
  #pragma unroll
  for (uint ni=0; ni<(NWI >> VWN_SHIFT); ++ni) {

    const uint NIxVWN = ni << VWN_SHIFT;

    #pragma unroll
    for (uint mi=0; mi<(MWI >> VWM_SHIFT); ++mi) {
      const realM aval = apm[mi];
      #if VWN == 1
        cpm[NIxVWN + 0][mi] = MultiplyAddVector(cpm[NIxVWN + 0][mi], aval, bpm[ni]);
      #elif VWN == 2
        cpm[NIxVWN + 0][mi] = MultiplyAddVector(cpm[NIxVWN + 0][mi], aval, bpm[ni].x);
        cpm[NIxVWN + 1][mi] = MultiplyAddVector(cpm[NIxVWN + 1][mi], aval, bpm[ni].y);
      #elif VWN == 4
        cpm[NIxVWN + 0][mi] = MultiplyAddVector(cpm[NIxVWN + 0][mi], aval, bpm[ni].x);
        cpm[NIxVWN + 1][mi] = MultiplyAddVector(cpm[NIxVWN + 1][mi], aval, bpm[ni].y);
        cpm[NIxVWN + 2][mi] = MultiplyAddVector(cpm[NIxVWN + 2][mi], aval, bpm[ni].z);
        cpm[NIxVWN + 3][mi] = MultiplyAddVector(cpm[NIxVWN + 3][mi], aval, bpm[ni].w);
      #elif VWN == 8
        cpm[NIxVWN + 0][mi] = MultiplyAddVector(cpm[NIxVWN + 0][mi], aval, bpm[ni].s0);
        cpm[NIxVWN + 1][mi] = MultiplyAddVector(cpm[NIxVWN + 1][mi], aval, bpm[ni].s1);
        cpm[NIxVWN + 2][mi] = MultiplyAddVector(cpm[NIxVWN + 2][mi], aval, bpm[ni].s2);
        cpm[NIxVWN + 3][mi] = MultiplyAddVector(cpm[NIxVWN + 3][mi], aval, bpm[ni].s3);
        cpm[NIxVWN + 4][mi] = MultiplyAddVector(cpm[NIxVWN + 4][mi], aval, bpm[ni].s4);
        cpm[NIxVWN + 5][mi] = MultiplyAddVector(cpm[NIxVWN + 5][mi], aval, bpm[ni].s5);
        cpm[NIxVWN + 6][mi] = MultiplyAddVector(cpm[NIxVWN + 6][mi], aval, bpm[ni].s6);
        cpm[NIxVWN + 7][mi] = MultiplyAddVector(cpm[NIxVWN + 7][mi], aval, bpm[ni].s7);
      #elif VWN == 16
        cpm[NIxVWN + 0 ][mi] = MultiplyAddVector(cpm[NIxVWN + 0 ][mi], aval, bpm[ni].s0);
        cpm[NIxVWN + 1 ][mi] = MultiplyAddVector(cpm[NIxVWN + 1 ][mi], aval, bpm[ni].s1);
        cpm[NIxVWN + 2 ][mi] = MultiplyAddVector(cpm[NIxVWN + 2 ][mi], aval, bpm[ni].s2);
        cpm[NIxVWN + 3 ][mi] = MultiplyAddVector(cpm[NIxVWN + 3 ][mi], aval, bpm[ni].s3);
        cpm[NIxVWN + 4 ][mi] = MultiplyAddVector(cpm[NIxVWN + 4 ][mi], aval, bpm[ni].s4);
        cpm[NIxVWN + 5 ][mi] = MultiplyAddVector(cpm[NIxVWN + 5 ][mi], aval, bpm[ni].s5);
        cpm[NIxVWN + 6 ][mi] = MultiplyAddVector(cpm[NIxVWN + 6 ][mi], aval, bpm[ni].s6);
        cpm[NIxVWN + 7 ][mi] = MultiplyAddVector(cpm[NIxVWN + 7 ][mi], aval, bpm[ni].s7);
        cpm[NIxVWN + 8 ][mi] = MultiplyAddVector(cpm[NIxVWN + 8 ][mi], aval, bpm[ni].s8);
        cpm[NIxVWN + 9 ][mi] = MultiplyAddVector(cpm[NIxVWN + 9 ][mi], aval, bpm[ni].s9);
        cpm[NIxVWN + 10][mi] = MultiplyAddVector(cpm[NIxVWN + 10][mi], aval, bpm[ni].sA);
        cpm[NIxVWN + 11][mi] = MultiplyAddVector(cpm[NIxVWN + 11][mi], aval, bpm[ni].sB);
        cpm[NIxVWN + 12][mi] = MultiplyAddVector(cpm[NIxVWN + 12][mi], aval, bpm[ni].sC);
        cpm[NIxVWN + 13][mi] = MultiplyAddVector(cpm[NIxVWN + 13][mi], aval, bpm[ni].sD);
        cpm[NIxVWN + 14][mi] = MultiplyAddVector(cpm[NIxVWN + 14][mi], aval, bpm[ni].sE);
        cpm[NIxVWN + 15][mi] = MultiplyAddVector(cpm[NIxVWN + 15][mi], aval, bpm[ni].sF);
      #endif
    }
  }
}

// Merges the results in Cpm with the global array in Cgm. This also performs the multiplication
// with the constants: Cgm = alpha*A*B + beta*Cgm = alpha*Cpm + beta*Cgm
inline void StoreResults(__global realM* cgm, realM cpm[NWI][MWI >> VWM_SHIFT], const int kSizeM,
                         const real alpha, const real beta) {


  const uint kSizeMxVWM  = kSizeM >> VWM_SHIFT;
  const uint LocalID1_M1 = get_local_id(1) << VWN_SHIFT;

  #if STRN == 0
    const uint LocalID1_M2 = (int) get_local_id(1) << NWI_SHIFT;
  #endif

  #if STRM == 0
    const uint LocalID0_M1 = (int) get_local_id(0) << ((MWI_SHIFT - VWM_SHIFT) >0 ? (MWI_SHIFT - VWM_SHIFT) : 0);
  #else
    const uint local_id0 = get_local_id(0);
  #endif

  #if USE_MAD24 == 1
    const uint GroupID1_M2 = mad24((uint) GetGroupID1() << NWG_SHIFT, kSizeMxVWM, (uint) GetGroupID0() << ((MWG_SHIFT - VWM_SHIFT) >0 ? (MWG_SHIFT - VWM_SHIFT) :0 ));
  #else
    const uint GroupID1_M2 = (GetGroupID1() << NWG_SHIFT) * kSizeMxVWM + (GetGroupID0() << ((MWG_SHIFT - VWM_SHIFT) >0 ? (MWG_SHIFT - VWM_SHIFT) : 0));
  #endif

  #pragma unroll
  for (uint ni=0; ni<NWI; ++ni) {

    #if STRN == 0
      const uint ng = ni + LocalID1_M2;
    #elif STRN == 1
      const uint ng = (ni - (ni & -VWN)) + LocalID1_M1 + ((ni >> VWN_SHIFT) << (VWN_SHIFT + NDIMC_SHIFT));
    #endif

    #if USE_MAD24 == 1
      const uint idn_M1 = mad24(ng, kSizeMxVWM , GroupID1_M2);
    #else
      const uint idn_M1 = ng*kSizeMxVWM + GroupID1_M2;
    #endif

    #pragma unroll
    for (uint mi=0; mi<(MWI >> VWM_SHIFT); ++mi) {

      #if STRM == 0
        const uint mg = mi + LocalID0_M1;
      #elif STRM == 1
        const uint mg = local_id0 + (mi << MDIMC_SHIFT);
      #endif


      // The final multiplication with alpha and the addition with beta*C
      const uint index = idn_M1 + mg;


      #if (USE_VECTOR_MAD == 1) && ((PRECISION == 32) || (PRECISION == 64)) 
        #if USE_CL_MAD == 1
          cgm[index]=mad(cpm[ni][mi],(realM)alpha , cgm[index]*(realM) beta);
        #elif USE_CL_FMA == 1
          cgm[index]=fma(cpm[ni][mi],(realM)alpha , cgm[index]*(realM) beta);
        #else
          cgm[index]=cpm[ni][mi] *(realM)alpha + cgm[index] * (realM) beta;
        #endif
      #else


      realM result;
      realM xval = cpm[ni][mi];
      realM yval = cgm[index];

      #if VWM == 1
        AXPBY(result, alpha, xval, beta, yval);
      #elif VWM == 2
        AXPBY(result.x, alpha, xval.x, beta, yval.x);
        AXPBY(result.y, alpha, xval.y, beta, yval.y);
      #elif VWM == 4
        AXPBY(result.x, alpha, xval.x, beta, yval.x);
        AXPBY(result.y, alpha, xval.y, beta, yval.y);
        AXPBY(result.z, alpha, xval.z, beta, yval.z);
        AXPBY(result.w, alpha, xval.w, beta, yval.w);
      #elif VWM == 8
        AXPBY(result.s0, alpha, xval.s0, beta, yval.s0);
        AXPBY(result.s1, alpha, xval.s1, beta, yval.s1);
        AXPBY(result.s2, alpha, xval.s2, beta, yval.s2);
        AXPBY(result.s3, alpha, xval.s3, beta, yval.s3);
        AXPBY(result.s4, alpha, xval.s4, beta, yval.s4);
        AXPBY(result.s5, alpha, xval.s5, beta, yval.s5);
        AXPBY(result.s6, alpha, xval.s6, beta, yval.s6);
        AXPBY(result.s7, alpha, xval.s7, beta, yval.s7);
      #elif VWM == 16
        AXPBY(result.s0, alpha, xval.s0, beta, yval.s0);
        AXPBY(result.s1, alpha, xval.s1, beta, yval.s1);
        AXPBY(result.s2, alpha, xval.s2, beta, yval.s2);
        AXPBY(result.s3, alpha, xval.s3, beta, yval.s3);
        AXPBY(result.s4, alpha, xval.s4, beta, yval.s4);
        AXPBY(result.s5, alpha, xval.s5, beta, yval.s5);
        AXPBY(result.s6, alpha, xval.s6, beta, yval.s6);
        AXPBY(result.s7, alpha, xval.s7, beta, yval.s7);
        AXPBY(result.s8, alpha, xval.s8, beta, yval.s8);
        AXPBY(result.s9, alpha, xval.s9, beta, yval.s9);
        AXPBY(result.sA, alpha, xval.sA, beta, yval.sA);
        AXPBY(result.sB, alpha, xval.sB, beta, yval.sB);
        AXPBY(result.sC, alpha, xval.sC, beta, yval.sC);
        AXPBY(result.sD, alpha, xval.sD, beta, yval.sD);
        AXPBY(result.sE, alpha, xval.sE, beta, yval.sE);
        AXPBY(result.sF, alpha, xval.sF, beta, yval.sF);
      #endif
      cgm[index] = result;
      #endif
    }
  }
}

// =================================================================================================

// Main body of the matrix-multiplication algorithm. It calls the (inlined) functions above.
inline void XgemmBody(const int kSizeM, const int kSizeN, const int kSizeK,
                      const __global realM* restrict agm, const __global realN* restrict bgm,
                      __global realM* cgm, realM cpm[NWI][MWI >> VWM_SHIFT]
                      #if SA == 1 && SB == 1
                        , __local realM* alm, __local realN* blm
                      #elif SA == 1
                        , __local realM* alm
                      #elif SB == 1
                        , __local realN* blm
                      #endif
                      ) {

  // Allocates workitem-private memory (registers)
  realM apm[MWI >> VWM_SHIFT];
  realN bpm[NWI >> VWN_SHIFT];

  // Combined thread identifier (volatile to disable caching)
  #if SA == 1 || SB == 1
    volatile uint tid = get_local_id(0) + (get_local_id(1) << MDIMC_SHIFT);
  #endif

  // Initializes the accumulation registers
  InitAccRegisters(cpm);

  // Loops over all workgroup tiles
  // #pragma unroll
  for (uint kwg=0; kwg<kSizeK; kwg+=KWG) {

    // Loads data: off-chip --> local (matrix A)
    #if SA == 1
      GlobalToLocalA(agm, alm, kSizeM, tid, kwg);
    #endif
    // Loads data: off-chip --> local (matrix B)
    #if SB == 1
      GlobalToLocalB(bgm, blm, kSizeN, tid, kwg);
    #endif
    #if SA == 1 || SB == 1
      barrier(CLK_LOCAL_MEM_FENCE);
    #endif

    // Loops over all workitem tiles, unrolled by a factor KWI
    // #pragma unroll
    for (uint pwi=0; pwi<KWG; pwi+=KWI) {

      #pragma unroll
      for (uint pit=0; pit<KWI; ++pit) {

        #if SA == 0 || SB == 0
          const uint idk = kwg + pwi + pit;
        #endif
        #if SA == 1 || SB == 1
          const uint kg = pwi+pit;
        #endif

        // Loads data: local --> private (matrix A)
        #if SA == 1
          LocalToPrivateA(alm, apm, kg);
        // Loads data: off-chip --> private (matrix A)
        #else
          GlobalToPrivateA(agm, apm, kSizeM, idk, kwg);
        #endif

        // Loads data: local --> private (matrix B)
        #if SB == 1
          LocalToPrivateB(blm, bpm, kg);
        // Loads data: off-chip --> private (matrix B)
        #else
          GlobalToPrivateB(bgm, bpm, kSizeN, idk);
        #endif

        // Performs the accumulation (Cpm += Apm * Bpm)
        MultiplyAccumulate(cpm, apm, bpm);
      }
    }
    #if SA == 1 || SB == 1
      barrier(CLK_LOCAL_MEM_FENCE);
    #endif
  }
  #if GLOBAL_MEM_FENCE == 1
    barrier(CLK_GLOBAL_MEM_FENCE);
  #endif
}

// =================================================================================================
// The upper-triangular and lower-triangular kernels are only used in special cases
#if defined(ROUTINE_SYRK) || defined(ROUTINE_HERK) || defined(ROUTINE_SYR2K) || defined(ROUTINE_HER2K)

// Main entry point of the kernel. This is the upper-triangular version.
__attribute__((reqd_work_group_size(MDIMC, NDIMC, 1)))
__kernel void XgemmUpper(const int kSizeN, const int kSizeK,
                         const __constant real* restrict arg_alpha,
                         const __constant real* restrict arg_beta,
                         const __global realM* restrict agm,
                         const __global realN* restrict bgm,
                         __global realM* cgm) {
  const real alpha = arg_alpha[0];
  const real beta = arg_beta[0];

  // Skip these threads if they do not contain threads contributing to the upper-triangle
  if (GetGroupID1()*NWG < GetGroupID0()*MWG) {
    return;
  }

  // Allocates workgroup-private memory (local memory)
  #if SA == 1
    __local realM alm[KWG * MWG/VWM];
  #endif
  #if SB == 1
    __local realN blm[KWG * NWG/VWN];
  #endif

  // Computes the matrix-multiplication and stores the result in register memory
  realM cpm[NWI][MWI/VWM];
  #if SA == 1 && SB == 1
    XgemmBody(kSizeN, kSizeN, kSizeK, agm, bgm, cgm, cpm, alm, blm);
  #elif SA == 1
    XgemmBody(kSizeN, kSizeN, kSizeK, agm, bgm, cgm, cpm, alm);
  #elif SB == 1
    XgemmBody(kSizeN, kSizeN, kSizeK, agm, bgm, cgm, cpm, blm);
  #else
    XgemmBody(kSizeN, kSizeN, kSizeK, agm, bgm, cgm, cpm);
  #endif

  // Stores an MWG * NWG tile of results and performs the multiplication with alpha and beta
  StoreResults(cgm, cpm, kSizeN, alpha, beta);
}

// Main entry point of the kernel. This is the lower-triangular version.
__attribute__((reqd_work_group_size(MDIMC, NDIMC, 1)))
__kernel void XgemmLower(const int kSizeN, const int kSizeK,
                         const __constant real* restrict arg_alpha,
                         const __constant real* restrict arg_beta,
                         const __global realM* restrict agm,
                         const __global realN* restrict bgm,
                         __global realM* cgm) {
  const real alpha = arg_alpha[0];
  const real beta = arg_beta[0];

  // Skip these threads if they do not contain threads contributing to the lower-triangle
  if (GetGroupID1()*NWG > GetGroupID0()*MWG) {
    return;
  }

  // Allocates workgroup-private memory (local memory)
  #if SA == 1
    __local realM alm[KWG * MWG/VWM];
  #endif
  #if SB == 1
    __local realN blm[KWG * NWG/VWN];
  #endif

  // Computes the matrix-multiplication and stores the result in register memory
  realM cpm[NWI][MWI/VWM];
  #if SA == 1 && SB == 1
    XgemmBody(kSizeN, kSizeN, kSizeK, agm, bgm, cgm, cpm, alm, blm);
  #elif SA == 1
    XgemmBody(kSizeN, kSizeN, kSizeK, agm, bgm, cgm, cpm, alm);
  #elif SB == 1
    XgemmBody(kSizeN, kSizeN, kSizeK, agm, bgm, cgm, cpm, blm);
  #else
    XgemmBody(kSizeN, kSizeN, kSizeK, agm, bgm, cgm, cpm);
  #endif

  // Stores an MWG * NWG tile of results and performs the multiplication with alpha and beta
  StoreResults(cgm, cpm, kSizeN, alpha, beta);
}

// =================================================================================================
// If not using a triangular version, include the regular kernel
#else

// Main entry point of the kernel. This is the regular full version.
__attribute__((reqd_work_group_size(MDIMC, NDIMC, 1)))
__kernel void Xgemm(const int kSizeM, const int kSizeN, const int kSizeK,
                    const __constant real* restrict arg_alpha,
                    const __constant real* restrict arg_beta,
                    const __global realM* restrict agm,
                    const __global realN* restrict bgm,
                    __global realM* cgm) {

  const real alpha = arg_alpha[0];
  const real beta = arg_beta[0];

  // Allocates workgroup-private memory (local memory)
  #if SA == 1
    __local realM alm[KWG << ((MWG_SHIFT - VWM_SHIFT) >0 ? (MWG_SHIFT - VWM_SHIFT) : 0)];
  #endif
  #if SB == 1
    __local realN blm[KWG << ((NWG_SHIFT - VWN_SHIFT) >0 ? (NWG_SHIFT - VWN_SHIFT) : 0)];
  #endif

  // Computes the matrix-multiplication and stores the result in register memory
  realM cpm[NWI][MWI >> VWM_SHIFT];
  #if SA == 1 && SB == 1
    XgemmBody(kSizeM, kSizeN, kSizeK, agm, bgm, cgm, cpm, alm, blm);
  #elif SA == 1
    XgemmBody(kSizeM, kSizeN, kSizeK, agm, bgm, cgm, cpm, alm);
  #elif SB == 1
    XgemmBody(kSizeM, kSizeN, kSizeK, agm, bgm, cgm, cpm, blm);
  #else
    XgemmBody(kSizeM, kSizeN, kSizeK, agm, bgm, cgm, cpm);
  #endif

  // Stores an MWG * NWG tile of results and performs the multiplication with alpha and beta
  StoreResults(cgm, cpm, kSizeM, alpha, beta);
}

#endif
// =================================================================================================

// End of the C++11 raw string literal
)"

// =================================================================================================

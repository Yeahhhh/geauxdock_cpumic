#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "geauxdock.h"
#include "size.h"
#include "toggle.h"




#if TARGET_DEVICE == TARGET_MIC
__attribute__ ((target (mic)))
#endif
static float
MyRand_d (int seed)
{
    unsigned int s = time (NULL) + seed;
    float randdd = (float) rand_r (&s) / RAND_MAX;
    //float randdd = 0.002f;
    return randdd;
}


#define MYRAND (MyRand_d (r))
//#define MYRAND 0.334f






// Intel MKL VSL Random Number Generator
//#include <mkl_vsl.h>


#if TARGET_DEVICE == TARGET_MIC
__attribute__ ((target (mic)))
#endif
static float
combine_linear (float *e, const float *a_para, const float *b_para, const float *w_para)
{

    for (int i = 0; i < MAXWEI - 1; ++i)
        e[i] = a_para[i] * e[i] + b_para[i];
    float etotal = 0.0f;

    for (int i = 0; i < MAXWEI - 1; ++i)
        etotal += w_para[i] * e[i];

    e[MAXWEI - 1] = etotal;

    return etotal;
}



#if TARGET_DEVICE == TARGET_MIC
__attribute__ ((target (mic)))
#endif
static float
combine_bayes (float *e, const float *para)
{
    float etotal = 0.00f;
    return etotal;
}






#if TARGET_DEVICE == TARGET_MIC
__attribute__ ((target (mic)))
#endif
void
print_energies (float *e)
{
    printf ("k1  ");
    for (int i = 0; i < MAXWEI; ++i)
        printf ("%.8f ", e[i]);
    printf ("\n");

    printf ("ko    ");
    for (int i = 0; i < MAXWEI; ++i)
        printf ("%.4f ", e[i]);
    printf ("\n");
}







void
MonteCarlo_d (Complex * complex, Record * rec, const int s1, const int s2max)
{

#if TARGET_DEVICE == TARGET_MIC
    const int reps = complex->rep_end - complex->rep_begin + 1;
#pragma offload target(mic) \
    inout (complex:length(1) alloc_if(1) free_if(1)), \
    out (rec:length(reps) alloc_if(1) free_if(1))
#endif

#if TARGET_DEVICE == TARGET_CPU
#pragma omp parallel for
#endif
#if TARGET_DEVICE == TARGET_MIC
#pragma omp parallel for schedule (dynamic,1)
#endif

    for (int r = complex->rep_begin; r <= complex->rep_end; ++r) {

        // constant
        // pointers
        ReplicaMC *rep = &complex->replica[r];
        const Ligand *const lig = &complex->lig[rep->idx_lig];
        const Protein *const prt = &complex->prt[rep->idx_prt];
        const Psp *const psp = &complex->psp;
        const Kde *const kde = &complex->kde;
        const Mcs_ELL *const mcs_ell = &complex->mcs_ell;
        const EnePara *const enepara = &complex->enepara;



        // temporary write-read
        float move_scale[6] MYALIGN; // translation x y z, rotation x y z
        float movematrix[6] MYALIGN; // translation x y z, rotation x y z
        float rot[3][3] MYALIGN; // rotz roty rotx

        // not sorted
        float lig_x2[MAXLIG] MYALIGN;
        float lig_y2[MAXLIG] MYALIGN;
        float lig_z2[MAXLIG] MYALIGN;

        // sorted by t
        float lig_x3[MAXLIG] MYALIGN;
        float lig_y3[MAXLIG] MYALIGN;
        float lig_z3[MAXLIG] MYALIGN;


        // kde is sorted by t
        int kde_begin_idx[MAXLIG] MYALIGN;
        int kde_end_idx[MAXLIG] MYALIGN;


        // constant
        // ligand
        int   lig_t3[MAXLIG] MYALIGN;
        float lig_c3[MAXLIG] MYALIGN;
        float lig_center[3] MYALIGN;


        // constant
        // protein
        float prt_pocket_center[3] MYALIGN;


        // constant
        // mcs
        int mcs_ncol[MAX_MCS_ROW] MYALIGN; // row length in the sparse matrix representation
        float mcs_tcc[MAX_MCS_ROW] MYALIGN;


        // constant
        // enepara
        float enepara_p1a[MAXTP2][MAXTP1] MYALIGN;
        float enepara_p2a[MAXTP2][MAXTP1] MYALIGN;
        float enepara_pmf0[MAXTP2][MAXTP1] MYALIGN;
        float enepara_pmf1[MAXTP2][MAXTP1] MYALIGN;
        float enepara_hdb0[MAXTP2][MAXTP1] MYALIGN;
        float enepara_hdb1[MAXTP2][MAXTP1] MYALIGN;
        float enepara_hpl0[MAXTP2] MYALIGN;
        float enepara_hpl1[MAXTP2] MYALIGN;
        float enepara_hpl2[MAXTP2] MYALIGN;
        float enepara_a_para[MAXWEI] MYALIGN;
        float enepara_b_para[MAXWEI] MYALIGN;
        float enepara_w[MAXWEI] MYALIGN;
        float enepara_lj0;
        float enepara_lj1;
        float enepara_el1;
        float enepara_el0;
        float enepara_a1;
        float enepara_b1;
        float enepara_kde2;
        float enepara_kde3_inv;


        // constant
        // scalars
        float sqrt_2_pi_inv;
        int lig_natom, prt_npoint, kde_npoint, mcs_nrow;


        // temporary write-read vars
        int is_accept_s;




        for (int i = 0; i < 6; ++i)
            move_scale[i] = complex->mcpara.move_scale[i];

        for (int l = 0; l < MAXTP2; ++l) {
            for (int p = 0; p < MAXTP1; ++p) {
                enepara_p1a[l][p] = enepara->p1a[l][p];
                enepara_p2a[l][p] = enepara->p2a[l][p];
                enepara_pmf0[l][p] = enepara->pmf0[l][p];
                enepara_pmf1[l][p] = enepara->pmf1[l][p];
                enepara_hdb0[l][p] = enepara->hdb0[l][p];
                enepara_hdb1[l][p] = enepara->hdb1[l][p];
            }
        }

        for (int l = 0; l < MAXTP2; ++l) {
            enepara_hpl0[l] = enepara->hpl0[l];
            enepara_hpl1[l] = enepara->hpl1[l];
            enepara_hpl2[l] = enepara->hpl2[l];
        }

        for (int i = 0; i < MAXWEI - 1; ++i) {
            enepara_a_para[i] = enepara->a_para[i];
            enepara_b_para[i] = enepara->b_para[i];
            enepara_w[i] = enepara->w[i];
        }

        enepara_lj0 = enepara->lj0;
        enepara_lj1 = enepara->lj1;
        enepara_el1 = enepara->el1;
        enepara_el0 = enepara->el0;
        enepara_a1 = enepara->a1;
        enepara_b1 = enepara->b1;
        enepara_kde2 = enepara->kde2;
        enepara_kde3_inv = 1.0f / enepara->kde3;

        sqrt_2_pi_inv = -1.0f / sqrtf (2.0f * PI);
        lig_natom = lig->lig_natom;
        prt_npoint = prt->prt_npoint;
        kde_npoint = complex->size.kde_npoint;
        mcs_nrow = complex->size.mcs_nrow;
        is_accept_s = rep->is_accept;

        rec[r - complex->rep_begin].next_entry = 0;	// reset the record's next_entry


        for (int l = 0; l < lig_natom; ++l) {
            lig_t3[l] = lig->t3[l];
            lig_c3[l] = lig->c3[l];
            kde_begin_idx[l] = lig->kde_begin_idx[l];
            kde_end_idx[l] = lig->kde_end_idx[l];
        }

        for (int i = 0; i < 3; ++i) {
            lig_center[i] = lig->center[i];
            prt_pocket_center[i] = prt->pocket_center[i];
        }

        for (int m = 0; m < mcs_nrow; ++m) {
            mcs_ncol[m] = mcs_ell->ncol[m];
            mcs_tcc[m] = mcs_ell->tcc[m];
        }







        for (int s2 = 0; s2 < s2max; ++s2) {


            /////////////////////////////////////////////////////////////////////////////
            // record old states
            if (is_accept_s == 1) {
                rep->step = s1 + s2;

                const int rr = r - complex->rep_begin;
                const int next_entry = rec[rr].next_entry;
                rec[rr].replica[next_entry] = *rep;
                rec[rr].next_entry = next_entry + 1;
            }






            /////////////////////////////////////////////////////////////////////////////
            // perturb

            //#pragma unroll (6)
            for (int i = 0; i < 6; ++i) {
#if IS_AWAY == 0
                const float fixed_var = 0.0f;
#elif IS_AWAY == 1
                const float fixed_var = 44.5f;
#endif
                float moveamount = s2max != 1 ? MYRAND : fixed_var;
                movematrix[i] = move_scale[i] * moveamount + rep->movematrix[i];
            }


            // http://en.wikipedia.org/wiki/Euler_angles
            // http://upload.wikimedia.org/math/e/9/c/e9cf817bce9c1780216921cd93233459.png
            // http://upload.wikimedia.org/math/f/4/e/f4e55dc2c9581007648967d29b15121e.png
            const float sin1 = sinf (movematrix[3]);
            const float cos1 = cosf (movematrix[3]);
            const float sin2 = sinf (movematrix[4]);
            const float cos2 = cosf (movematrix[4]);
            const float sin3 = sinf (movematrix[5]);
            const float cos3 = cosf (movematrix[5]);
            rot[0][0] = cos1 * cos2;
            rot[0][1] = cos1 * sin2 * sin3 - cos3 * sin1;
            rot[0][2] = sin1 * sin3 + cos1 * cos3 * sin2;
            rot[1][0] = cos2 * sin1;
            rot[1][1] = cos1 * cos3 + sin1 * sin2 * sin3;
            rot[1][2] = cos3 * sin1 * sin2 - cos1 * sin3;
            rot[2][0] = -1 * sin2;
            rot[2][1] = cos2 * sin3;
            rot[2][2] = cos2 * cos3;


            // rotation, translation, coordinate system transformation
            for (int l = 0; l < lig_natom; ++l) {
                const float x2 = lig->x2[l];
                const float y2 = lig->y2[l];
                const float z2 = lig->z2[l];
                const float x3 = lig->x3[l];
                const float y3 = lig->y3[l];
                const float z3 = lig->z3[l];
                lig_x2[l] = rot[0][0] * x2 + rot[0][1] * y2 + rot[0][2] * z2 + movematrix[0] + lig_center[0];
                lig_y2[l] = rot[1][0] * x2 + rot[1][1] * y2 + rot[1][2] * z2 + movematrix[1] + lig_center[1];
                lig_z2[l] = rot[2][0] * x2 + rot[2][1] * y2 + rot[2][2] * z2 + movematrix[2] + lig_center[2];
                lig_x3[l] = rot[0][0] * x3 + rot[0][1] * y3 + rot[0][2] * z3 + movematrix[0] + lig_center[0];
                lig_y3[l] = rot[1][0] * x3 + rot[1][1] * y3 + rot[1][2] * z3 + movematrix[1] + lig_center[1];
                lig_z3[l] = rot[2][0] * x3 + rot[2][1] * y3 + rot[2][2] * z3 + movematrix[2] + lig_center[2];
            }





            /////////////////////////////////////////////////////////////////////////////
            // calcenergy

            float evdw = 0.0f;
            float eele = 0.0f;
            float epmf = 0.0f;
            float epsp = 0.0f;
            float ehdb = 0.0f;
            float ehpc = 0.0f;


            for (int l = 0; l < lig_natom; ++l) { // lig loop, ~30
                float ehpc1 = 0.0f;
                const int lig__t = lig_t3[l];

#pragma simd reduction(+:evdw, eele, epmf, epsp, ehdb, ehpc, ehpc1)
                for (int p = 0; p < prt_npoint; ++p) { // prt loop, ~300
                    const int prt__t = prt->t[p];

                    const float dx = lig_x3[l] - prt->x[p];
                    const float dy = lig_y3[l] - prt->y[p];
                    const float dz = lig_z3[l] - prt->z[p];
                    const float dst_pow2 = dx * dx + dy * dy + dz * dz;
                    const float dst_pow4 = dst_pow2 * dst_pow2;
                    const float dst = sqrtf (dst_pow2);


                    /* hydrophobic potential */
                    if (prt->cdc[p] == 1 && dst_pow2 <= 81.0f) {
                        ehpc1 += prt->hpp[p] *
                            (1.0f - (3.5f / 81.0f * dst_pow2 -
                                     4.5f / 81.0f / 81.0f * dst_pow4 +
                                     2.5f / 81.0f / 81.0f / 81.0f * dst_pow4 * dst_pow2 -
                                     0.5f / 81.0f / 81.0f / 81.0f / 81.0f * dst_pow4 * dst_pow4));
                    }

                    /* L-J potential */
                    const float p1 =
                        enepara_p1a[lig__t][prt__t] / (dst_pow4 * dst_pow4 * dst);
                    const float p2 =
                        enepara_p2a[lig__t][prt__t] / (dst_pow4 * dst_pow2);
                    const float p4 =
                        p1 * enepara_lj0 * (1.0f + enepara_lj1 * dst_pow2) + 1.0f;
                    evdw += (p1 - p2) / p4;


                    /* electrostatic potential */
                    const float s1 = enepara_el1 * dst;
                    float g1;
                    if (s1 < 1.0f)
                        g1 = enepara_el0 + enepara_a1 * s1 * s1 + enepara_b1 * s1 * s1 * s1;
                    else
                        g1 = 1.0f / s1;
                    eele += lig_c3[l] * prt->ele[p] * g1;


                    /* contact potential */
                    const float dst_minus_pmf0 = dst - enepara_pmf0[lig__t][prt__t];
                    epmf +=
                        enepara_pmf1[lig__t][prt__t] /
                        (1.0f + expf ((-0.5f * dst + 6.0f) * dst_minus_pmf0));


                    /* pocket-specific potential */
                    // the senmatics do not match with the original program:
                    // if (found psp[][])
                    //   accumulate to epsp
                    // else
                    //   do nothing
                    if (prt->c[p] == 2 && dst_minus_pmf0 <= 0) {
                        const int i1 = prt->seq3r[p];
                        epsp += psp->psp[lig__t][i1];
                    }

                    /* hydrogen bond potential */
                    const float hdb0 = enepara_hdb0[lig__t][prt__t];
                    if (hdb0 > 0.1f) {
                        const float hdb1 = enepara_hdb1[lig__t][prt__t];
                        const float hdb3 = (dst - hdb0) * hdb1;
                        ehdb += hdb1 * expf (-0.5f * hdb3 * hdb3);
                    }

                } // prt loop


                /* hydrophobic restraits */
                const float hpc2 = (ehpc1 - enepara_hpl0[lig__t]) / enepara_hpl1[lig__t];
                ehpc += 0.5f * hpc2 * hpc2 - enepara_hpl2[lig__t];
            } // lig loop

            float e_s[MAXWEI] MYALIGN;
            e_s[0] = evdw; // 0 - vdw
            e_s[1] = eele; // 1 - ele
            e_s[2] = epmf; // 2 - pmf (CP)
            e_s[3] = epsp; // 3 - psp (PS CP)
            e_s[4] = ehdb * sqrt_2_pi_inv; // 4 - hdb (HB)
            e_s[5] = ehpc; // 5 - hpc (HP)





            // KDE
            float ekde = 0.0f;
            for (int l = 0; l < lig_natom; ++l) { // lig loop, ~30
                float ekde1 = 0.0f;
                const int lig__t = lig_t3[l];

                const int begin = kde_begin_idx[l];
                const int end = kde_end_idx[l];

                for (int k = begin; k < end; ++k) { // kde loop, ~400
                    const float dx = lig_x3[l] - kde->x[k];
                    const float dy = lig_y3[l] - kde->y[k];
                    const float dz = lig_z3[l] - kde->z[k];
                    const float kde_dst_pow2 = dx * dx + dy * dy + dz * dz;
                    ekde1 += expf (enepara_kde2 * kde_dst_pow2);
                }

                float kde_sz = (float) (end - begin);
                ekde += ekde1 / kde_sz;
            }
            e_s[6] = ekde * enepara_kde3_inv;




            // MCS, ELLPACK format sparse matrix
            float elhm = 0.0f;
            for (int m = 0; m < mcs_nrow; ++m) {
                float elhm1 = 0.0f;
                const int ncol = mcs_ncol[m];
#pragma simd reduction(+:elhm1)
#pragma vector aligned
                for (int i = 0; i < ncol; ++i) {
                    const int l = mcs_ell->i[m][i];
                    const float dx = lig_x2[l] - mcs_ell->x[m][i];
                    const float dy = lig_y2[l] - mcs_ell->y[m][i];
                    const float dz = lig_z2[l] - mcs_ell->z[m][i];
                    elhm1 += dx * dx + dy * dy + dz * dz;
                }
                elhm += mcs_tcc[m] * sqrtf (elhm1 / (float) ncol);
            }
            e_s[7] = logf (elhm / mcs_nrow);




            // EDST
            float dst = 0;
#pragma simd reduction(+:dst)
            for (int i = 0; i < 3; ++i) {
                float d = lig_center[i] + movematrix[i] - prt_pocket_center[i];
                dst += d * d;
            }
            e_s[8] = sqrtf (dst);






            for (int i = 0; i < 7; ++i)
                e_s[i] = e_s[i] / lig_natom;



            const float etotal = combine_linear (e_s, enepara_a_para, enepara_b_para, enepara_w);

#if 0
            if (r == 0 && s2max == 1)
                print_energies (e_s);
#endif




            ////////////////////////////////////////////////////////////////////////
            // accept
            const float delta_energy = etotal - rep->energy[MAXWEI - 1];
            const float beta = complex->temp[rep->idx_tmp].minus_beta;
            float rand = s2max != 1 ? MYRAND : 0.0f; // force to accept if s2max == 1
            is_accept_s = (rand < expf (delta_energy * beta));	// mybeta < 0

            if (is_accept_s == 1) {
                for (int i = 0; i < MAXWEI; ++i)
                    rep->energy[i] = e_s[i];
                for (int i = 0; i < 6; ++i)
                    rep->movematrix[i] = movematrix[i];
            }


        } // s2 loop


        rep->is_accept = is_accept_s;

    } // replica loop

}




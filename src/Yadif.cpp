/*
**   Yadif for VapourSynth
**   Copyright (C) 2016  Kusaki
**
**   Yadif C-plugin for Avisynth 2.5 - Yet Another DeInterlacing Filter
**   Copyright (C) 2007  Alexander G. Balakhnin aka Fizick  http://avisynth.org.ru
**   Port of YADIF filter from MPlayer
**   Copyright (C) 2006  Michael Niedermayer <michaelni@gmx.at>
**
**   This program is free software; you can redistribute it and/or modify
**   it under the terms of the GNU General Public License as published by
**   the Free Software Foundation; either version 2 of the License, or
**   (at your option) any later version.
**
**   This program is distributed in the hope that it will be useful,
**   but WITHOUT ANY WARRANTY; without even the implied warranty of
**   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**   GNU General Public License for more details.
**
**   You should have received a copy of the GNU General Public License along
**   with this program; if not, write to the Free Software Foundation, Inc.,
**   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <vapoursynth/VapourSynth.h>
#include <vapoursynth/VSHelper.h>

#ifdef VS_TARGET_CPU_X86
#define MAX_VECTOR_SIZE 256
#include "vectorclass/vectorclass.h"
#endif

#define YADIF_VERSION "1"

struct YadifData
{
    VSNodeRef *node;
    VSVideoInfo vi;
    const VSVideoInfo *vi_src;
    
    int order;
    int field;
    int mode;
};

static void VS_CC yadifInit(VSMap *in, VSMap *out, void **instanceData, VSNode *node, VSCore *core, const VSAPI *vsapi)
{
    YadifData *d = static_cast<YadifData *>(*instanceData);
    vsapi->setVideoInfo(&d->vi, 1, node);
}

#ifdef VS_TARGET_CPU_X86

template<typename T1, typename T2, typename T3, typename V1, typename V2, typename V2b, int vectorSize>
static void filter_line(const T1 *prev, const T1 *cur, const T1 *next, T1 * VS_RESTRICT dst,
                        T2 *dp, T2 *diffp, T2 *dir0p, T2 *dir1p, T2 *dir2p, T2 *dir3p, T2 *dir4p, T3 *pred_dirp,
                        const int w, const int stride, const int parity, const int mode)
{
    const T1 *prev2 = parity ? prev : cur;
    const T1 *next2 = parity ? cur  : next;
    const V2 vzero(0), vone(1), vtwo(2), vthree(3), vfour(4);
    T1 *prednp = vs_aligned_malloc<T1>(16, 32);
    T1 *predpp = vs_aligned_malloc<T1>(16, 32);

    for (int x = 0; x < w; x += vectorSize)
    {
        const V1 vprevn1 = V1().load_a(prev - stride + x);
        const V1 vprevp1 = V1().load_a(prev + stride + x);
        const V1 vprev2 = V1().load_a(prev2 + x);
        const V1 vprev2n2 = V1().load_a(prev2 - 2 * stride + x);
        const V1 vprev2p2 = V1().load_a(prev2 + 2 * stride + x);
        const V1 vcurn1 = V1().load_a(cur - stride + x);
        const V1 vcurp1 = V1().load_a(cur + stride + x);
        const V1 vnextn1 = V1().load_a(next - stride + x);
        const V1 vnextp1 = V1().load_a(next + stride + x);
        const V1 vnext2 = V1().load_a(next2 + x);
        const V1 vnext2n2 = V1().load_a(next2 - 2 * stride + x);
        const V1 vnext2p2 = V1().load_a(next2 + 2 * stride + x);
        
        const V2 c = V2(extend_low(vcurn1), extend_high(vcurn1));
        const V2 d = (V2(extend_low(vprev2), extend_high(vprev2)) + V2(extend_low(vnext2), extend_high(vnext2))) >> 1;
        const V2 e = V2(extend_low(vcurp1), extend_high(vcurp1));
        const V2 tdiff0 = abs(V2(extend_low(vprev2), extend_high(vprev2)) - V2(extend_low(vnext2), extend_high(vnext2))) >> 1;
        const V2 tdiff1 = (abs(V2(extend_low(vprevn1), extend_high(vprevn1)) - c) + abs(V2(extend_low(vprevp1), extend_high(vprevp1)) - e)) >> 1;
        const V2 tdiff2 = (abs(V2(extend_low(vnextn1), extend_high(vnextn1)) - c) + abs(V2(extend_low(vnextp1), extend_high(vnextp1)) - e)) >> 1;
        V2 diff = max(max(tdiff0, tdiff1), tdiff2);

        // check
        if (mode < 2)
        {
            const V2 b = (V2(extend_low(vprev2n2), extend_high(vprev2n2)) + V2(extend_low(vnext2n2), extend_high(vnext2n2))) >> 1;
            const V2 f = (V2(extend_low(vprev2p2), extend_high(vprev2p2)) + V2(extend_low(vnext2p2), extend_high(vnext2p2))) >> 1;
            const V2 maxs = max(max(d - e, d - c), min(b - c, f - e));
            const V2 mins = min(min(d - e, d - c), max(b - c, f - e));
            diff = max(max(diff, mins), -maxs);
        }
        
        d.store_a(dp + x);
        diff.store_a(diffp + x);
        
        const V1 vdir0n1 = V1().load(cur - stride + x - 2);
        const V1 vdir0p1 = V1().load(cur + stride + x + 2);
        const V1 vdir1n1 = V1().load(cur - stride + x - 1);
        const V1 vdir1p1 = V1().load(cur + stride + x + 1);
        const V1 vdir3n1 = V1().load(cur - stride + x + 1);
        const V1 vdir3p1 = V1().load(cur + stride + x - 1);
        const V1 vdir4n1 = V1().load(cur - stride + x + 2);
        const V1 vdir4p1 = V1().load(cur + stride + x - 2);
        
        const V2 vdir0 = abs(V2(extend_low(vdir0n1), extend_high(vdir0n1)) - V2(extend_low(vdir0p1), extend_high(vdir0p1)));
        const V2 vdir1 = abs(V2(extend_low(vdir1n1), extend_high(vdir1n1)) - V2(extend_low(vdir1p1), extend_high(vdir1p1)));
        const V2 vdir2 = abs(c - e);
        const V2 vdir3 = abs(V2(extend_low(vdir3n1), extend_high(vdir3n1)) - V2(extend_low(vdir3p1), extend_high(vdir3p1)));
        const V2 vdir4 = abs(V2(extend_low(vdir4n1), extend_high(vdir4n1)) - V2(extend_low(vdir4p1), extend_high(vdir4p1)));
        
        vdir0.store_a(dir0p + x);
        vdir1.store_a(dir1p + x);
        vdir2.store_a(dir2p + x);
        vdir3.store_a(dir3p + x);
        vdir4.store_a(dir4p + x);
    }
    
    for (int x = 0; x < w; x += vectorSize)
    {
        const V2 vdir0l = V2().load(dir0p + x - 1);
        const V2 vdir0m = V2().load_a(dir0p + x);
        const V2 vdir0r = V2().load(dir0p + x + 1);
        const V2 vdir1l = V2().load(dir1p + x - 1);
        const V2 vdir1m = V2().load_a(dir1p + x);
        const V2 vdir1r = V2().load(dir1p + x + 1);
        const V2 vdir2l = V2().load(dir2p + x - 1);
        const V2 vdir2m = V2().load_a(dir2p + x);
        const V2 vdir2r = V2().load(dir2p + x + 1);
        const V2 vdir3l = V2().load(dir3p + x - 1);
        const V2 vdir3m = V2().load_a(dir3p + x);
        const V2 vdir3r = V2().load(dir3p + x + 1);
        const V2 vdir4l = V2().load(dir4p + x - 1);
        const V2 vdir4m = V2().load_a(dir4p + x);
        const V2 vdir4r = V2().load(dir4p + x + 1);
        
        const V2 vscore_dir0 = vdir0l + vdir0m + vdir0r;
        const V2 vscore_dir1 = vdir1l + vdir1m + vdir1r;
        const V2 vscore_dir2 = vdir2l + vdir2m + vdir2r;
        const V2 vscore_dir3 = vdir3l + vdir3m + vdir3r;
        const V2 vscore_dir4 = vdir4l + vdir4m + vdir4r;
        
        
        V2b vscoreb = vscore_dir1 < vscore_dir0;
        V2 vscore = select(vscoreb, vscore_dir1, vscore_dir0);
        V2 vpred_dir = select(vscoreb, vone, vzero);
        
        vscoreb = vscore_dir2 < vscore;
        vscore = select(vscoreb, vscore_dir2, vscore);
        vpred_dir = select(vscoreb, vtwo, vpred_dir);
        
        vscoreb = vscore_dir3 < vscore;
        vscore = select(vscoreb, vscore_dir3, vscore);
        vpred_dir = select(vscoreb, vthree, vpred_dir);
        
        vscoreb = vscore_dir4 < vscore;
        vscore = select(vscoreb, vscore_dir4, vscore);
        vpred_dir = select(vscoreb, vfour, vpred_dir);
        
        // spatial_pred = (cur[-stride+x-2] + cur[stride+x+2]) >> 1; - dir0
        vpred_dir.store_a(pred_dirp);
        for (int i = 0; i < vectorSize; i++)
        {
            prednp[i] = cur[-stride + x + i + pred_dirp[i] - 2];
        }

        for (int i = 0; i < vectorSize; i++)
        {
            predpp[i] = cur[stride + x + i + 2 - pred_dirp[i]];
        }
        
        const V1 vdirnn1 = V1().load_a(prednp);
        const V1 vdirnp1 = V1().load_a(predpp);
        const V2 vspatial_pred = (V2(extend_low(vdirnn1), extend_high(vdirnn1)) + V2(extend_low(vdirnp1), extend_high(vdirnp1))) >> 1;
        
        const V2 d = V2().load_a(dp + x);
        const V2 diff = V2().load_a(diffp + x);
        const V2 vpred = max(min(vspatial_pred, d + diff), d - diff);
        V1(compress(vpred.get_low(), vpred.get_high())).store_a(dst + x);
    }
    
    vs_aligned_free(prednp);
    vs_aligned_free(predpp);
}

template<>
void filter_line<float, float, int32_t, Vec8f, Vec8f, Vec8fb, 8>(const float *prev, const float *cur, const float *next, float * VS_RESTRICT dst,
                                                                 float *dp, float *diffp, float *dir0p, float *dir1p, float *dir2p, float *dir3p, float *dir4p, int32_t *pred_dirp,
                                                                 const int w, const int stride, const int parity, const int mode)
{
    const float *prev2 = parity ? prev : cur;
    const float *next2 = parity ? cur  : next;
    const Vec8f vpointFive(0.5f);
    const Vec8i vzero(0), vone(1), vtwo(2), vthree(3), vfour(4);
    
    for (int x = 0; x < w; x += 8)
    {
        const Vec8f vprevn1 = Vec8f().load_a(prev - stride + x);
        const Vec8f vprevp1 = Vec8f().load_a(prev + stride + x);
        const Vec8f vprev2 = Vec8f().load_a(prev2 + x);
        const Vec8f vprev2n2 = Vec8f().load_a(prev2 - 2 * stride + x);
        const Vec8f vprev2p2 = Vec8f().load_a(prev2 + 2 * stride + x);
        const Vec8f vnextn1 = Vec8f().load_a(next - stride + x);
        const Vec8f vnextp1 = Vec8f().load_a(next + stride + x);
        const Vec8f vnext2 = Vec8f().load_a(next2 + x);
        const Vec8f vnext2n2 = Vec8f().load_a(next2 - 2 * stride + x);
        const Vec8f vnext2p2 = Vec8f().load_a(next2 + 2 * stride + x);
                
        const Vec8f c = Vec8f().load_a(cur - stride + x);
        const Vec8f d = (vprev2 + vnext2) * vpointFive;
        const Vec8f e = Vec8f().load_a(cur + stride + x);
        const Vec8f tdiff0 = abs(vprev2 - vnext2) * vpointFive;
        const Vec8f tdiff1 = (abs(vprevn1 - c) + abs(vprevp1 - e)) * vpointFive;
        const Vec8f tdiff2 = (abs(vnextn1 - c) + abs(vnextp1 - e)) * vpointFive;
        Vec8f diff = max(max(tdiff0, tdiff1), tdiff2);

        // check
        if (mode < 2)
        {
            const Vec8f b = (vprev2n2 + vnext2n2) * vpointFive;
            const Vec8f f = (vprev2p2 + vnext2p2) * vpointFive;
            const Vec8f maxs = max(max(d - e, d - c), min(b - c, f - e));
            const Vec8f mins = min(min(d - e, d - c), max(b - c, f - e));
            diff = max(max(diff, mins), -maxs);
        }
        
        d.store_a(dp + x);
        diff.store_a(diffp + x);
        
        const Vec8f vdir0n1 = Vec8f().load(cur - stride + x - 2);
        const Vec8f vdir0p1 = Vec8f().load(cur + stride + x + 2);
        const Vec8f vdir1n1 = Vec8f().load(cur - stride + x - 1);
        const Vec8f vdir1p1 = Vec8f().load(cur + stride + x + 1);
        const Vec8f vdir3n1 = Vec8f().load(cur - stride + x + 1);
        const Vec8f vdir3p1 = Vec8f().load(cur + stride + x - 1);
        const Vec8f vdir4n1 = Vec8f().load(cur - stride + x + 2);
        const Vec8f vdir4p1 = Vec8f().load(cur + stride + x - 2);
        
        const Vec8f vdir0 = abs(vdir0n1 - vdir0p1);
        const Vec8f vdir1 = abs(vdir1n1 - vdir1p1);
        const Vec8f vdir2 = abs(c - e);
        const Vec8f vdir3 = abs(vdir3n1 - vdir3p1);
        const Vec8f vdir4 = abs(vdir4n1 - vdir4p1);
        
        vdir0.store_a(dir0p + x);
        vdir1.store_a(dir1p + x);
        vdir2.store_a(dir2p + x);
        vdir3.store_a(dir3p + x);
        vdir4.store_a(dir4p + x);
    }
    
    for (int x = 0; x < w; x += 8)
    {
        const Vec8f vdir0l = Vec8f().load(dir0p + x - 1);
        const Vec8f vdir0m = Vec8f().load_a(dir0p + x);
        const Vec8f vdir0r = Vec8f().load(dir0p + x + 1);
        const Vec8f vdir1l = Vec8f().load(dir1p + x - 1);
        const Vec8f vdir1m = Vec8f().load_a(dir1p + x);
        const Vec8f vdir1r = Vec8f().load(dir1p + x + 1);
        const Vec8f vdir2l = Vec8f().load(dir2p + x - 1);
        const Vec8f vdir2m = Vec8f().load_a(dir2p + x);
        const Vec8f vdir2r = Vec8f().load(dir2p + x + 1);
        const Vec8f vdir3l = Vec8f().load(dir3p + x - 1);
        const Vec8f vdir3m = Vec8f().load_a(dir3p + x);
        const Vec8f vdir3r = Vec8f().load(dir3p + x + 1);
        const Vec8f vdir4l = Vec8f().load(dir4p + x - 1);
        const Vec8f vdir4m = Vec8f().load_a(dir4p + x);
        const Vec8f vdir4r = Vec8f().load(dir4p + x + 1);
        
        const Vec8f vscore_dir0 = vdir0l + vdir0m + vdir0r;
        const Vec8f vscore_dir1 = vdir1l + vdir1m + vdir1r;
        const Vec8f vscore_dir2 = vdir2l + vdir2m + vdir2r;
        const Vec8f vscore_dir3 = vdir3l + vdir3m + vdir3r;
        const Vec8f vscore_dir4 = vdir4l + vdir4m + vdir4r;
        
        
        Vec8fb vscoreb = vscore_dir1 < vscore_dir0;
        Vec8f vscore = select(vscoreb, vscore_dir1, vscore_dir0);
        Vec8i vpred_dir = select(Vec8ib(vscoreb), vone, vzero);
        
        vscoreb = vscore_dir2 < vscore;
        vscore = select(vscoreb, vscore_dir2, vscore);
        vpred_dir = select(Vec8ib(vscoreb), vtwo, vpred_dir);
        
        vscoreb = vscore_dir3 < vscore;
        vscore = select(vscoreb, vscore_dir3, vscore);
        vpred_dir = select(Vec8ib(vscoreb), vthree, vpred_dir);
        
        vscoreb = vscore_dir4 < vscore;
        vscore = select(vscoreb, vscore_dir4, vscore);
        vpred_dir = select(Vec8ib(vscoreb), vfour, vpred_dir);
        
        // spatial_pred = (cur[-stride+x-2] + cur[stride+x+2]) * 0.5f; - dir0
        vpred_dir.store_a(pred_dirp);
        const Vec8f vdirnn1(cur[-stride + x + 0 + pred_dirp[0] - 2], cur[-stride + x + 1 + pred_dirp[1] - 2],
                            cur[-stride + x + 2 + pred_dirp[2] - 2], cur[-stride + x + 3 + pred_dirp[3] - 2],
                            cur[-stride + x + 4 + pred_dirp[4] - 2], cur[-stride + x + 5 + pred_dirp[5] - 2],
                            cur[-stride + x + 6 + pred_dirp[6] - 2], cur[-stride + x + 7 + pred_dirp[7] - 2]);
        const Vec8f vdirnp1(cur[stride + x + 0 + 2 - pred_dirp[0]], cur[stride + x + 1 + 2 - pred_dirp[1]],
                            cur[stride + x + 2 + 2 - pred_dirp[2]], cur[stride + x + 3 + 2 - pred_dirp[3]],
                            cur[stride + x + 4 + 2 - pred_dirp[4]], cur[stride + x + 5 + 2 - pred_dirp[5]],
                            cur[stride + x + 6 + 2 - pred_dirp[6]], cur[stride + x + 7 + 2 - pred_dirp[7]]);
        const Vec8f vspatial_pred = (vdirnn1 + vdirnp1) * vpointFive;
        
        const Vec8f d = Vec8f().load_a(dp + x);
        const Vec8f diff = Vec8f().load_a(diffp + x);
        const Vec8f vpred = max(min(vspatial_pred, d + diff), d - diff);
        vpred.store_a(dst + x);
    }
}

#else

template<typename T1>
static void filter_line(const T1 *prev, const T1 *cur, const T1 *next, T1 * VS_RESTRICT dst,
                        const int w, const int stride, const int parity, const int mode)
{
    const T1 *prev2 = parity ? prev : cur;
    const T1 *next2 = parity ? cur  : next;

    for (int x = 3; x < w - 3; x++)
    {
        const int c = cur[-stride+x];
        const int d = (prev2[x] + next2[x]) >> 1;
        const int e = cur[stride+x];
        const int tdiff0 = std::abs(prev2[x] - next2[x]) >> 1;
        const int tdiff1 = (std::abs(prev[-stride+x] - c) + std::abs(prev[stride+x] - e)) >> 1;
        const int tdiff2 = (std::abs(next[-stride+x] - c) + std::abs(next[stride+x] - e)) >> 1;
        int diff = std::max(std::max(tdiff0, tdiff1), tdiff2);

        // edi
        int spatial_pred = (c + e) >> 1;
        int spatial_score = std::abs(cur[-stride+x-1] - cur[stride+x-1]) + std::abs(c - e) + 
                            std::abs(cur[-stride+x+1] - cur[stride+x+1]) - 1;

        int score = std::abs(cur[-stride+x-2] - cur[stride+x]) +
                    std::abs(cur[-stride+x-1] - cur[stride+x+1]) +
                    std::abs(cur[-stride+x] - cur[stride+x+2]);
        if (score < spatial_score)
        {
            spatial_score = score;
            spatial_pred = (cur[-stride+x-1] + cur[stride+x+1]) >> 1;
            score = std::abs(cur[-stride+x-3] - cur[stride+x+1]) +
                    std::abs(cur[-stride+x-2] - cur[stride+x+2]) +
                    std::abs(cur[-stride+x-1] - cur[stride+x+3]);
            if (score < spatial_score)
            {
                spatial_score = score;
                spatial_pred = (cur[-stride+x-2] + cur[stride+x+2]) >> 1;
            }
        }

        score = std::abs(cur[-stride+x] - cur[stride+x-2]) +
                std::abs(cur[-stride+x+1] - cur[stride+x-1]) +
                std::abs(cur[-stride+x+2] - cur[stride+x]);
        if (score < spatial_score)
        {
            spatial_score = score;
            spatial_pred = (cur[-stride+x+1] + cur[stride+x-1]) >> 1;
            score = std::abs(cur[-stride+x+1] - cur[stride+x-3]) +
                    std::abs(cur[-stride+x+2] - cur[stride+x-2]) +
                    std::abs(cur[-stride+x+3] - cur[stride+x-1]);
            if (score < spatial_score)
            {
                spatial_score = score;
                spatial_pred = (cur[-stride+x+2] + cur[stride+x-2]) >> 1;
            }
        }

        // check
        if (mode < 2)
        {
            const int b = (prev2[-2*stride+x] + next2[-2*stride+x]) >> 1;
            const int f = (prev2[2*stride+x] + next2[2*stride+x]) >> 1;
            const int maxs = std::max(std::max(d - e, d - c), std::min(b - c, f - e));
            const int mins = std::min(std::min(d - e, d - c), std::max(b - c, f - e));
            diff = std::max(std::max(diff, mins), -maxs);
        }

        dst[x] = std::max(std::min<int>(spatial_pred, d + diff), d - diff);
    }
}

template<>
void filter_line<float>(const float *prev, const float *cur, const float *next, float * VS_RESTRICT dst,
                        const int w, const int stride, const int parity, const int mode)
{
    const float *prev2 = parity ? prev : cur;
    const float *next2 = parity ? cur  : next;

    for (int x = 3; x < w - 3; x++)
    {
        const float c = cur[-stride+x];
        const float d = (prev2[x] + next2[x]) * 0.5f;
        const float e = cur[stride+x];
        const float tdiff0 = std::abs(prev2[x] - next2[x]) * 0.5f;
        const float tdiff1 = (std::abs(prev[-stride+x] - c) + std::abs(prev[stride+x] - e)) * 0.5f;
        const float tdiff2 = (std::abs(next[-stride+x] - c) + std::abs(next[stride+x] - e)) * 0.5f;
        float diff = std::max(std::max(tdiff0, tdiff1), tdiff2);

        // edi
        float spatial_pred = (c + e) * 0.5f;
        float spatial_score = std::abs(cur[-stride+x-1] - cur[stride+x-1]) + std::abs(c - e) + 
                              std::abs(cur[-stride+x+1] - cur[stride+x+1]);

        float score = std::abs(cur[-stride+x-2] - cur[stride+x]) +
                      std::abs(cur[-stride+x-1] - cur[stride+x+1]) +
                      std::abs(cur[-stride+x] - cur[stride+x+2]);
        if (score < spatial_score)
        {
            spatial_score = score;
            spatial_pred = (cur[-stride+x-1] + cur[stride+x+1]) * 0.5f;
            score = std::abs(cur[-stride+x-3] - cur[stride+x+1]) +
                    std::abs(cur[-stride+x-2] - cur[stride+x+2]) +
                    std::abs(cur[-stride+x-1] - cur[stride+x+3]);
            if (score < spatial_score)
            {
                spatial_score = score;
                spatial_pred = (cur[-stride+x-2] + cur[stride+x+2]) * 0.5f;
            }
        }

        score = std::abs(cur[-stride+x] - cur[stride+x-2]) +
                std::abs(cur[-stride+x+1] - cur[stride+x-1]) +
                std::abs(cur[-stride+x+2] - cur[stride+x]);
        if (score < spatial_score)
        {
            spatial_score = score;
            spatial_pred = (cur[-stride+x+1] + cur[stride+x-1]) * 0.5f;
            score = std::abs(cur[-stride+x+1] - cur[stride+x-3]) +
                    std::abs(cur[-stride+x+2] - cur[stride+x-2]) +
                    std::abs(cur[-stride+x+3] - cur[stride+x-1]);
            if (score < spatial_score)
            {
                spatial_score = score;
                spatial_pred = (cur[-stride+x+2] + cur[stride+x-2]) * 0.5f;
            }
        }

        // check
        if (mode < 2)
        {
            const float b = (prev2[-2*stride+x] + next2[-2*stride+x]) * 0.5f;
            const float f = (prev2[2*stride+x] + next2[2*stride+x]) * 0.5f;
            const float maxs = std::max(std::max(d - e, d - c), std::min(b - c, f - e));
            const float mins = std::min(std::min(d - e, d - c), std::max(b - c, f - e));
            diff = std::max(std::max(diff, mins), -maxs);
        }

        dst[x] = std::max(std::min(spatial_pred, d + diff), d - diff);
    }
}

#endif

#ifdef VS_TARGET_CPU_X86

template<typename T1, typename V1, typename V2, int vectorSize>
static void interpolate(const T1 *cur0, const T1 *cur2, T1 * VS_RESTRICT dst, const int w)
{
    const V2 vone(1);
    
    for (int x = 0; x < w; x += vectorSize)
    {
        const V1 vcur0 = V1().load_a(cur0 + x);
        const V1 vcur2 = V1().load_a(cur2 + x);
        
        const V2 vpred = (V2(extend_low(vcur0), extend_high(vcur0)) + V2(extend_low(vcur2), extend_high(vcur2)) + vone) >> 1;
        V1(compress(vpred.get_low(), vpred.get_high())).store_a(dst + x);
    }
}

template<>
void interpolate<float, Vec8f, Vec8f, 8>(const float *cur0, const float *cur2, float * VS_RESTRICT dst, const int w)
{
    const Vec8f vpointFive(0.5f);
    
    for (int x = 0; x < w; x += 8)
    {
        const Vec8f vcur0 = Vec8f().load_a(cur0 + x);
        const Vec8f vcur2 = Vec8f().load_a(cur2 + x);
        
        const Vec8f vpred = (vcur0 + vcur2) * vpointFive;
        vpred.store_a(dst + x);
    }
}

#else

template<typename T1>
static void interpolate(const T1 *cur0, const T1 *cur2, T1 * VS_RESTRICT dst, const int w)
{
    for (int x = 0; x < w; x++)
    {
        dst[x] = (cur0[x] + cur2[x] + 1) >> 1;
    }
}

template<>
void interpolate<float>(const float *cur0, const float *cur2, float * VS_RESTRICT dst, const int w)
{
    for (int x = 0; x < w; x++)
    {
        dst[x] = (cur0[x] + cur2[x]) * 0.5f;
    }
}

#endif

template<typename T1>
static void interpolate_edges(const T1 *cur0, const T1 *cur2, T1 * VS_RESTRICT dst, const int start, const int end)
{
    for (int x = start; x < end; x++)
    {
        dst[x] = (cur0[x] + cur2[x] + 1) >> 1;
    }
}

template<>
void interpolate_edges<float>(const float *cur0, const float *cur2, float * VS_RESTRICT dst, const int start, const int end)
{
    for (int x = start; x < end; x++)
    {
        dst[x] = (cur0[x] + cur2[x]) * 0.5f;
    }
}

template<typename T1, typename T2, typename T3, typename V1, typename V2, typename V2b, int vectorSize>
static void filter_plane(const T1 *prev0, const T1 *cur0, const T1 *next0, T1 * VS_RESTRICT dst0,
                         const int w, const int h, const int stride,
                         const int tff, const int parity, const int mode)
{
#ifdef VS_TARGET_CPU_X86
    const int guardSize = 32;
    const int offset = guardSize / sizeof(T2);
    T2 *dp = vs_aligned_malloc<T2>(stride * sizeof(T2), 32);
    T2 *diffp = vs_aligned_malloc<T2>(stride * sizeof(T2), 32);
    T2 *dir0p = vs_aligned_malloc<T2>(stride * sizeof(T2) + guardSize, 32); // -2
    T2 *dir1p = vs_aligned_malloc<T2>(stride * sizeof(T2) + guardSize, 32);
    T2 *dir2p = vs_aligned_malloc<T2>(stride * sizeof(T2) + guardSize, 32);
    T2 *dir3p = vs_aligned_malloc<T2>(stride * sizeof(T2) + guardSize, 32);
    T2 *dir4p = vs_aligned_malloc<T2>(stride * sizeof(T2) + guardSize, 32);
    T3 *pred_dirp = vs_aligned_malloc<T3>(32, 32); // 32
#endif
    for (int y = parity; y < h; y += 2)
    {
        memcpy(dst0 + y * stride, cur0 + y * stride, w * sizeof(T1));
    }
    
    if (parity)
    {
        memcpy(dst0, cur0 + stride, w * sizeof(T1)); // 1 -> 0
    }
    else
    {
#ifdef VS_TARGET_CPU_X86
        interpolate<T1, V1, V2, vectorSize>(cur0, cur0 + 2 * stride, dst0 + stride, w);
#else
        interpolate<T1>(cur0, cur0 + 2 * stride, dst0 + stride, w); // interpolate 0 and 2
#endif
    }
    
    for (int y = 2 + (parity ^ 1); y < h - 2; y += 2)
    {
        const T1 *prev = prev0 + y * stride;
        const T1 *cur  = cur0  + y * stride;
        const T1 *next = next0 + y * stride;
        T1 * VS_RESTRICT dst = dst0 + y * stride;
        
#ifdef VS_TARGET_CPU_X86
        filter_line<T1, T2, T3, V1, V2, V2b, vectorSize>(prev, cur, next, dst, dp, diffp, dir0p + offset, dir1p + offset, dir2p + offset, dir3p + offset, dir4p + offset, pred_dirp, w, stride, (parity ^ tff), mode);
#else
        filter_line<T1>(prev, cur, next, dst, w, stride, (parity ^ tff), mode);
#endif
        interpolate_edges<T1>(cur - stride, cur + stride, dst, 0, 3);
        interpolate_edges<T1>(cur - stride, cur + stride, dst, (w - 3), w);
    }
    
    if ((h & 1) ^ parity)
    {
#ifdef VS_TARGET_CPU_X86
        interpolate<T1, V1, V2, vectorSize>(cur0 + (h - 3) * stride, cur0 + (h - 1) * stride, dst0 + (h - 2) * stride, w);
#else
        interpolate<T1>(cur0 + (h - 3) * stride, cur0 + (h - 1) * stride, dst0 + (h - 2) * stride, w); // interpolate h-3 and h-1
#endif
    }
    else
    {
        memcpy(dst0 + (h - 1) * stride, cur0 + (h - 2) * stride, w * sizeof(T1)); // h-2 -> h-1
    }
#ifdef VS_TARGET_CPU_X86
    vs_aligned_free(dp);
    vs_aligned_free(diffp);
    vs_aligned_free(dir0p);
    vs_aligned_free(dir1p);
    vs_aligned_free(dir2p);
    vs_aligned_free(dir3p);
    vs_aligned_free(dir4p);
    vs_aligned_free(pred_dirp);
#endif
}

// T1 - uint8_t, uint16_t, float
template<typename T1, typename T2, typename T3, typename V1, typename V2, typename V2b, int vectorSize>
static void filter_frame(const VSFrameRef *prev, const VSFrameRef *src, const VSFrameRef *next, VSFrameRef *dst,
                         const int tff, const int parity, const YadifData *d, const VSAPI *vsapi)
{
    for (int plane = 0; plane < d->vi.format->numPlanes; plane++)
    {
        const int width = vsapi->getFrameWidth(src, plane);
        const int height = vsapi->getFrameHeight(src, plane);
        const int stride = vsapi->getStride(src, plane) / sizeof(T1);
        const T1 *prevp = reinterpret_cast<const T1 *>(vsapi->getReadPtr(prev, plane));
        const T1 *srcp  = reinterpret_cast<const T1 *>(vsapi->getReadPtr(src, plane));
        const T1 *nextp = reinterpret_cast<const T1 *>(vsapi->getReadPtr(next, plane));
        T1 * VS_RESTRICT dstp = reinterpret_cast<T1 *>(vsapi->getWritePtr(dst, plane));

        filter_plane<T1, T2, T3, V1, V2, V2b, vectorSize>(prevp, srcp, nextp, dstp, width, height, stride, tff, parity, d->mode);
    }
}

static const VSFrameRef *VS_CC yadifGetFrame(int n, int activationReason, void **instanceData, void **frameData, VSFrameContext *frameCtx, VSCore *core, const VSAPI *vsapi)
{
    const YadifData *d = static_cast<const YadifData *>(*instanceData);

    if (activationReason == arInitial)
    {
        if (d->mode & 1)
            n /= 2;

        if (n > 0)
            vsapi->requestFrameFilter(n - 1, d->node, frameCtx);
        vsapi->requestFrameFilter(n, d->node, frameCtx);
        if (n < d->vi_src->numFrames - 1)
            vsapi->requestFrameFilter(n + 1, d->node, frameCtx);
    }
    else if (activationReason == arAllFramesReady)
    {
        const int ndst = n;
        if (d->mode & 1)
            n /= 2;

        const VSFrameRef *prev = vsapi->getFrameFilter((n > 0) ? (n - 1) : ((d->vi_src->numFrames > 1) ? 1 : 0), d->node, frameCtx);
        const VSFrameRef *src  = vsapi->getFrameFilter(n, d->node, frameCtx);
        const VSFrameRef *next = vsapi->getFrameFilter((n < d->vi_src->numFrames - 1) ? (n + 1) : ((d->vi_src->numFrames > 1) ? (d->vi_src->numFrames - 1) : 0), d->node, frameCtx);

        VSFrameRef *dst = vsapi->newVideoFrame(d->vi.format, d->vi.width, d->vi.height, src, core);
        
        int err;
        const int fieldBased = int64ToIntS(vsapi->propGetInt(vsapi->getFramePropsRO(src), "_FieldBased", 0, &err));
        const int tff = (fieldBased == 1) ? 0 : ((fieldBased == 2) ? 1 : d->order);

        const int parity = (d->mode & 1) ? ((ndst & 1) ^ (tff ^ 1)) : ((d->field == -1) ? (tff ^ 1) : (d->field ^ 1));

        if (d->vi.format->sampleType == stInteger)
        {
            if (d->vi.format->bitsPerSample == 8)
            {
#ifdef VS_TARGET_CPU_X86
                filter_frame<uint8_t, int16_t, int16_t, Vec16uc, Vec16s, Vec16sb, 16>(prev, src, next, dst, tff, parity, d, vsapi);
#else
                filter_frame<uint8_t, void, void, void, void, void, 0>(prev, src, next, dst, tff, parity, d, vsapi);
#endif
            }
            else
            {
#ifdef VS_TARGET_CPU_X86
                filter_frame<uint16_t, int32_t, int32_t, Vec8us, Vec8i, Vec8ib, 8>(prev, src, next, dst, tff, parity, d, vsapi);
#else
                filter_frame<uint16_t, void, void, void, void, void, 0>(prev, src, next, dst, tff, parity, d, vsapi);
#endif
            }
        }
        else
        {
#ifdef VS_TARGET_CPU_X86
            filter_frame<float, float, int32_t, Vec8f, Vec8f, Vec8fb, 8>(prev, src, next, dst, tff, parity, d, vsapi);
#else
            filter_frame<float, void, void, void, void, void, 8>(prev, src, next, dst, tff, parity, d, vsapi);
#endif
        }

        VSMap *props = vsapi->getFramePropsRW(dst);
        vsapi->propSetInt(props, "_FieldBased", 0, paReplace);

        if (d->mode & 1)
        {
            int errNum, errDen;
            int64_t durationNum = vsapi->propGetInt(props, "_DurationNum", 0, &errNum);
            int64_t durationDen = vsapi->propGetInt(props, "_DurationDen", 0, &errDen);
            if (!errNum && !errDen)
            {
                muldivRational(&durationNum, &durationDen, 1, 2);
                vsapi->propSetInt(props, "_DurationNum", durationNum, paReplace);
                vsapi->propSetInt(props, "_DurationDen", durationDen, paReplace);
            }
        }

        vsapi->freeFrame(prev);
        vsapi->freeFrame(src);
        vsapi->freeFrame(next);
        return dst;
    }

    return nullptr;
}

static void VS_CC yadifFree(void *instanceData, VSCore *core, const VSAPI *vsapi)
{
    YadifData *d = static_cast<YadifData *>(instanceData);
    vsapi->freeNode(d->node);
    delete d;
}

static void VS_CC yadifCreate(const VSMap *in, VSMap *out, void *userData, VSCore *core, const VSAPI *vsapi)
{
    YadifData d;
    int err;

    d.order = int64ToIntS(vsapi->propGetInt(in, "order", 0, nullptr));

    d.field = int64ToIntS(vsapi->propGetInt(in, "field", 0, &err));
    if (err)
        d.field = -1;

    d.mode = int64ToIntS(vsapi->propGetInt(in, "mode", 0, &err));

    if (d.order < 0 || d.order > 1)
    {
        vsapi->setError(out, "Yadif: order must be 0 or 1");
        return;
    }

    if (d.field < -1 || d.field > 1)
    {
        vsapi->setError(out, "Yadif: field must be -1, 0, or 1");
        return;
    }

    if (d.mode < 0 || d.mode > 3)
    {
        vsapi->setError(out, "Yadif: mode must be 0, 1, 2, or 3");
        return;
    }


    d.node = vsapi->propGetNode(in, "clip", 0, nullptr);
    d.vi = *vsapi->getVideoInfo(d.node);
    d.vi_src = vsapi->getVideoInfo(d.node);

    if (!isConstantFormat(&d.vi))
    {
        vsapi->setError(out, "Yadif: only constant format input supported");
        vsapi->freeNode(d.node);
        return;
    }

    if ((d.vi.format->sampleType == stInteger && d.vi.format->bitsPerSample > 16) ||
        (d.vi.format->sampleType == stFloat && d.vi.format->bitsPerSample != 32))
    {
        vsapi->setError(out, "Yadif: only 8-16 bit integer or 32 bit float input supported");
        vsapi->freeNode(d.node);
        return;
    }

    if (d.mode & 1)
    {
        if (d.vi.numFrames > INT_MAX / 2)
        {
            vsapi->setError(out, "Yadif: output clip would be too long");
            vsapi->freeNode(d.node);
            return;
        }

        d.vi.numFrames *= 2;

        if (d.vi.fpsNum && d.vi.fpsDen)
            muldivRational(&d.vi.fpsNum, &d.vi.fpsDen, 2, 1);
    }

    YadifData *data = new YadifData(d);
    
    vsapi->createFilter(in, out, "Yadif", yadifInit, yadifGetFrame, yadifFree, fmParallel, 0, data, core);
}


////////////////////////////////////////////

VS_EXTERNAL_API(void) VapourSynthPluginInit(VSConfigPlugin configFunc, VSRegisterFunction registerFunc, VSPlugin *plugin)
{
    configFunc("com.kusakinomori.yadif", "yadif", "Yadif for VapourSynth, r" YADIF_VERSION, VAPOURSYNTH_API_VERSION, 1, plugin);
    registerFunc("Yadif",
                 "clip:clip;"
                 "order:int;"
                 "field:int:opt;"
                 "mode:int:opt;",
                 yadifCreate, nullptr, plugin);
}

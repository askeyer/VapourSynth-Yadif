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

template<typename T>
static void filter_line(const T *prev, const T *cur, const T *next, T * VS_RESTRICT dst,
                        const int w, const int stride, const int parity, const int mode)
{
    const T *prev2 = parity ? prev : cur;
    const T *next2 = parity ? cur  : next;

    for (int x = 3; x < w - 3; x++)
    {
        const int c = cur[-stride];
        const int d = (prev2[0] + next2[0]) >> 1;
        const int e = cur[stride];
        const int tdiff0 = std::abs(prev2[0] - next2[0]) >> 1;
        const int tdiff1 = (std::abs(prev[-stride] - c) + std::abs(prev[stride] - e)) >> 1;
        const int tdiff2 = (std::abs(next[-stride] - c) + std::abs(next[stride] - e)) >> 1;
        int diff = std::max(std::max(tdiff0, tdiff1), tdiff2);

        // edi
        int spatial_pred = (c + e) >> 1;
        int spatial_score = std::abs(cur[-stride-1] - cur[stride-1]) + std::abs(c - e) + 
                            std::abs(cur[-stride+1] - cur[stride+1]) - 1;

        int score = std::abs(cur[-stride-2] - cur[stride]) +
                    std::abs(cur[-stride-1] - cur[stride+1]) +
                    std::abs(cur[-stride] - cur[stride+2]);
        if (score < spatial_score)
        {
            spatial_score = score;
            spatial_pred = (cur[-stride-1] + cur[stride+1]) >> 1;
            score = std::abs(cur[-stride-3] - cur[stride+1]) +
                    std::abs(cur[-stride-2] - cur[stride+2]) +
                    std::abs(cur[-stride-1] - cur[stride+3]);
            if (score < spatial_score)
            {
                spatial_score = score;
                spatial_pred = (cur[-stride-2] + cur[stride+2]) >> 1;
            }
        }

        score = std::abs(cur[-stride] - cur[stride-2]) +
                std::abs(cur[-stride+1] - cur[stride-1]) +
                std::abs(cur[-stride+2] - cur[stride]);
        if (score < spatial_score)
        {
            spatial_score = score;
            spatial_pred = (cur[-stride+1] + cur[stride-1]) >> 1;
            score = std::abs(cur[-stride+1] - cur[stride-3]) +
                    std::abs(cur[-stride+2] - cur[stride-2]) +
                    std::abs(cur[-stride+3] - cur[stride-1]);
            if (score < spatial_score)
            {
                spatial_score = score;
                spatial_pred = (cur[-stride+2] + cur[stride-2]) >> 1;
            }
        }

        // check
        if (mode < 2)
        {
            const int b = (prev2[-2*stride] + next2[-2*stride]) >> 1;
            const int f = (prev2[2*stride] + next2[2*stride]) >> 1;
            const int maxs = std::max(std::max(d - e, d - c), std::min(b - c, f - e));
            const int mins = std::min(std::min(d - e, d - c), std::max(b - c, f - e));
            diff = std::max(std::max(diff, mins), -maxs);
        }

        dst[0] = std::max(std::min<int>(spatial_pred, d + diff), d - diff);

        dst++;
        cur++;
        prev++;
        next++;
        prev2++;
        next2++;
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
        const float c = cur[-stride];
        const float d = (prev2[0] + next2[0]) * 0.5f;
        const float e = cur[stride];
        const float tdiff0 = std::abs(prev2[0] - next2[0]) * 0.5f;
        const float tdiff1 = (std::abs(prev[-stride] - c) + std::abs(prev[stride] - e)) * 0.5f;
        const float tdiff2 = (std::abs(next[-stride] - c) + std::abs(next[stride] - e)) * 0.5f;
        float diff = std::max(std::max(tdiff0, tdiff1), tdiff2);

        // edi
        float spatial_pred = (c + e) * 0.5f;
        float spatial_score = std::abs(cur[-stride-1] - cur[stride-1]) + std::abs(c - e) + 
                              std::abs(cur[-stride+1] - cur[stride+1]);

        float score = std::abs(cur[-stride-2] - cur[stride]) +
                      std::abs(cur[-stride-1] - cur[stride+1]) +
                      std::abs(cur[-stride] - cur[stride+2]);
        if (score < spatial_score)
        {
            spatial_score = score;
            spatial_pred = (cur[-stride-1] + cur[stride+1]) * 0.5f;
            score = std::abs(cur[-stride-3] - cur[stride+1]) +
                    std::abs(cur[-stride-2] - cur[stride+2]) +
                    std::abs(cur[-stride-1] - cur[stride+3]);
            if (score < spatial_score)
            {
                spatial_score = score;
                spatial_pred = (cur[-stride-2] + cur[stride+2]) * 0.5f;
            }
        }

        score = std::abs(cur[-stride] - cur[stride-2]) +
                std::abs(cur[-stride+1] - cur[stride-1]) +
                std::abs(cur[-stride+2] - cur[stride]);
        if (score < spatial_score)
        {
            spatial_score = score;
            spatial_pred = (cur[-stride+1] + cur[stride-1]) * 0.5f;
            score = std::abs(cur[-stride+1] - cur[stride-3]) +
                    std::abs(cur[-stride+2] - cur[stride-2]) +
                    std::abs(cur[-stride+3] - cur[stride-1]);
            if (score < spatial_score)
            {
                spatial_score = score;
                spatial_pred = (cur[-stride+2] + cur[stride-2]) * 0.5f;
            }
        }

        // check
        if (mode < 2)
        {
            const float b = (prev2[-2*stride] + next2[-2*stride]) * 0.5f;
            const float f = (prev2[2*stride] + next2[2*stride]) * 0.5f;
            const float maxs = std::max(std::max(d - e, d - c), std::min(b - c, f - e));
            const float mins = std::min(std::min(d - e, d - c), std::max(b - c, f - e));
            diff = std::max(std::max(diff, mins), -maxs);
        }

        dst[0] = std::max(std::min(spatial_pred, d + diff), d - diff);

        dst++;
        cur++;
        prev++;
        next++;
        prev2++;
        next2++;
    }
}

template<typename T>
static void interpolate(const T *cur0, const T *cur2, T * VS_RESTRICT dst, const int start, const int end)
{
    for (int x = start; x < end; x++)
    {
        dst[x] = (cur0[x] + cur2[x] + 1) >> 1;
    }
}

template<>
void interpolate<float>(const float *cur0, const float *cur2, float * VS_RESTRICT dst, const int start, const int end)
{
    for (int x = start; x < end; x++)
    {
        dst[x] = (cur0[x] + cur2[x]) * 0.5f;
    }
}

template<typename T>
static void filter_plane(const T *prev0, const T *cur0, const T *next0, T * VS_RESTRICT dst,
                         const int w, const int h, const int stride,
                         const int tff, const int parity, const int mode)
{
    int y = 0;
    if ((y ^ parity) & 1)
    {
        memcpy(dst, cur0 + stride, w * sizeof(T)); // 1 -> 0
    }
    else
    {
        memcpy(dst, cur0, w * sizeof(T));
    }

    y = 1;
    if ((y ^ parity) & 1)
    {
        interpolate<T>(cur0, cur0 + 2 * stride, dst + stride, 0, w); // interpolate 0 and 2
    }
    else
    {
        memcpy(dst + stride, cur0 + stride, w * sizeof(T));
    }

    for (y = 2; y < h - 2; y++)
    {
        if ((y ^ parity) & 1)
        {
            const T *prev = prev0 + y * stride + 3;
            const T *cur  = cur0  + y * stride + 3;
            const T *next = next0 + y * stride + 3;
            T * VS_RESTRICT dst2 = dst + y * stride + 3;
            
            interpolate<T>(cur0 + (y - 1) * stride, cur0 + (y + 1) * stride, dst + y * stride, 0, 3);
            filter_line<T>(prev, cur, next, dst2, w, stride, (parity ^ tff), mode);
            interpolate<T>(cur0 + (y - 1) * stride, cur0 + (y + 1) * stride, dst + y * stride, (w - 3), w);
        }
        else
        {
            memcpy(dst + y * stride, cur0 + y * stride, w * sizeof(T));
        }
    }

    y = h - 2;
    if ((y ^ parity) & 1)
    {
        interpolate<T>(cur0 + (y - 1) * stride, cur0 + (y + 1) * stride, dst + y * stride, 0, w); // interpolate h-3 and h-1
    }
    else
    {
        memcpy(dst + y * stride, cur0 + y * stride, w * sizeof(T));
    }

    y = h - 1;
    if ((y ^ parity) & 1)
    {
        memcpy(dst + y * stride, cur0 + (y - 1) * stride, w * sizeof(T)); // h-2 -> h-1
    }
    else
    {
        memcpy(dst + y * stride, cur0 + y * stride, w * sizeof(T));
    }
}

// uint8_t, uint16_t, float
template<typename T>
static void filter_frame(const VSFrameRef *prev, const VSFrameRef *src, const VSFrameRef *next, VSFrameRef *dst,
                         const int tff, const int parity, const YadifData *d, const VSAPI *vsapi)
{
    for (int plane = 0; plane < d->vi.format->numPlanes; plane++)
    {
        const int width = vsapi->getFrameWidth(src, plane);
        const int height = vsapi->getFrameHeight(src, plane);
        const int stride = vsapi->getStride(src, plane) / sizeof(T);
        const T *prevp = reinterpret_cast<const T *>(vsapi->getReadPtr(prev, plane));
        const T *srcp  = reinterpret_cast<const T *>(vsapi->getReadPtr(src, plane));
        const T *nextp = reinterpret_cast<const T *>(vsapi->getReadPtr(next, plane));
        T * VS_RESTRICT dstp = reinterpret_cast<T *>(vsapi->getWritePtr(dst, plane));

        filter_plane<T>(prevp, srcp, nextp, dstp, width, height, stride, tff, parity, d->mode);
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
                filter_frame<uint8_t>(prev, src, next, dst, tff, parity, d, vsapi);
            }
            else
            {
                filter_frame<uint16_t>(prev, src, next, dst, tff, parity, d, vsapi);
            }
        }
        else
        {
            filter_frame<float>(prev, src, next, dst, tff, parity, d, vsapi);
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

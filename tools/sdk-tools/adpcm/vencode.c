#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vadpcm.h"


/**
 * Number of bits per sample for a frame size: 9 bytes -> 4 bits, 7 -> 3, 5 -> 2.
 */
s32 vadpcm_bits_for_frame_size(s32 frameSize)
{
    return frameSize == 9 ? 4 : frameSize == 7 ? 3 : 2;
}

/**
 * Write the 16 quantized samples of a frame (each in [-2^(bits-1), 2^(bits-1)-1])
 * after the header byte, in the layout the RSP decoders expect:
 *   4-bit (8 bytes): two samples per byte, first sample in the high nibble
 *   2-bit (4 bytes): four samples per byte, first sample in the top two bits
 *   3-bit (6 bytes): the upper two bits of each sample packed like a 2-bit frame,
 *                    then one byte per half holding the low bits (bit n = sample n)
 */
void vadpcm_pack_frame(FILE *ofile, s32 *ix, s32 frameSize)
{
    s32 i;
    u8 c;

    if (frameSize == 9)
    {
        for (i = 0; i < 16; i += 2)
        {
            c = (ix[i] << 4) | (ix[i + 1] & 0xf);
            fwrite(&c, 1, 1, ofile);
        }
        return;
    }
    for (i = 0; i < 16; i += 4)
    {
        s32 sh = (frameSize == 7) ? 1 : 0;
        c = (((ix[i] >> sh) & 3) << 6) | (((ix[i + 1] >> sh) & 3) << 4) | (((ix[i + 2] >> sh) & 3) << 2) | ((ix[i + 3] >> sh) & 3);
        fwrite(&c, 1, 1, ofile);
    }
    if (frameSize == 7)
    {
        for (i = 0; i < 16; i += 8)
        {
            s32 j;
            c = 0;
            for (j = 0; j < 8; j++)
            {
                c |= (ix[i + j] & 1) << j;
            }
            fwrite(&c, 1, 1, ofile);
        }
    }
}

/**
 * Closed-loop encoder used for 2-bit and 3-bit frames: try every predictor and
 * every scale, decode each candidate exactly like the RSP does (prediction with
 * floor(x / 2048), outputs clamped to 16 bits) and keep the candidate with
 * the smallest reconstruction error. With only 4 or 8 levels per sample the
 * SDK's open-loop heuristics (predictor from the unquantized error, scale from
 * the peak error) clip and mistrack badly, and the search space is tiny.
 */
static s32 clamp16(s32 x)
{
    return x < -0x8000 ? -0x8000 : x > 0x7fff ? 0x7fff : x;
}

static void vencodeframe_closedloop(FILE *ofile, s16 *inBuffer, s32 *state, s32 ***coefTable, s32 order, s32 npredictors, s32 nsam, s32 frameSize)
{
    s32 encBits = vadpcm_bits_for_frame_size(frameSize);
    s32 maxScale = 16 - encBits;
    s32 llevel = -(1 << (encBits - 1));
    s32 ulevel = -llevel - 1;
    s32 inVector[16];
    s32 out[16];
    s32 ix[16];
    s32 bestIx[16];
    s32 bestOut[16];
    s32 bestP = 0;
    s32 bestScale = 0;
    s64 bestErr = -1;
    s32 p;
    s32 scale;
    s32 half;
    s32 i;
    u8 header;

    for (i = nsam; i < 16; i++)
    {
        inBuffer[i] = 0;
    }

    for (p = 0; p < npredictors; p++)
    {
        for (scale = 0; scale <= maxScale; scale++)
        {
            s64 err = 0;

            for (half = 0; half < 2 && (bestErr < 0 || err < bestErr); half++)
            {
                for (i = 0; i < order; i++)
                {
                    inVector[i] = (half == 0) ? state[16 - order + i] : out[8 - order + i];
                }
                for (i = 0; i < 8; i++)
                {
                    s32 n = half * 8 + i;
                    s32 prediction = inner_product(order + i, coefTable[p][i], inVector);
                    s32 q = qsample((f32) inBuffer[n] - (f32) prediction, 1 << scale);
                    s32 d;

                    q = clip(q, llevel, ulevel);
                    ix[n] = q;
                    inVector[i + order] = q * (1 << scale);
                    out[n] = clamp16(prediction + inVector[i + order]);
                    d = inBuffer[n] - out[n];
                    err += (s64) d * d;
                }
            }

            if (half == 2 && (bestErr < 0 || err < bestErr))
            {
                bestErr = err;
                bestP = p;
                bestScale = scale;
                for (i = 0; i < 16; i++)
                {
                    bestIx[i] = ix[i];
                    bestOut[i] = out[i];
                }
            }
        }
    }

    for (i = 0; i < 16; i++)
    {
        state[i] = bestOut[i];
    }

    header = (bestScale << 4) | (bestP & 0xf);
    fwrite(&header, 1, 1, ofile);
    vadpcm_pack_frame(ofile, bestIx, frameSize);
}

/**
 * Encode one frame of up to 16 samples. frameSize selects the format:
 *   9: standard VADPCM, 4-bit samples (1 header byte + 8 data bytes)
 *   7: 3-bit samples (1 header byte + 6 data bytes)
 *   5: "small" VADPCM, 2-bit samples (1 header byte + 4 data bytes)
 */
void vencodeframe(FILE *ofile, s16 *inBuffer, s32 *state, s32 ***coefTable, s32 order, s32 npredictors, s32 nsam, s32 frameSize)
{
    if (frameSize != 9)
    {
        vencodeframe_closedloop(ofile, inBuffer, state, coefTable, order, npredictors, nsam, frameSize);
        return;
    }
    {
    s16 ix[16];
    s32 prediction[16];
    s32 inVector[16];
    s32 saveState[16];
    s32 optimalp;
    s32 scale;
    s32 llevel;
    s32 ulevel;
    s32 i;
    s32 j;
    s32 k;
    s32 ie[16];
    s32 nIter;
    s32 max;
    s32 cV;
    s32 maxClip;
    u8 header;
    s32 encBits = vadpcm_bits_for_frame_size(frameSize);
    s32 maxScale = 16 - encBits; // 12 for 4-bit samples, 13 for 3-bit, 14 for 2-bit
    f32 e[16];
    f32 se;
    f32 min;

    // We are only given 'nsam' samples; pad with zeroes to 16.
    for (i = nsam; i < 16; i++)
    {
        inBuffer[i] = 0;
    }

    llevel = -(1 << (encBits - 1)); // [-8, 7] for 4-bit samples, [-4, 3] for 3-bit, [-2, 1] for 2-bit
    ulevel = -llevel - 1;

    // Determine the best-fitting predictor.
    min = 1e30;
    optimalp = 0;
    for (k = 0; k < npredictors; k++)
    {
        // Copy over the last 'order' samples from the previous output.
        for (i = 0; i < order; i++)
        {
            inVector[i] = state[16 - order + i];
        }

        // For 8 samples...
        for (i = 0; i < 8; i++)
        {
            // Compute a prediction based on 'order' values from the old state,
            // plus previous errors in this chunk, as an inner product with the
            // coefficient table.
            prediction[i] = inner_product(order + i, coefTable[k][i], inVector);
            // Record the error in inVector (thus, its first 'order' samples
            // will contain actual values, the rest will be error terms), and
            // in floating point form in e (for no particularly good reason).
            inVector[i + order] = inBuffer[i] - prediction[i];
            e[i] = (f32) inVector[i + order];
        }

        // For the next 8 samples, start with 'order' values from the end of
        // the previous 8-sample chunk of inBuffer. (The code is equivalent to
        // inVector[i] = inBuffer[8 - order + i].)
        for (i = 0; i < order; i++)
        {
            inVector[i] = prediction[8 - order + i] + inVector[8 + i];
        }

        // ... and do the same thing as before to get predictions.
        for (i = 0; i < 8; i++)
        {
            prediction[8 + i] = inner_product(order + i, coefTable[k][i], inVector);
            inVector[i + order] = inBuffer[8 + i] - prediction[8 + i];
            e[8 + i] = (f32) inVector[i + order];
        }

        // Compute the L2 norm of the errors; the lowest norm decides which
        // predictor to use.
        se = 0.0f;
        for (j = 0; j < 16; j++)
        {
            se += e[j] * e[j];
        }

        if (se < min)
        {
            min = se;
            optimalp = k;
        }
    }

    // Do exactly the same thing again, for real.
    for (i = 0; i < order; i++)
    {
        inVector[i] = state[16 - order + i];
    }

    for (i = 0; i < 8; i++)
    {
        prediction[i] = inner_product(order + i, coefTable[optimalp][i], inVector);
        inVector[i + order] = inBuffer[i] - prediction[i];
        e[i] = (f32) inVector[i + order];
    }

    for (i = 0; i < order; i++)
    {
        inVector[i] = prediction[8 - order + i] + inVector[8 + i];
    }

    for (i = 0; i < 8; i++)
    {
        prediction[8 + i] = inner_product(order + i, coefTable[optimalp][i], inVector);
        inVector[i + order] = inBuffer[8 + i] - prediction[8 + i];
        e[8 + i] = (f32) inVector[i + order];
    }

    // Clamp the errors to 16-bit signed ints, and put them in ie.
    clamp(16, e, ie, 16);

    // Find a value with highest absolute value.
    // @bug If this first finds -2^n and later 2^n, it should set 'max' to the
    // latter, which needs a higher value for 'scale'.
    max = 0;
    for (i = 0; i < 16; i++)
    {
        if (fabs(ie[i]) > fabs(max))
        {
            max = ie[i];
        }
    }

    // Compute which power of two we need to scale down by in order to make
    // all values representable as 'encBits'-bit signed integers (i.e. be in
    // [-8, 7] or [-2, 1]). The worst-case 'max' is -2^15, so this will be at
    // most maxScale (12 or 14).
    for (scale = 0; scale <= maxScale; scale++)
    {
        if (max <= ulevel && max >= llevel)
        {
            goto out;
        }
        max /= 2;
    }
out:;

    for (i = 0; i < 16; i++)
    {
        saveState[i] = state[i];
    }

    // Try with the computed scale, but if it turns out we don't fit in encBits
    // bits (if some |cV| >= 2), use scale + 1 instead (i.e. downscaling by
    // another factor of 2).
    scale--;
    nIter = 0;
    do
    {
        nIter++;
        maxClip = 0;
        scale++;
        if (scale > maxScale)
        {
            scale = maxScale;
        }

        // Copy over the last 'order' samples from the previous output.
        for (i = 0; i < order; i++)
        {
            inVector[i] = saveState[16 - order + i];
        }

        // For 8 samples...
        for (i = 0; i < 8; i++)
        {
            // Compute a prediction based on 'order' values from the old state,
            // plus previous *quantized* errors in this chunk (because that's
            // all the decoder will have available).
            prediction[i] = inner_product(order + i, coefTable[optimalp][i], inVector);

            // Compute the error, and divide it by 2^scale, rounding to the
            // nearest integer. This should ideally result in an encBits-bit integer.
            se = (f32) inBuffer[i] - (f32) prediction[i];
            ix[i] = qsample(se, 1 << scale);

            // Clamp the error to an encBits-bit signed integer, and record what
            // delta was needed for that.
            cV = (s16) clip(ix[i], llevel, ulevel) - ix[i];
            if (maxClip < abs(cV))
            {
                maxClip = abs(cV);
            }
            ix[i] += cV;

            // Record the quantized error in inVector for later predictions,
            // and the quantized (decoded) output in state (for use in the next
            // batch of 8 samples).
            inVector[i + order] = ix[i] * (1 << scale);
            state[i] = prediction[i] + inVector[i + order];
        }

        // Copy over the last 'order' decoded samples from the above chunk.
        for (i = 0; i < order; i++)
        {
            inVector[i] = state[8 - order + i];
        }

        // ... and do the same thing as before.
        for (i = 0; i < 8; i++)
        {
            prediction[8 + i] = inner_product(order + i, coefTable[optimalp][i], inVector);
            se = (f32) inBuffer[8 + i] - (f32) prediction[8 + i];
            ix[8 + i] = qsample(se, 1 << scale);
            cV = (s16) clip(ix[8 + i], llevel, ulevel) - ix[8 + i];
            if (maxClip < abs(cV))
            {
                maxClip = abs(cV);
            }
            ix[8 + i] += cV;
            inVector[i + order] = ix[8 + i] * (1 << scale);
            state[8 + i] = prediction[8 + i] + inVector[i + order];
        }
    }
    while (maxClip >= 2 && nIter < 2);

    // The scale and the predictor index are 4-bit numbers, and the 16 computed
    // outputs are encBits-bit numbers. Write them out as 1 + 8 bytes (4-bit
    // samples), 1 + 6 bytes (3-bit) or 1 + 4 bytes (2-bit).
    header = (scale << 4) | (optimalp & 0xf);
    fwrite(&header, 1, 1, ofile);
    {
        s32 ix32[16];
        for (i = 0; i < 16; i++)
        {
            ix32[i] = ix[i];
        }
        vadpcm_pack_frame(ofile, ix32, frameSize);
    }
    }
}

#include <stdio.h>
#include "vadpcm.h"

/**
 * Decode one 16-sample frame. frameSize is 9 for 4-bit samples, 7 for 3-bit
 * samples and 5 for 2-bit ("small" VADPCM) samples. See vadpcm_pack_frame in
 * vencode.c for the layouts.
 */
void vdecodeframe(FILE *ifile, s32 *outp, s32 order, s32 ***coefTable, s32 frameSize)
{
    s32 optimalp;
    s32 scale;
    s32 maxlevel;
    s32 i;
    s32 j;
    s32 in_vec[16];
    s32 ix[16];
    u8 header;
    u8 c;

    fread(&header, 1, 1, ifile);
    scale = 1 << (header >> 4);
    optimalp = header & 0xf;

    if (frameSize == 5 || frameSize == 7)
    {
        s32 bits = (frameSize == 7) ? 3 : 2;

        maxlevel = (1 << (bits - 1)) - 1;
        for (i = 0; i < 16; i += 4)
        {
            fread(&c, 1, 1, ifile);
            ix[i] = c >> 6;
            ix[i + 1] = (c >> 4) & 3;
            ix[i + 2] = (c >> 2) & 3;
            ix[i + 3] = c & 3;
        }
        if (frameSize == 7)
        {
            for (i = 0; i < 16; i += 8)
            {
                fread(&c, 1, 1, ifile);
                for (j = 0; j < 8; j++)
                {
                    ix[i + j] = (ix[i + j] << 1) | ((c >> j) & 1);
                }
            }
        }
        for (i = 0; i < 16; i++)
        {
            if (ix[i] <= maxlevel)
            {
                ix[i] *= scale;
            }
            else
            {
                ix[i] = (-(1 << bits) - -ix[i]) * scale;
            }
        }
    }
    else
    {
        maxlevel = 7;
        for (i = 0; i < 16; i += 2)
        {
            fread(&c, 1, 1, ifile);
            ix[i] = c >> 4;
            ix[i + 1] = c & 0xf;

            if (ix[i] <= maxlevel)
            {
                ix[i] *= scale;
            }
            else
            {
                ix[i] = (-0x10 - -ix[i]) * scale;
            }

            if (ix[i + 1] <= maxlevel)
            {
                ix[i + 1] *= scale;
            }
            else
            {
                ix[i + 1] = (-0x10 - -ix[i + 1]) * scale;
            }
        }
    }

    for (j = 0; j < 2; j++)
    {
        for (i = 0; i < 8; i++)
        {
            in_vec[i + order] = ix[j * 8 + i];
        }

        if (j == 0)
        {
            for (i = 0; i < order; i++)
            {
                in_vec[i] = outp[16 - order + i];
            }
        }
        else
        {
            for (i = 0; i < order; i++)
            {
                in_vec[i] = outp[j * 8 - order + i];
            }
        }

        for (i = 0; i < 8; i++)
        {
            outp[i + j * 8] = inner_product(order + 8, coefTable[optimalp][i], in_vec);
        }
    }
}

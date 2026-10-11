#include "macros.inc"

.section .data

.balign 16
glabel gSoundDataADSR
.incbin "sound/sound_data.ctl"

.balign 16
glabel gSoundDataRaw
.incbin "sound/sound_data.tbl"

.balign 16
glabel gMusicData
.incbin "sound/sequences.bin"

#if !defined(VERSION_SH) && !defined(VERSION_CN)
.balign 16
glabel gBankSetsData
.incbin "sound/bank_sets"
#endif

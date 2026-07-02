#!/bin/bash
# Phase 12 CREST ensembles -> confirm exposed polar SASA in 3D (QM tier).
# Reserve GFN2 for finalists; GFN-FF screening tier otherwise.

crest p12_01_notail_polaraxis_ax8_seryl.xyz --gfnff --alpb water --quick -ewin 6 --T 8 > p12_01.out # reward=0.4763
crest p12_02_notail_polaraxis_ax9_polyhydroxy.xyz --gfnff --alpb water --quick -ewin 6 --T 8 > p12_02.out # reward=0.4715
crest p12_03_notail_polaraxis_ax7_succinamoyl.xyz --gfnff --alpb water --quick -ewin 6 --T 8 > p12_03.out # reward=0.4628
crest p12_04_notail_polaraxis_ax4_glycolyl.xyz --gfnff --alpb water --quick -ewin 6 --T 8 > p12_04.out # reward=0.459
crest p12_05_notail_polaraxis_ax6_carboxyethanoyl.xyz --gfnff --alpb water --quick -ewin 6 --T 8 > p12_05.out # reward=0.4525
crest p12_06_notail_polaraxis_ax11_sulfobenzoyl.xyz --gfnff --alpb water --quick -ewin 6 --T 8 > p12_06.out # reward=0.4454
crest p12_07_notail_polaraxis_ax5_diglycolyl.xyz --gfnff --alpb water --quick -ewin 6 --T 8 > p12_07.out # reward=0.4328
crest p12_08_notail_polaraxis_ax3_4OH_benzoyl.xyz --gfnff --alpb water --quick -ewin 6 --T 8 > p12_08.out # reward=0.4252
crest p12_09_notail_polaraxis_ax10_aminomethylbenzoyl.xyz --gfnff --alpb water --quick -ewin 6 --T 8 > p12_09.out # reward=0.4163
crest p12_10_notail_design_naphthoyl_6OH.xyz --gfnff --alpb water --quick -ewin 6 --T 8 > p12_10.out # reward=0.4042
crest p12_11_notail_polaraxis_ax2_benzoyl.xyz --gfnff --alpb water --quick -ewin 6 --T 8 > p12_11.out # reward=0.402
crest p12_12_notail_design_quinolinecarbonyl.xyz --gfnff --alpb water --quick -ewin 6 --T 8 > p12_12.out # reward=0.3915

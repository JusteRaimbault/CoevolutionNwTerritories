
# figure 1
montage figuresraw/ex_synthvirtual_1_t0.png figuresraw/ex_synthvirtual_1_t30.png figuresraw/ex_synthphysical_1_t0.png figuresraw/ex_synthphysical_1_t30.png figuresraw/ex_realphysical_1975_t0.png figuresraw/ex_realphysical_1975_tf.png -resize 500x -tile 2x3 -geometry +5+5 figures/Fig1.png

# figure 2
montage figuresraw/hierarchiesPopAlpha_nwExp1_wG0_001_xgravityDecay_colgravityGamma_facetsynthRankSize-nwThresholdQuantile.png figuresraw/hierarchiesClosenessAlpha_nwExp1_wG0_001_xgravityDecay_colgravityGamma_facetsynthRankSize-nwThresholdQuantile.png figuresraw/rankCorrsPopCloseness_nwExp1_wG0_001_xgravityDecay_colgravityGamma_facetsynthRankSize-nwThresholdQuantile.png figuresraw/segHierarchiesClosenessPsi_nwExp1_wG0_001_xgravityDecay_colgravityGamma_facetsynthRankSize-nwThresholdQuantile.png  -resize 500x -tile 2x2 -geometry +5+5 figures/Fig2.png

# figure 2 pdf
montage figuresraw/hierarchiesPopAlpha_nwExp1_wG0_001_xgravityDecay_colgravityGamma_facetsynthRankSize-nwThresholdQuantile.pdf figuresraw/hierarchiesClosenessAlpha_nwExp1_wG0_001_xgravityDecay_colgravityGamma_facetsynthRankSize-nwThresholdQuantile.pdf figuresraw/rankCorrsPopCloseness_nwExp1_wG0_001_xgravityDecay_colgravityGamma_facetsynthRankSize-nwThresholdQuantile.pdf figuresraw/segHierarchiesClosenessPsi_nwExp1_wG0_001_xgravityDecay_colgravityGamma_facetsynthRankSize-nwThresholdQuantile.pdf  -resize 500x -tile 2x2 -geometry +5+5 figures/Fig2.pdf


# figure 3
montage figuresraw/physical_hierarchiesPopAlpha_nwExp1_wG0_001_xgravityDecay_colgravityGamma_facetsynthRankSize-nwThresholdQuantile.png figuresraw/physical_hierarchiesClosenessAlpha_nwExp1_wG0_001_xgravityDecay_colgravityGamma_facetsynthRankSize-nwThresholdQuantile.png figuresraw/physical_rankCorrsPopCloseness_nwExp1_wG0_001_xgravityDecay_colgravityGamma_facetsynthRankSize-nwThresholdQuantile.png  figuresraw/physical_segHierarchiesClosenessPsi_nwExp1_wG0_001_xgravityDecay_colgravityGamma_facetsynthRankSize-nwThresholdQuantile.png -resize 500x -tile 2x2 -geometry +5+5 figures/Fig3.png

# figure 3 pdf
montage figuresraw/physical_hierarchiesPopAlpha_nwExp1_wG0_001_xgravityDecay_colgravityGamma_facetsynthRankSize-nwThresholdQuantile.pdf figuresraw/physical_hierarchiesClosenessAlpha_nwExp1_wG0_001_xgravityDecay_colgravityGamma_facetsynthRankSize-nwThresholdQuantile.pdf figuresraw/physical_rankCorrsPopCloseness_nwExp1_wG0_001_xgravityDecay_colgravityGamma_facetsynthRankSize-nwThresholdQuantile.pdf  figuresraw/physical_segHierarchiesClosenessPsi_nwExp1_wG0_001_xgravityDecay_colgravityGamma_facetsynthRankSize-nwThresholdQuantile.pdf -tile 2x2 figures/Fig3.pdf


# figure 4
cp figuresraw/scatterdeltaobjs_colorgravityDecay_samples10.png figures/Fig4.png



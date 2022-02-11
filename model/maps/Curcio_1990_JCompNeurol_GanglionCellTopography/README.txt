Data from this paper:

Curcio, Christine A., and Kimberly A. Allen. "Topography of ganglion cells in human retina." Journal of comparative Neurology 300.1 (1990): 5-25.

The file 'curcioRGCDensityPerSqMm_average.mat' contains the mean RGC density values (and associated support in deg eccentricity) reported in this paper. The values themselves were derived from the file 'Curcio_JCompNeurol1990_GCtop_FS.xls', which was downloaded from https://info.cis.uab.edu/curcio/GanglionCellTopography/ on 6/28/2017.

Additional files contain the density measurements from each eye studied in this paper. The file 'DENSITY5_gc.xlsx' was provided by Christine Curcio by email to gkaguirre@me.com on 10/24/2017. The values in this file were then transposed and arranged and saved as 'DENSITY5_gc_resorted.csv'. The routine 'makeIndividualRGCDensityFiles.m' was used to read the rearranged individual subject data file and then save separate .m files with the measurements from each subject in a standard format.

Curcio's measurements were made on the retinal surface in units counts per mm^2. In discussions with Christine, we determined that the original units of measurement of position were in degrees of retinal angle in the eyeball. Specifically, globe diameter was measured with calipers and the site of density measurement was expressed as distance from the fovea in degrees mapped to a 360° sphere of 360°.

We confirmed that the average values reported in the main paper are reconstructed by first assigning to subject 29986 the average values from that subject's left and right eye, and then taking the average of all subject values.
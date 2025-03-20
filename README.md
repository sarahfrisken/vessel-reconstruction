# vessel-reconstruction
A C++ visual studio application for reconstructing a 3D vesselness volume from 2D DSA images. It is assumed the DSA images have been segmented using the VESCL contouring library to generate annotated 2D images.

The value of each voxel in a vesselness volume represents the probability that that voxel is inside a vessel. A probability iso-surface could be used to build a model of the vascular structure. However, more processing of the vesselness volume to find vessel struture will likely provide better models. We intend to investigate such processing methods in the future.

Please cite the following paper if you make use of this software for any purpose. S. Frisken et al. "Spatiotemporally Constrained 3D Reconstruction from Biplanar Digital Subtraction Angiography", Journal of Computer Assisted Radiology and Surgery, 2025.

This directory enables the complete implementation of “template CoMFA”, as introduced in the following references.

1. Cramer RD. Template CoMFA Generates Single 3D-QSAR Models that, for Twelve of Twelve Biological Targets, Predict All ChEMBL-Tabulated Affinities. PloS one 10.6 (2015): e0129307. doi.org/10.1371/journal.pone.0129307
2. Cramer RD, Wendt B. Template CoMFA: The 3D-QSAR Grail? J Chem Inf Model. 2014; 54:660–671 doi: 10.1021/ci400696v. pmid:24437630 
3. Cramer RD. Template CoMFA applied to 114 Biological Targets. J Chem Inf Model. 2014; 57:2147–2156 
4. Wendt B, Cramer RD. Challenging the Gold Standard for 3D-QSAR: Template CoMFA versus X-ray Alignment. J. Comp-Aided Mol Des. 2014; 28:803–824 doi: 10.1007/s10822-014-9761-z. pmid:24934658

The “tcfa” program performs the template CoMFA alignments themselves. Although its alignment algorithm significantly differs from the one described in the references, its output produces very similar models.

The “comfa” program performs the same default workflow that Sybyl formerly provided, including Gasteiger-Huckel charge calculation, steric and electrostatic field sampling, SAMPLS-assisted model generation, prediction, and visualization support.

Both programs are command-line utilities, thus deployable in any environment, and are provided in C++ source code form, to be built as described in Open Eye’s documentation. Access to Open Eye’s C++ API is of course required, for both building and running these programs. Both codes have been built in OSX and Linux (Ubuntu) environments, but not tried in Windows. The contents of the “data” directory are provided for validating these builds.

Please direct any feedback or questions to rdcramer3cfa@gmail.com. 

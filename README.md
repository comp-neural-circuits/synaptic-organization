Emergence of local and global synaptic organization on cortical dendrites
====================================================================

Matlab code for the manuscript "Emergence of local and global synaptic organization on cortical dendrites" by Jan H. Kirchner and Julijana Gjorgjieva

Getting Started
---------------

Download archive and extract into a folder on your system. When using a Windows machine, try to avoid path names with spaces or special characters.

### Prerequisites

We wrote and executed the code with Matlab 2019a, but previous versions of Matlab should also be compatible. Please download the toolboxes and functions listed under "Acknowledgments" and save them in the subfolder "tools", otherwise the code will require substantial rewriting to run properly.

### Reproducing the figures

To the matlab figure that was used to generate a panel of a figure, open the appropriately named script in the corresponding folder in matlab and execute.

Example:

To generate panel D from Figure 4, open the script PanelD.m in the subfolder Fig4 in matlab and press F5 (or click "Run").

### Reproducing the simulations

Some panels combine output from several simulations. To generate these panels, either use the output we provide (subfolder "sims") or run the corresponding batch_process_\[\].m script to run and store the required simulations (potentially very long runtimes depending on system).

Example:

To generate panel G from Figure 4, open the script PanelG.m in the subfolder Fig4 in matlab and press F5 (or click "Run"). If the message "run batch_process_ferret.m first" is displayed, open batch_process_ferret.m and press F5 (or click "Run"). For this panel, expect a runtime of ~20 minutes per simulation on a normal desktop machine, i.e. 100 hours for 300 simulations in total.

License
-------

MIT License

Copyright (c) \[2019\] \[Jan H Kirchner\]

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

Acknowledgments
---------------

*   We thank the authors of
    *   the TREES toolbox (https://www.treestoolbox.org/),
    *   the corPop toolbox (https://de.mathworks.com/matlabcentral/fileexchange/20591-sampling-from-multivariate-correlated-binary-and-poisson-random-variables?focused=5103226&tab=function)
    *   the circ\_stat toolbox (http://bethgelab.org/software/circstat/)
    *   the cbrewer function (https://de.mathworks.com/matlabcentral/fileexchange/34087-cbrewer-colorbrewer-schemes-for-matlab)
    *   the gramm toolbox (https://de.mathworks.com/matlabcentral/fileexchange/54465-gramm-complete-data-visualization-toolbox-ggplot2-r-like)
    *   the alphamask function (https://de.mathworks.com/matlabcentral/fileexchange/34936-alphamask-semi-transparent-image-overlay)
    *   the rdir function (https://de.mathworks.com/matlabcentral/fileexchange/47125-rdir-m)
    *   the rgb function (https://de.mathworks.com/matlabcentral/fileexchange/24497-rgb-triple-of-color-name-version-2?focused=5124709&tab=function)
    *   the scatstat1 function (https://de.mathworks.com/matlabcentral/fileexchange/57945-scatstat1-local-statistics-of-scattered-1d-data)
    *   the trandn function (https://de.mathworks.com/matlabcentral/fileexchange/53180-truncated-normal-generator)

[Access input movies and data points extracted from publications](https://www.dropbox.com/sh/rbtum08gmmbufl3/AACQ-mHfJYoEN6g8wi1xdukKa?dl=0)

[Access to precomputed simulation output](https://www.dropbox.com/sh/0jjjrfit5qs5zen/AACLGdXuP8J6d8Ki5u16CI0ba?dl=0)

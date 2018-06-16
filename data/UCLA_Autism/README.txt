This is fMRI and DTI-based data downloaded from http://umcd.humanconnectomeproject.org/,
from the UCLA Autism study. The study is described in Rudie et al., 
Altered functional and structural brain network organization in autism,
NeuroImage Clinical, 2013

The graphs are in the 264-region atlas from Power et al, Neuron, 2011


1) Resting-state fMRI dataset

- Subjects: 42 Autism Spectrum Disorder (36 male)
+ 37 Typically Developing children (31 male) from 8-17 y.o.

- Acquisition: Siemens 3T, 6 minutes resting state, TR=3s, TE=28ms, 3x3x4mm voxels	

- Preprocessing:
Motion correction, brain extraction, bandpass filtering (.1-.01Hz)
regression: motion parameters, motion spikes, WM, CSF, global signals
Extraction of timeseries from regions in transformed to subject native space

- Findings: 
ASD -> reduced short and long-range connectivity within functional systems,
stronger connectivity between functional systems,  particularly in default
and higher-order visual regions. Network level reductions in modularity and
clustering (local efficiency), but shorter characteristic path lengths
(higher global efficiency). 


2) DTI dataset

- Subjects: 51 Autism Spectrum Disorder (45 male)
+ 43 Typically Developing (36 male)

- Acquisition: Siemens 3T

- Preprocessing: 
brain extraction, eddy correction/motion correction
deterministic tractography (5mm fibers excluded), Edges are numbers of fibers
between regions in subject native space 

- Findings: lower levels of white matter integrity yet higher numbers of fibers

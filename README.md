# BE-562-Toehold-Project

### Project Description
Synthetic biologists are increasingly looking to develop models that predict the functional performance of toehold switches, which are RNA-based regulatory devices that can both bind to cognate mRNAs and gate protein translation. Recent efforts to develop a deep neural network to predict toehold switch performance led to a publicly available dataset of toehold variants, consisting of over 100,000 RNA sequences each with annotated variable and constant design regions as well as a corresponding level of functional performance. Given this resource, we developed three separate models with distinct architectures in an attempt to discern any latent features found among high- or low-performing toehold switches, as well as predict the performance of toehold switches as a function of their sequence. We identified enriched instances of A and T nucleotides in the hairpin region of high functioning toehold sequences by implementing a classic one occurrence per sequence (OOPS) model, and then developed a preliminary classifier to distinguish between sequences as ‘high’ or ‘low’ performing given its sequence. We then implemented a biologically-inspired neural network architecture to model the unfolding phenomena associated with leaky, poorly-performing toehold switches; this enabled us to generate OFF predictions that rival those published in the literature. Our findings demonstrate that the biological principles of thermodynamics, kinetics, and molecular diffusion (e.g., spontaneous second structure unfolding governed by a Poisson distribution) can be incorporated into models for predicting toehold switch performance. These efforts complement other recent works to predict the functional performance of toeholds with deep neural networks, and future directions will likely benefit from using a combination of these techniques in parallel.

### Data Sets/References
[1]A. A. Green, P. A. Silver, J. J. Collins, and P. Yin, “Toehold Switches: De-Novo-Designed Regulators of Gene Expression,” Cell, vol. 159, no. 4, pp. 925–939, Nov. 2014, doi: 10.1016/j.cell.2014.10.002.

[2]A. A. Green, J. Kim, D. Ma, P. A. Silver, J. J. Collins, and P. Yin, “Complex cellular logic computation using ribocomputing devices,” Nature, vol. 548, no. 7665, pp. 117–121, Aug. 2017, doi: 10.1038/nature23271.

[3] T. H. T. Chau and E. Y. Lee, “Development of cell-free platform-based toehold switch system for detection of IP-10 mRNA, an indicator for acute kidney allograft rejection diagnosis,” Clin Chim Acta, vol. 510, pp. 619–624, Nov. 2020, doi: 10.1016/j.cca.2020.08.034.

[4] N. M. Angenent-Mari, A. S. Garruss, L. R. Soenksen, G. Church, and J. J. Collins, “A deep learning approach to programmable RNA switches,” Nat Commun, vol. 11, no. 1, p. 5057, Oct. 2020, doi: 10.1038/s41467-020-18677-1.

[5] J. A. Valeri et al., “Sequence-to-function deep learning frameworks for engineered riboregulators,” Nat Commun, vol. 11, no. 1, p. 5058, Oct. 2020, doi: 10.1038/s41467-020-18676-2.

[6] M. A. English, R. V. Gayet, and J. J. Collins, “Designing Biological Circuits: Synthetic Biology Within the Operon Model and Beyond,” Annu. Rev. Biochem., vol. 90, no. 1, pp. 221–244, Jun. 2021, doi: 10.1146/annurev-biochem-013118-111914.

[7] T. L. Bailey, “DREME: motif discovery in transcription factor ChIP-seq data,” Bioinformatics, vol. 27, no. 12, pp. 1653–1659, Jun. 2011, doi: 10.1093/bioinformatics/btr261.

### Installation
### Usage
### References
### Contact Information
Aiden Riley; atriley@bu.edu\n
Juan Montezo; jjjarami@bu.edu\n
Eric South; esouth@bu.edu\n
<!-- New-HybridZone-Model README (HTML) -->

<h1>New-HybridZone-Model</h1>
<p><strong>Creating and implementing new mathematical model for gene flow across two species populations</strong></p>

<hr>
<h2>Disclaimer</h2>
<p>
  This project was done under the mentorship and collaboration of Pyron Lab. All data used was provided
  by Pyron Lab. 
</p>

<h2>TLDR</h2>
<p>
  <strong>New-HybridZone-Model</strong> creates a new mathematical model for the flow of genes across
  a hybrid zone where two species meet. It jointly models previously existing spatial cline and genomic 
  cline. The resulting model is a gradient-descent-optimized joint logistic regression approach, with additional
  tail options, that models gene flow based on spatial and genetic aspects to determine selection being
  performed against certain alleles. Our test case usage for the creation and assessment of this model was
  the Cheaha-Monticola Salamander hybrid zones across the eastern US.
</p>

<h2>Highlights</h2>
<ul>
  <li><strong>Mathematical Modeling</strong>: Developed a new </strong>joint spatial-genomic cline</strong> model to describe gene flow across
                                              hybrid zones</li>
  <li><strong>Machine Learning/Optimization</strong>: Implemented </strong>gradient descent optimization<strong> to fit a biologically informed </strong>logistic model<strong></li>
  <li><strong>Computational Skills</strong>: Implemented custom likelihood functions and optimizations, visualized clines</li>
  <li><strong>Creative Statistical Thinking</strong>: Applied logistic regression in a non-traditional, biologically motivated context</li>
</ul>

<h2>Data Visualization: Some examples!</h2>
<p>
  <img src="figures/extreme_identification.png" width="400" alt="Identification of Extreme Carbon Fluxes">
  <img src="figures/VIMP_interactions.png" width="400" alt="Random Forest Variable Importances and Interactions">
</p>

<h2>Skills</h2>
<ul>
  <li><strong>Languages</strong>: R</li>
  <li><strong>Libraries</strong>: ggplot2, optimx, plotly, parallel</li>
  <li><strong>Compute</strong>: Optimization and parallel processing</li>
</ul>

<h2>Why this project?</h2>
<p>
  Understanding the patterns of gene flow across the hybrid zone where two species meet and create hybrids provides insight
  into the types of selection being performed (e.g., directional or disruptive selection) on certain alleles. Previously,
  our ability to understand selection across hybrid zones was limited to patterns in allele frequencies across a spatial
  gradient, or across a gradient of genetic ancestry. We jointly modeled these approaches to provide a more robust understanding
  of hybrid gene flow and possible genetic selection on both a spatial and genetic ancestry level!
</p>

<h2>Workflow Overview</h2>
<pre>
carbonflux/
├── data/                   # Data on salamander species and their allele presence
├── analysis.R              # Implementation and analysis of previous clines and new joint model
├── jointcline.R            # The joint model and its optimization
├── aux_functions.R         # Previous models and their optimizations
├── figures/                # Some example figures of model implementation
└── README.md
</pre>

<h2>Products and Communication</h2>
<ul>
  <li><strong>Peer-Reviewed Publication</strong>: Manuscript in writing, stay tuned!</li>
  <li><strong>R Package for Joint Cline Implementation</strong>: Package in the works, stay tuned! </li>
</ul>

<h2>Contact</h2>

<p>
  For questions or collaboration: <strong>marleeyork2025@gmail.com</strong><br>
<p>


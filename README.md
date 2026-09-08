# Ichthyosaura_alpestris_life-cycle
Data and R scripts associated with the article: Fournier M, Louppe V, Toussaint-Lardé I, Pierson C, Chatar N, Tseng ZJ, Ferreira G, Kyriakouli C, Taillades M, Poloni L, Clavel J and Fabre AC. 2026.  Intespecific variation in life cycle impacts the morphology of Alpine newts. Integrative and Comparative Biology.

## Description of the data and file structure

# Data:

list_of_specimens: list of specimens used in this study, including information on subspecies, life cycle category, locality, date of collection, and presence of hyoid bones.

Linear_measurements: file containing linear measurement data used to quantify body morphology.
	> Linear_measurements_Ichthyosaura_alpestris_lc: raw data of linear measurements taken with a digital caliper (Mitutoyo Absolute AOS Digimatic)
	> df_svl: log10-transformed data for SVL only, used for the univariate analyses.
	> df_linear_filtered: log10-transformed data filtered to remove NAs, used for the multivariate analyses.
	> df_linear_mp: log10-transformed data without specimens undergoing metamorphosis, used for the disparity and modularity/integration analyses.

Landmarks: 3D landmarks quantifying the shape of the cranium and mandible.
	> Consensus_specimens: reference specimen meshes used to wrap the deformation along PC axes.
	> Ident: list of specimens with category information and centroid size values.
	> pts: Rdata files containing 3D landmark coordinates after GPA for each structure and for each bone separately.

# Scripts: scripts used to run the analyses (used in R statistical software v4.5.2).

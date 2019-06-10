This folder contains the computer program called Platypus accompanying the book "The Scaled Boundary Finite Element Method: Theory and Implementation" by Chongmin Song.

No technical support is provided. Comments on Platypus are welcome (Email: c.song@unsw.edu.au).

The computer program Platypus is protected by copyright, whether or not a copyright notice appears on the particular screen where the material is displayed.

The content of this folder, besides this file (Readme1st.txt), is as follows:
- subfolder "src": functions in Chapters 2 and 3 for the scaled boundary finite element analysis
- subfolder "examples": functions for reproducing the examples in Chapters 2 and 3
- subfolder "mesher": functions for automatic mesh generation and for reproducing the examples in Chapter 4. The public domain software PolyMesher (https://paulino.ce.gatech.edu/software.html) and DistMesh (http://persson.berkeley.edu/distmesh/) are required.

To use this computer program, all subfolders have to be added to the path of MATLAB. 

After a problem definitin file is executed to input the model, the analysis follows the flowchart below:
1. Set up global system of equations ([K]{u} = {P}) 
	SBFEMAssembly.m
		SElementSlnEigenMethod.m
			SElementCoeffMtx.m
			EleCoeff2NodeEle.m
			IsoElasMtrx.m
		AddNodalForces.m
		addSurfTraction.m
	
2. Solution using one or more of the 3 function
	SolverMode.m
	SolverNewmark.m
	SolverStatics.m
	
3. Post-processing
	SElementIntgConst.m
	SElementStrainMode2NodeEle.m
	SElementInDispStrain.m

Additionally, there are two plotting functions in the subfolder "src":
	PlotDeformedMesh.m
	PlotSBFEMesh.m


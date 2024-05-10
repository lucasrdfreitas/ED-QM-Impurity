(* ::Package:: *)

(* ::Section:: *)
(*Exact diagonalization for Quantum Magnets with an impurity in a honeycomb lattice *)


(* ::Text:: *)
(*The code is organized in two section, one for Adatom case and other for the Substitutional case. *)
(*The Hamiltonian and other definition are in the  "definitions.wl"  file, and here I wrote the code for computing the eigenvalues and the *)
(**)
(*I wrote a code for fixed magnetic field h and varying Kondo coupling JK (which is easily modified for fix JK and varying h)*)
(*A strategy to not compute the Hamiltonian matrix for each JK is to compute the Hamiltonian for the Impurity Himp for JK=1 and the remain H0 (e.g. XXZ+Zeeman),*)
(*	we compute once Himp, H0, and at each iteration of the JK loop we just compute the sum  H=H0+JK Himp (much more efficient).*)
(**)
(*Within each section of this file (Adatom or Substitution)  we specifically compute*)
(*	- the eigenvalues for first "k" levels*)
(*	- the eigenvector for the lowest energy level  -> for computing spin operator projection in the ground state near the impurity*)
(*	- the data  is stored in the folder "data" under the subfolder with the specific model (e.g. Kitaev or XXZ) and label by the parameters (e.g. JK values)*)
(*	- to see how to load the data see the notebook "figure.nb"*)
(*	*)
(*Coupling:  XXZ*)
(*dataName: KondoCoupling Range & eigenvalues k levels*)


(* ::Text:: *)
(*first load functions from "definitions.wl"  file*)


If[ \[Not]($FrontEnd===Null), SetDirectory[NotebookDirectory[]] ];
$FileName=If[$FrontEnd === Null, $InputFileName, NotebookFileName[] ];
Get[ FileNameJoin[{Directory[],"definitions.wl" }] ]


(* ::Section::Closed:: *)
(*Launching Kernel *)


If[ ($FrontEnd===Null),
	Print["Before Starting Kernels"];
	Needs["ClusterIntegration`"];
	(*kernels = LaunchKernels[SGE["micro4", 10]];*)
	Quiet[kernels = LaunchKernels[]];
	Print["Starting Kernels"];
];


(* ::Section:: *)
(*Substitution*)


(* ::Subsection:: *)
(*Preamble*)


If[ \[Not]($FrontEnd===Null), SetDirectory[NotebookDirectory[]] ];
$FileName=If[$FrontEnd === Null, $InputFileName, NotebookFileName[] ];
Get[ FileNameJoin[{Directory[],"definitions.wl" }] ]


(* ::Text:: *)
(*couplings and parameters for the system *)


systemDimensions = {{2,2}};
kitaev           = {0{-1,-1,-1}};
heisenberg       = {-{1,1,1}};
anisotropy       = {-.1 cvec}; 
bulkSpin         = {3/2};
impuritySpin     = {1/2};
gs               = {1};
KondoCoupling    = {.5};
hRange           = {0,4.1,.02};
parameters       = N@Tuples[{systemDimensions,kitaev,heisenberg,anisotropy,KondoCoupling,impuritySpin,gs,bulkSpin} ]; 
klevels          = 15; 
HamCoupling="XXZ_FM_SUB";
dataName=Module[{i,f,\[Delta],k,jk}, 	{i,f,\[Delta]}=hRange;  k=klevels;    jk=KondoCoupling;  {i,f,\[Delta],k,jk}=ToString/@{i,f,\[Delta],k,jk};					
StringReplace["h=Range[i,f,d]_JK=jk_k=k0", {"i"->i,"f"->f,"d"->\[Delta],"k0"->k,"jk"->jk}]  ];
dataName2=Module[{i,f,\[Delta],k,jk}, 	{i,f,\[Delta]}=hRange;  k=klevels;  jk=KondoCoupling;  {i,f,\[Delta],k,jk}=ToString/@{i,f,\[Delta],k,jk};				
StringReplace["h=Range[i,f,d]_JK=jk", {"i"->i,"f"->f,"d"->\[Delta],"k0"->k,"jk"->jk}]  ];

hRange=Range@@hRange;  
(*hRange=Sort[hRange~Join~{3.65,3.7,3.75,4.85,3.9,3.95,4.05,4.1,4.15,4.25,4.3,4.35,4.5}];   *) 
Print["h Range length= ",Length@hRange];



(* ::Subsection:: *)
(*Code -- eigenvalues*)


Print[];Print["Computing Eigenvalues"];Print[];


Module[{Lx,Ly,J,\[Lambda]n,Simp,K,JK,g,bulkSpin,H0,HJ,HI,HK,HZ,eValues,pathToMatrices,info,datapath},
	{{Lx,Ly},K,J,\[Lambda]n,JK,Simp,g,bulkSpin}=parameters[[1]]; {Lx,Ly}=Round@{Lx,Ly};eValues={};	
	datapath=dataPathTXT[dataName,HamCoupling,Simp,{Lx,Ly},dataFolder];	Print["Data path : ",datapath];

		Print["Computing Hamiltonian matrices "]; 		
		Print["    H Kitaev timing      =  ",round@AbsoluteTiming[HK=If[Norm[K]==0,0,N@substitutionKitaev[K,Simp,Lx,Ly,bulkSpin]] ][[1]]," sec" ];
		Print["    H Heisenberg timing  = ", round@AbsoluteTiming[HJ=If[Norm[J]==0,0,N@substitutionHeisenberg[J,\[Lambda]n,Simp,Lx,Ly,bulkSpin]]][[1]]," sec"  ];
		Print["    H Zeeman timing      =  ",round@AbsoluteTiming[HZ=N@substitutionZeeman[cvec,Simp,g,2 Lx Ly,bulkSpin]  ][[1]]," sec"  ];
		Print["    H imp timing         =  ",round@AbsoluteTiming[HI=N@substitutionImp[JK,Simp, Lx, Ly,bulkSpin]][[1]]," sec"  ];
		H0=HK+HJ+HI; Clear[HK,HJ,HI];
		Print[" "];Print["    Memory in use:  ",round@N[10^-9  MemoryInUse[]]," GB"  ]; Print[" "];
	
	(* Loop for compute the k-first eigenvalues of H0+JK HI  *)
	Print["Starting JK Loop"];Print["Loop timing=",NumberForm[round[AbsoluteTiming[	
	Do[Module[{Himp,ev,h },
		h=hRange[[j]];
		Himp=N[(H0+h HZ)];
		Print["    Eigenvalue timing=",NumberForm[round[AbsoluteTiming[ 
		ev=Sort@(-Eigenvalues[-Himp, klevels,
		Method -> {"Arnoldi","Criteria"->"RealPart","MaxIterations"->2000,"Tolerance"->10^-8}]);  ][[1]]/60],{\[Infinity],3}]," min -- saving data for  j=",j,"/",Length@hRange, "; h=",h "; " ];		(*Print["    Memory in use:  ",N[10^-9  MemoryInUse[] ] ," GB" ];*)
		dataAppend[datapath,{Norm[h],ev}]; 
		
],{j,1,Length@hRange}]  ][[1]]/60],{\[Infinity],3}]," min " ];
]


(* ::Subsection:: *)
(*Code -- Spin projection*)


Module[{Lx,Ly,J,\[Lambda]n,Simp,K,JK,g,H0,HJ,HI,HK,HZ,eValues,pathToMatrices,info,datapath,bulkSpin},
	{{Lx,Ly},K,J,\[Lambda]n,JK,Simp,g,bulkSpin}=parameters[[1]]; {Lx,Ly}=Round@{Lx,Ly};	
	datapath=dataPathTXT[StringJoin[dataName2,"_spin_components"],HamCoupling,Simp,{Lx,Ly},dataFolder];	Print["Data path : ",datapath]; 
	
	Print["Computing Hamiltonian matrices "]; 		
		Print["    H Kitaev timing      =  ",round@AbsoluteTiming[HK=If[Norm[K]==0,0,N@substitutionKitaev[K,Simp,Lx,Ly,bulkSpin]] ][[1]]," sec" ];
		Print["    H Heisenberg timing  = ", round@AbsoluteTiming[HJ=If[Norm[J]==0,0,N@substitutionHeisenberg[J,\[Lambda]n,Simp,Lx,Ly,bulkSpin]]][[1]]," sec"  ];
		Print["    H Zeeman timing      =  ",round@AbsoluteTiming[HZ=N@substitutionZeeman[cvec,Simp,g,2 Lx Ly,bulkSpin]  ][[1]]," sec"  ];
		Print["    H imp timing         =  ",round@AbsoluteTiming[HI=N@substitutionImp[JK,Simp, Lx, Ly,bulkSpin]][[1]]," sec"  ];
		H0=HK+HJ+HI; 
		Clear[HK,HJ,HI];

	Print["Starting JK Loop"];
	Print["Loop timing=",NumberForm[round[AbsoluteTiming[
	
Do[Module[{Himp,evec,h,simp,s1,s2 },
		h=hRange[[j]];
		Himp=N[(H0+h HZ)];
		Print["    Eigenvector timing=",NumberForm[round[AbsoluteTiming[
		evec=(-Eigensystem[-Himp, 1,
		Method -> {"Arnoldi","Criteria"->"RealPart","MaxIterations"->2000,"Tolerance"->10^-9}])[[2,1]];  ][[1]]/60],{\[Infinity],3}], " min ;   j=",j,"/",Length@hRange];
		simp = Table[Conjugate[evec] . spinImpOp[Simp,2 Lx Ly-1 ,bulkSpin][[\[Gamma]]] . evec,{\[Gamma],1,3}]; 
		s1   = Table[Conjugate[evec] . spinOp[Simp,1,2 Lx Ly -1,bulkSpin][[\[Gamma]]] . evec,{\[Gamma],1,3}];
		s2   = Table[Conjugate[evec] . spinOp[Simp,2,2 Lx Ly -1,bulkSpin][[\[Gamma]]] . evec,{\[Gamma],1,3}];
		dataAppend[datapath,
		Chop@{h,{simp . avec,simp . bvec,simp . cvec},{s1 . avec,s1 . bvec,s1 . cvec},{s2 . avec,s2 . bvec,s2 . cvec}}    ];
		 
],{j,1,Length@hRange}]  ][[1]]/60],{\[Infinity],3}]," min " ];

(*AbsoluteTiming@dataWrite[datapath,eValues];*)

]


(* ::Section::Closed:: *)
(*Closing Kernels*)


If[ ($FrontEnd===Null), 
	Print[];
	Print["Closing Kernels"];
	CloseKernels[];
];

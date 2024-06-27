(* ::Package:: *)

(* ::Section:: *)
(*Exact diagonalization for Quantum Magnets with an impurity in a honeycomb lattice *)


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
(*Adatom*)


(* ::Subsection:: *)
(*Preamble*)


If[ \[Not]($FrontEnd===Null), SetDirectory[NotebookDirectory[]] ];
$FileName=If[$FrontEnd === Null, $InputFileName, NotebookFileName[] ];
Get[ FileNameJoin[{Directory[],"definitions.wl" }] ]


(* ::Text:: *)
(*couplings and parameters for the system *)

Print["Mathematica version= ",$VersionNumber]; Print[]; 

systemDimensions = {{3,3}};
kitaev           = {{-1,-1,-1}};
heisenberg       = {0{1,1,1}};
anisotropy       = {0 cvec}; 
bulkSpin         = {1/2};
impuritySpin     = {1/2};
gs               = {1.5};
KondoCoupling    = {1.50011};  (* 0.5*)
hRange           = {0.01,1.51,.04};
parameters       = N@Tuples[{systemDimensions,kitaev,heisenberg,anisotropy,KondoCoupling,impuritySpin,gs,bulkSpin} ]; 
klevels          = 150; 
HamCoupling="XXZ_FM_ADA";
dataName=Module[{i,f,\[Delta],k,jk}, 	{i,f,\[Delta]}=hRange;  k=klevels;    jk=KondoCoupling;  {i,f,\[Delta],k,jk}=ToString/@{i,f,\[Delta],k,jk};					
StringReplace["h=Range[i,f,d]_JK=jk_k=k0", {"i"->i,"f"->f,"d"->\[Delta],"k0"->k,"jk"->jk}]  ];
dataName2=Module[{i,f,\[Delta],k,jk}, 	{i,f,\[Delta]}=hRange;  k=klevels;  jk=KondoCoupling;  {i,f,\[Delta],k,jk}=ToString/@{i,f,\[Delta],k,jk};				
StringReplace["h=Range[i,f,d]_JK=jk", {"i"->i,"f"->f,"d"->\[Delta],"k0"->k,"jk"->jk}]  ];

hRange=Range@@hRange;   
hRange=hRange[[17;;-1]];
Print[hRange];
Print["h Range length= ",Length@hRange];




(* ::Subsection:: *)
(*Code -- eigenvalues*)


Module[{Lx,Ly,J,\[Lambda]n,Simp,K,JK,g,bulkSpin,H0,HJ,HI,HK,HZ,eValues,pathToMatrices,info,datapath},
	{{Lx,Ly},K,J,\[Lambda]n,JK,Simp,g,bulkSpin}=parameters[[1]]; {Lx,Ly}=Round@{Lx,Ly};eValues={};	
	datapath=dataPathTXT[dataName,HamCoupling,Simp,{Lx,Ly},dataFolder];	Print["Data path : ",datapath];

	
		Print["Computing Hamiltonian matrices "]; 		
		Print["    H Kitaev timing      =  ",round@AbsoluteTiming[HK=If[Norm[K]==0,0,N@AdatomKitaev[K,Simp,Lx,Ly,bulkSpin]] ][[1]]," sec" ];
		Print["    H Heisenberg timing  = ", round@AbsoluteTiming[HJ=If[Norm[J]==0,0,N@AdatomHeisenberg[J,\[Lambda]n,Simp,Lx,Ly,bulkSpin]]][[1]]," sec"  ];
		Print["    H Zeeman timing      =  ",round@AbsoluteTiming[HZ=N@AdatomZeeman[cvec,Simp,g,2 Lx Ly,bulkSpin]  ][[1]]," sec"  ];
		Print["    H imp timing         =  ",round@AbsoluteTiming[HI=N@AdatomImp[JK,Simp,2 Lx Ly,bulkSpin]][[1]]," sec"  ];
		H0=HK+HJ+HI; Clear[HK,HJ,HI];
		
	(* Loop for compute the k-first eigenvalues of H0+JK HI  *)
	Print["Starting JK Loop"];Print["Loop timing=",NumberForm[round[AbsoluteTiming[	
	Do[Module[{Himp,ev,h },
		h=hRange[[j]];
		Himp=N[(H0+h HZ)];
		Print["    Eigenvalue timing=",NumberForm[round[AbsoluteTiming[ 
		ev=Sort@(-Eigenvalues[-Himp, klevels,
		Method -> {"Arnoldi","Criteria"->"RealPart","MaxIterations"->6000,"Tolerance"->10^-4}]);  ][[1]]/60],{\[Infinity],3}]," min -- saving data for  j=",j,"/",Length@hRange, "; h=",h ,"; " ];		(*Print["    Memory in use:  ",N[10^-9  MemoryInUse[] ] ," GB" ];*)
		dataAppend[datapath,{Norm[h],ev}]; 
		
],{j,1,Length@hRange}]  ][[1]]/60],{\[Infinity],3}]," min " ];

];


(* ::Subsection:: *)
(*Code --  eigenvectors for   Spin projection*)


Print[];Print["Computing Eigenvectors"];Print[];


(*Module[{Lx,Ly,J,\[Lambda]n,Simp,K,JK,bulkSpin,g,H0,HJ,HI,HK,HZ,eValues,pathToMatrices,info,datapath},
	{{Lx,Ly},K,J,\[Lambda]n,JK,Simp,g,bulkSpin}=parameters[[1]]; {Lx,Ly}=Round@{Lx,Ly};eValues={};	
	datapath=dataPathTXT[StringJoin[dataName2,"_spin_components"],HamCoupling,Simp,{Lx,Ly},dataFolder];	Print["Data path : ",datapath];

		Print["Computing Hamiltonian matrices "]; 		
		Print["    H Kitaev timing      =  ",round@AbsoluteTiming[HK=If[Norm[K]==0,0,N@AdatomKitaev[K,Simp,Lx,Ly,bulkSpin]] ][[1]]," sec" ];
		Print["    H Heisenberg timing  = ", round@AbsoluteTiming[HJ=If[Norm[J]==0,0,N@AdatomHeisenberg[J,\[Lambda]n,Simp,Lx,Ly,bulkSpin]]][[1]]," sec"  ];
		Print["    H Zeeman timing      =  ",round@AbsoluteTiming[HZ=N@AdatomZeeman[cvec,Simp,g,2 Lx Ly,bulkSpin]  ][[1]]," sec"  ];
		Print["    H imp timing         =  ",round@AbsoluteTiming[HI=N@AdatomImp[JK,Simp,2 Lx Ly,bulkSpin]][[1]]," sec"  ];
		H0=HK+HJ+HI; 		Clear[HK,HJ,HI];
		
	Print["Starting JK Loop"];
	Print["Loop timing=",NumberForm[round[AbsoluteTiming[
	
	Do[Module[{Himp,evec,h,simp,s1,s2,s3,s4 },
		h=hRange[[j]];
		Himp=N[(H0+h HZ)];
		Print["    Eigenvector timing=",NumberForm[round[AbsoluteTiming[
		evec=(-Eigensystem[-Himp, 1,
		Method -> {"Arnoldi","Criteria"->"RealPart","MaxIterations"->2000,"Tolerance"->10^-9}])[[2,1]];  ][[1]]/60],{\[Infinity],3}], " min ;   j=",j,"/",Length@hRange];
		simp = Table[Conjugate[evec] . spinImpOp[Simp,2 Lx Ly ,bulkSpin][[\[Gamma]]] . evec,{\[Gamma],1,3}]; 
		s1   = Table[Conjugate[evec] . spinOp[Simp,1,2 Lx Ly,bulkSpin ][[\[Gamma]]] . evec,{\[Gamma],1,3}];
		s2   = Table[Conjugate[evec] . spinOp[Simp,2,2 Lx Ly,bulkSpin ][[\[Gamma]]] . evec,{\[Gamma],1,3}];
		s3   = Table[Conjugate[evec] . spinOp[Simp,3,2 Lx Ly,bulkSpin ][[\[Gamma]]] . evec,{\[Gamma],1,3}];
		s4   = Table[Conjugate[evec] . spinOp[Simp,3,2 Lx Ly,bulkSpin ][[\[Gamma]]] . evec,{\[Gamma],1,3}];
		dataAppend[datapath,
		Chop@{h,{simp . avec,simp . bvec,simp . cvec},{s1 . avec,s1 . bvec,s1 . cvec},{s2 . avec,s2 . bvec,s2 . cvec},{s3 . avec,s3 . bvec,s3 . cvec},{s4 . avec,s4 . bvec,s4 . cvec}     }    ];
		 
],{j,1,Length@hRange}]  ][[1]]/60],{\[Infinity],3}]," min " ];

(*AbsoluteTiming@dataWrite[datapath,eValues];*)

]*)


(* ::Section::Closed:: *)
(*Closing Kernels*)


If[ ($FrontEnd===Null), 
	Print[];
	Print["Closing Kernels"];
	CloseKernels[];
];

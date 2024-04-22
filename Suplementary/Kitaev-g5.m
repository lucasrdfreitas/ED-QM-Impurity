(* ::Package:: *)

(* ::Section:: *)
(*Exact diagonalization for Quantum Magnets with an impurity in a honeycomb lattice *)


(* ::Text:: *)
(*first load functions from "definitions.wl"  file*)


Quit[]


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


systemDimensions = {{3,3}};
kitaev           = {{-1,-1,-1}};
heisenberg       = {0{1,1,1}};
anisotropy       = {0 cvec}; 
impuritySpin     = {3/2};
gs               = {5};
KondoCoupling    = {.5};
hRange           = {0,1.2,1.2/40};
parameters       = N@Tuples[{systemDimensions,kitaev,heisenberg,anisotropy,KondoCoupling,impuritySpin,gs} ]; 
klevels          = 10; 
HamCoupling="Kitaev_FM_ADA";
dataName=Module[{i,f,\[Delta],k,jk}, 	{i,f,\[Delta]}=hRange;  k=klevels;    jk=KondoCoupling;  {i,f,\[Delta],k,jk}=ToString/@{i,f,\[Delta],k,jk};					
StringReplace["h=Range[i,f,d]_JK=jk_k=k0", {"i"->i,"f"->f,"d"->\[Delta],"k0"->k,"jk"->jk}]  ];
dataName2=Module[{i,f,\[Delta],k,jk}, 	{i,f,\[Delta]}=hRange;  k=klevels;  jk=KondoCoupling;  {i,f,\[Delta],k,jk}=ToString/@{i,f,\[Delta],k,jk};				
StringReplace["h=Range[i,f,d]_JK=jk", {"i"->i,"f"->f,"d"->\[Delta],"k0"->k,"jk"->jk}]  ];

hRange=Range@@hRange;      Print["h Range length= ",Length@hRange];


(* ::Subsection::Closed:: *)
(*Code -- save matrices*)


(* ::Text:: *)
(*to use the cluster in a older Mathematica version, I save the Hamiltonian matrices on my laptop*)


(*Module[{Lx,Ly,J,\[Lambda]n,Simp,K,h,g,H0,HK,HJ,HZ,HI,eValues,info,path,Huncompressed}, 
{{Lx,Ly},K,J,\[Lambda]n,h,Simp,g}=parameters[[1]]; 
{Lx,Ly}=Round@{Lx,Ly};
info=StringReplace["simp=X_h=Y",{"X"->ToString@Simp,"Y"->ToString@N[Round[1000 Norm@h]/1000]}];
path=dataPath[#,HamCoupling,Simp,{Lx,Ly},dataFolder]&@(StringJoin@{{"Hmatrices_"},{info}}); 

If[ FindFile[StringJoin[path,".zip"]]===$Failed,
	Print["Compressed data not found - Computing Hamiltonian matrix "];
	Print[];
	Print["    Memory in use:  ",     N[10^-9  MemoryInUse[] ]  ];
	Print["H Kitaev timing=",     AbsoluteTiming[HK=If[Norm[K]==0,0,N@AdatomKitaev[K,Simp,Lx,Ly]]; Dimensions@HK ] ];
	Print["H Heisenberg timing=", AbsoluteTiming[HJ=If[Norm[J]==0,0,N@AdatomHeisenberg[J,\[Lambda]n,Simp,Lx,Ly]];Dimensions@HJ] ];
	Print["H Zeeman timing=",     AbsoluteTiming[HZ=If[Norm[h]==0,0,N@AdatomZeeman[h,Simp,g,2 Lx Ly]];Dimensions@HZ] ];
	Print["H imp timing=",        AbsoluteTiming[HI=N@AdatomImp[1,Simp,2 Lx Ly]; Dimensions@HI] ];
	H0=HK+HJ+HZ;
	Clear[HK,HJ,HZ];
	Print["    Memory in use:  ",     N[10^-9  MemoryInUse[] ]  ];
	Print["Compress time=",       AbsoluteTiming[ dataZipExport[path,{H0,HI}] ;]     ];
	,
	Print["Compressed data found - Skipping matrix calculation"]
];

];*)


(* ::Subsection:: *)
(*Code -- eigenvalues*)


Pause[200];


Print[];Print["Computing Eigenvalues"];Print[];


Module[{Lx,Ly,J,\[Lambda]n,Simp,K,JK,g,H0,HJ,HI,HK,HZ,eValues,pathToMatrices,info,datapath},
	{{Lx,Ly},K,J,\[Lambda]n,JK,Simp,g}=parameters[[1]]; {Lx,Ly}=Round@{Lx,Ly};eValues={};	
	datapath=dataPathTXT[dataName,HamCoupling,Simp,{Lx,Ly},dataFolder];	Print["Data path : ",datapath];

	(* If the Hamiltonian matrix were already computed then load it, otherwise compute it*)
	info=StringReplace["simp=X_JK=Y",{"X"->ToString@Simp,"Y"->ToString@N[Round[1000 Norm@JK]/1000]}];		pathToMatrices=dataPath[#,HamCoupling,Simp,{Lx,Ly},dataFolder]&@(StringJoin@{{"Hmatrices_"},{info}});
	If[ Length@FileNames[pathToMatrices]!= 0,	
		Print["Loading matrices -- Uncompress time=",AbsoluteTiming[ 
		{H0,HI}=dataZipImport[pathToMatrices]; ]];   	
		,
		Print["Computing Hamiltonian matrices "]; 		
		Print["    H Kitaev timing      =  ",round@AbsoluteTiming[HK=If[Norm[K]==0,0,N@AdatomKitaev[K,Simp,Lx,Ly]] ][[1]]," sec" ];
		Print["    H Heisenberg timing  = ", round@AbsoluteTiming[HJ=If[Norm[J]==0,0,N@AdatomHeisenberg[J,\[Lambda]n,Simp,Lx,Ly]]][[1]]," sec"  ];
		Print["    H Zeeman timing      =  ",round@AbsoluteTiming[HZ=N@AdatomZeeman[cvec,Simp,g,2 Lx Ly]  ][[1]]," sec"  ];
		Print["    H imp timing         =  ",round@AbsoluteTiming[HI=N@AdatomImp[JK,Simp,2 Lx Ly]][[1]]," sec"  ];
		H0=HK+HJ+HI; Clear[HK,HJ,HI];
	]; Print[" "];Print["    Memory in use:  ",round@N[10^-9  MemoryInUse[]]," GB"  ]; Print[" "];
	
	(* Loop for compute the k-first eigenvalues of H0+JK HI  *)
	Print["Starting JK Loop"];Print["Loop timing=",NumberForm[round[AbsoluteTiming[	
	Do[Module[{Himp,ev,h },
		h=hRange[[j]];
		Himp=N[(H0+h HZ)];
		Print["    Eigenvalue timing=",NumberForm[round[AbsoluteTiming[ 
		ev=Sort@(-Eigenvalues[-Himp, klevels,
		Method -> {"Arnoldi","Criteria"->"RealPart","MaxIterations"->4000,"Tolerance"->10^-9}]);  ][[1]]/60],{\[Infinity],3}]," min -- saving data for  j=",j,"/",Length@hRange, "; h=",h "; " ];		(*Print["    Memory in use:  ",N[10^-9  MemoryInUse[] ] ," GB" ];*)
		dataAppend[datapath,{Norm[h],ev}]; 
		
],{j,1,Length@hRange}]  ][[1]]/60],{\[Infinity],3}]," min " ];
]


(* ::Subsection::Closed:: *)
(*Code --  eigenvectors for   Spin projection*)


Module[{Lx,Ly,J,\[Lambda]n,Simp,K,h,g,H0,HJ,HI,HK,HZ,eValues,pathToMatrices,info,datapath},
	{{Lx,Ly},K,J,\[Lambda]n,h,Simp,g}=parameters[[1]]; {Lx,Ly}=Round@{Lx,Ly};
	
	datapath=dataPathTXT[StringJoin[dataName2,"_spin_components"],HamCoupling,Simp,{Lx,Ly},dataFolder];	Print["Data path : ",datapath];

	(* If the Hamiltonian matrix we already computed, then load it, otherwise compute it*)
	info=StringReplace["simp=X_h=Y",{"X"->ToString@Simp,"Y"->ToString@N[Round[1000 Norm@h]/1000]}];		pathToMatrices=dataPath[#,HamCoupling,Simp,{Lx,Ly},dataFolder]&@(StringJoin@{{"Hmatrices_"},{info}});
	
	If[ Length@FileNames[pathToMatrices]!= 0,Print["Loading matrices -- Uncompress time=",AbsoluteTiming[ {H0,HI}=dataZipImport[pathToMatrices]; ]               ];   ,	Print["Computing Hamiltonian matrices "]; 		
		Print["    H Kitaev timing      =  ",     round@AbsoluteTiming[HK=If[Norm[K]==0,0,N@AdatomKitaev[K,Simp,Lx,Ly]] ][[1]]," sec" ];
		Print["    H Heisenberg timing  = ", round@AbsoluteTiming[HJ=If[Norm[J]==0,0,N@AdatomHeisenberg[J,\[Lambda]n,Simp,Lx,Ly]]][[1]]," sec"  ];
		Print["    H Zeeman timing      =  ",     round@AbsoluteTiming[HZ=If[Norm[h]==0,0,N@AdatomZeeman[h,Simp,g,2 Lx Ly]]][[1]]," sec"  ];
		Print["    H imp timing         =  ",        round@AbsoluteTiming[HI=N@AdatomImp[1,Simp,2 Lx Ly]][[1]]," sec"  ];
		H0=HK+HJ+HZ; Clear[HK,HJ,HZ];
	]; Print[" "];Print["    Memory in use:  ",round@N[10^-9  MemoryInUse[]]," GB"  ]; Print[" "];
	

	Print["Starting JK Loop"];
	Print["Loop timing=",NumberForm[round[AbsoluteTiming[
	
	Do[Module[{Himp,evec,JK,simp,s1,s2 },
		JK=KondoCouplings[[j]];
		Himp=N[  (H0+JK HI)];
		Print["    Eigenvector timing=",NumberForm[round[AbsoluteTiming[
		evec=(-Eigensystem[-Himp, 1,
		Method -> {"Arnoldi","Criteria"->"RealPart","MaxIterations"->2000,"Tolerance"->10^-9}])[[2,1]];  ][[1]]/60],{\[Infinity],3}], " min ;   j=",j,"/",Length@KondoCouplings];
		simp = Table[Conjugate[evec] . spinImpOp[Simp,2 Lx Ly ][[\[Gamma]]] . evec,{\[Gamma],1,3}]; 
		s1   = Table[Conjugate[evec] . spinOp[Simp,1,2 Lx Ly ][[\[Gamma]]] . evec,{\[Gamma],1,3}];
		s2   = Table[Conjugate[evec] . spinOp[Simp,2,2 Lx Ly ][[\[Gamma]]] . evec,{\[Gamma],1,3}];
		dataAppend[datapath,
		Chop@{JK,{simp . avec,simp . bvec,simp . cvec},{s1 . avec,s1 . bvec,s1 . cvec},{s2 . avec,s2 . bvec,s2 . cvec}}    ];
		 
],{j,1,Length@KondoCouplings}]  ][[1]]/60],{\[Infinity],3}]," min " ];

(*AbsoluteTiming@dataWrite[datapath,eValues];*)

]


(* ::Section:: *)
(*Substitution*)


(* ::Subsection:: *)
(*Preamble*)


Pause[200];


If[ \[Not]($FrontEnd===Null), SetDirectory[NotebookDirectory[]] ];
$FileName=If[$FrontEnd === Null, $InputFileName, NotebookFileName[] ];
Get[ FileNameJoin[{Directory[],"definitions.wl" }] ]


(* ::Text:: *)
(*couplings and parameters for the system *)


systemDimensions = {{3,3}};
kitaev           = {{-1,-1,-1}};
heisenberg       = {0{1,1,1}};
anisotropy       = {0 cvec}; 
impuritySpin     = {3/2};
gs               = {5};
KondoCoupling    = {.5};
hRange           = {0,1.2,1.2/40};
parameters       = N@Tuples[{systemDimensions,kitaev,heisenberg,anisotropy,KondoCoupling,impuritySpin,gs} ]; 
klevels          = 10; 
HamCoupling="Kitaev_FM_SUB";
dataName=Module[{i,f,\[Delta],k,jk}, 	{i,f,\[Delta]}=hRange;  k=klevels;    jk=KondoCoupling;  {i,f,\[Delta],k,jk}=ToString/@{i,f,\[Delta],k,jk};					
StringReplace["h=Range[i,f,d]_JK=jk_k=k0", {"i"->i,"f"->f,"d"->\[Delta],"k0"->k,"jk"->jk}]  ];
dataName2=Module[{i,f,\[Delta],k,jk}, 	{i,f,\[Delta]}=hRange;  k=klevels;  jk=KondoCoupling;  {i,f,\[Delta],k,jk}=ToString/@{i,f,\[Delta],k,jk};				
StringReplace["h=Range[i,f,d]_JK=jk", {"i"->i,"f"->f,"d"->\[Delta],"k0"->k,"jk"->jk}]  ];

hRange=Range@@hRange;      Print["h Range length= ",Length@hRange];


(* ::Subsection::Closed:: *)
(*Code -- save matrices*)


(* ::Text:: *)
(*to use the cluster in a older Mathematica version, I save the Hamiltonian matrices on my laptop*)


(*Module[{Lx,Ly,J,\[Lambda]n,Simp,K,h,g,H0,HK,HJ,HZ,HI,eValues,info,path,Huncompressed}, 
{{Lx,Ly},K,J,\[Lambda]n,h,Simp,g}=parameters[[1]]; 
{Lx,Ly}=Round@{Lx,Ly};
info=StringReplace["simp=X_h=Y",{"X"->ToString@Simp,"Y"->ToString@N[Round[1000 Norm@h]/1000]}];
path=dataPath[#,HamCoupling,Simp,{Lx,Ly},dataFolder]&@(StringJoin@{{"Hmatrices_"},{info}}); 

If[ FindFile[StringJoin[path,".zip"]]===$Failed,
	Print["Compressed data not found - Computing Hamiltonian matrix "];
	Print[];
	Print["    Memory in use:  ",     N[10^-9  MemoryInUse[] ]  ];
	Print["H Kitaev timing=",     AbsoluteTiming[HK=If[Norm[K]==0,0,N@substitutionKitaev[K,Simp,Lx,Ly]]; Dimensions@HK ] ];
	Print["H Heisenberg timing=", AbsoluteTiming[HJ=If[Norm[J]==0,0,N@substitutionHeisenberg[J,\[Lambda]n,Simp,Lx,Ly]];Dimensions@HJ] ];
	Print["H Zeeman timing=",     AbsoluteTiming[HZ=If[Norm[h]==0,0,N@substitutionZeeman[h,Simp,g,2 Lx Ly]];Dimensions@HZ] ];
	Print["H imp timing=",        AbsoluteTiming[HI=N@substitutionImp[1,Simp,2 Lx Ly]; Dimensions@HI] ];
	H0=HK+HJ+HZ;
	Clear[HK,HJ,HZ];
	Print["    Memory in use:  ",     N[10^-9  MemoryInUse[] ]  ];
	Print["Compress time=",       AbsoluteTiming[ dataZipExport[path,{H0,HI}] ;]     ];
	,
	Print["Compressed data found - Skipping matrix calculation"]
];

];*)


(* ::Subsection:: *)
(*Code -- eigenvalues*)


Print[];Print["Computing Eigenvalues"];Print[];


Module[{Lx,Ly,J,\[Lambda]n,Simp,K,JK,g,H0,HJ,HI,HK,HZ,eValues,pathToMatrices,info,datapath},
	{{Lx,Ly},K,J,\[Lambda]n,JK,Simp,g}=parameters[[1]]; {Lx,Ly}=Round@{Lx,Ly};eValues={};	
	datapath=dataPathTXT[dataName,HamCoupling,Simp,{Lx,Ly},dataFolder];	Print["Data path : ",datapath];

	(* If the Hamiltonian matrix were already computed then load it, otherwise compute it*)
	info=StringReplace["simp=X_JK=Y",{"X"->ToString@Simp,"Y"->ToString@N[Round[1000 Norm@JK]/1000]}];		pathToMatrices=dataPath[#,HamCoupling,Simp,{Lx,Ly},dataFolder]&@(StringJoin@{{"Hmatrices_"},{info}});
	If[ Length@FileNames[pathToMatrices]!= 0,	
		Print["Loading matrices -- Uncompress time=",AbsoluteTiming[ 
		{H0,HI}=dataZipImport[pathToMatrices]; ]];   	
		,
		Print["Computing Hamiltonian matrices "]; 		
		Print["    H Kitaev timing      =  ",round@AbsoluteTiming[HK=If[Norm[K]==0,0,N@substitutionKitaev[K,Simp,Lx,Ly]] ][[1]]," sec" ];
		Print["    H Heisenberg timing  = ", round@AbsoluteTiming[HJ=If[Norm[J]==0,0,N@substitutionHeisenberg[J,\[Lambda]n,Simp,Lx,Ly]]][[1]]," sec"  ];
		Print["    H Zeeman timing      =  ",round@AbsoluteTiming[HZ=N@substitutionZeeman[cvec,Simp,g,2 Lx Ly]  ][[1]]," sec"  ];
		Print["    H imp timing         =  ",round@AbsoluteTiming[HI=N@substitutionImp[JK,Simp, Lx, Ly]][[1]]," sec"  ];
		H0=HK+HJ+HI; Clear[HK,HJ,HI];
	]; Print[" "];Print["    Memory in use:  ",round@N[10^-9  MemoryInUse[]]," GB"  ]; Print[" "];
	
	(* Loop for compute the k-first eigenvalues of H0+JK HI  *)
	Print["Starting JK Loop"];Print["Loop timing=",NumberForm[round[AbsoluteTiming[	
	Do[Module[{Himp,ev,h },
		h=hRange[[j]];
		Himp=N[(H0+h HZ)];
		Print["    Eigenvalue timing=",NumberForm[round[AbsoluteTiming[ 
		ev=Sort@(-Eigenvalues[-Himp, klevels,
		Method -> {"Arnoldi","Criteria"->"RealPart","MaxIterations"->3000,"Tolerance"->10^-8}]);  ][[1]]/60],{\[Infinity],3}]," min -- saving data for  j=",j,"/",Length@hRange, "; h=",h "; " ];		(*Print["    Memory in use:  ",N[10^-9  MemoryInUse[] ] ," GB" ];*)
		dataAppend[datapath,{Norm[h],ev}]; 
		
],{j,1,Length@hRange}]  ][[1]]/60],{\[Infinity],3}]," min " ];
]


(* ::Subsection:: *)
(*Code -- Spin projection*)


Module[{Lx,Ly,J,\[Lambda]n,Simp,K,h,g,H0,HJ,HI,HK,HZ,eValues,pathToMatrices,info,datapath},
	{{Lx,Ly},K,J,\[Lambda]n,h,Simp,g}=parameters[[1]]; {Lx,Ly}=Round@{Lx,Ly};
	
	datapath=dataPathTXT[StringJoin[dataName2,"_spin_components"],HamCoupling,Simp,{Lx,Ly},dataFolder];	Print["Data path : ",datapath];

	(* If the Hamiltonian matrix we already computed, then load it, otherwise compute it*)
	info=StringReplace["simp=X_h=Y",{"X"->ToString@Simp,"Y"->ToString@N[Round[1000 Norm@h]/1000]}];		pathToMatrices=dataPath[#,HamCoupling,Simp,{Lx,Ly},dataFolder]&@(StringJoin@{{"Hmatrices_"},{info}});
	
	If[ Length@FileNames[pathToMatrices]!= 0,Print["Loading matrices -- Uncompress time=",AbsoluteTiming[ {H0,HI}=dataZipImport[pathToMatrices]; ]               ];   ,	Print["Computing Hamiltonian matrices "]; 		
		Print["    H Kitaev timing      =  ",     round@AbsoluteTiming[HK=If[Norm[K]==0,0,N@substitutionKitaev[K,Simp,Lx,Ly]] ][[1]]," sec" ];
		Print["    H Heisenberg timing  = ", round@AbsoluteTiming[HJ=If[Norm[J]==0,0,N@substitutionHeisenberg[J,\[Lambda]n,Simp,Lx,Ly]]][[1]]," sec"  ];
		Print["    H Zeeman timing      =  ",     round@AbsoluteTiming[HZ=If[Norm[h]==0,0,N@substitutionZeeman[h,Simp,g,2 Lx Ly]]][[1]]," sec"  ];
		Print["    H imp timing         =  ",        round@AbsoluteTiming[HI=N@substitutionImp[1,Simp, Lx, Ly]][[1]]," sec"  ];
		H0=HK+HJ+HZ; Clear[HK,HJ,HZ];
	]; Print[" "];Print["    Memory in use:  ",round@N[10^-9  MemoryInUse[]]," GB"  ]; Print[" "];
	

	Print["Starting JK Loop"];
	Print["Loop timing=",NumberForm[round[AbsoluteTiming[
	
	Do[Module[{Himp,evec,JK,simp,s1,s2 },
		JK=KondoCouplings[[j]];
		Himp=N[  (H0+JK HI)];
		Print["    Eigenvector timing=",NumberForm[round[AbsoluteTiming[
		evec=(-Eigensystem[-Himp, 1,
		Method -> {"Arnoldi","Criteria"->"RealPart","MaxIterations"->2000,"Tolerance"->10^-9}])[[2,1]];  ][[1]]/60],{\[Infinity],3}], " min ;   j=",j,"/",Length@KondoCouplings];
		simp = Table[Conjugate[evec] . spinImpOp[Simp,2 Lx Ly-1 ][[\[Gamma]]] . evec,{\[Gamma],1,3}]; 
		s1   = Table[Conjugate[evec] . spinOp[Simp,1,2 Lx Ly -1][[\[Gamma]]] . evec,{\[Gamma],1,3}];
		s2   = Table[Conjugate[evec] . spinOp[Simp,2,2 Lx Ly -1][[\[Gamma]]] . evec,{\[Gamma],1,3}];
		dataAppend[datapath,
		Chop@{JK,{simp . avec,simp . bvec,simp . cvec},{s1 . avec,s1 . bvec,s1 . cvec},{s2 . avec,s2 . bvec,s2 . cvec}}    ];
		 
],{j,1,Length@KondoCouplings}]  ][[1]]/60],{\[Infinity],3}]," min " ];

(*AbsoluteTiming@dataWrite[datapath,eValues];*)

]


(* ::Section::Closed:: *)
(*Closing Kernels*)


If[ ($FrontEnd===Null), 
	Print[];
	Print["Closing Kernels"];
	CloseKernels[];
];

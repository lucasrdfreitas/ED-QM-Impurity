(* ::Package:: *)

(* ::Text:: *)
(*Adatom and Substitution *)
(**)
(*-- varying Kondo coupling JK, to compute*)
(*	- the eigenvalues for first "k" levels*)
(*	- the eigenvector for the ground state -> spin operator projection*)
(*	*)
(*Coupling:  XXZ*)
(*dataName: KondoCoupling Range & eigenvalues k levels*)


(* ::Subsection:: *)
(*preamble*)


(* ::Text:: *)
(*load functions from "definitions.wl"  file*)


If[ \[Not]($FrontEnd===Null), SetDirectory[NotebookDirectory[]] ];
$FileName=If[$FrontEnd === Null, $InputFileName, NotebookFileName[] ];
Get[ FileNameJoin[{Directory[],"definitions.wl" }] ]


(* ::Text:: *)
(*couplings and parameters for the system *)


systemDimensions = {{3,3}};
kitaev           = {0{-1,-1,-1}};
heisenberg       = {-{1,1,1}};
anisotropy       = {-.1{1,1,1}};
hfields          = {0.5 cvec};
impuritySpin     = {1/2};
gs               = {1};
KondoCouplings   = {0,1,.02};
parameters       = N@Tuples[{systemDimensions,kitaev,heisenberg,anisotropy,hfields,impuritySpin,gs} ];

(* klevels = number of eigenvalues to compute*)

klevels          = 2;

(* name to save the files : *)

HamCoupling="XXZ_FM";
dataName=Module[{i,f,\[Delta],k}, 				{i,f,\[Delta]}=KondoCouplings;  k=klevels;    {i,f,\[Delta],k}=ToString/@{i,f,\[Delta],k};				StringReplace["JK=Range[i,f,d]", {"i"->i,"f"->f,"d"->\[Delta],"k0"->k}]  ];

KondoCouplings=Range@@KondoCouplings;      Print["\!\(\*SubscriptBox[\(J\), \(K\)]\) couplings length= ",Length@KondoCouplings];


(* ::Subsection::Closed:: *)
(*Launching Kernel *)


If[ ($FrontEnd===Null),
	Print["Before Starting Kernels"];
	Needs["ClusterIntegration`"];
	(*kernels = LaunchKernels[SGE["micro4", 10]];*)
	Quiet[kernels = LaunchKernels[]];
	Print["Starting Kernels"];
];


(* ::Subsection::Closed:: *)
(*code -- save matrices*)


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


Print[];Print["Computing Eigenvalues"];Print[];


Module[{Lx,Ly,J,\[Lambda]n,Simp,K,h,g,H0,HJ,HI,HK,HZ,eValues,pathToMatrices,info,datapath},
{{Lx,Ly},K,J,\[Lambda]n,h,Simp,g}=parameters[[1]]; 
{Lx,Ly}=Round@{Lx,Ly};eValues={};
	
	datapath=dataPath[dataName,HamCoupling,Simp,{Lx,Ly},dataFolder];
	Print["Data path : ",datapath];

	(* If the Hamiltonian matrix we already computed, then load it, otherwise compute it*)
	info=StringReplace["simp=X_h=Y",{"X"->ToString@Simp,"Y"->ToString@N[Round[1000 Norm@h]/1000]}];		pathToMatrices=dataPath[#,HamCoupling,Simp,{Lx,Ly},dataFolder]&@(StringJoin@{{"Hmatrices_"},{info}});
	
	If[ Length@FileNames[pathToMatrices]!= 0,
		Print["Loading matrices -- Uncompress time=",AbsoluteTiming[ 
		{H0,HI}=dataZipImport[path]; ]               ];   
	,
		Print["Computing Hamiltonian matrices "]; 		
		Print["H Kitaev timing=",     AbsoluteTiming[HK=If[Norm[K]==0,0,N@AdatomKitaev[K,Simp,Lx,Ly]]; Dimensions@HK ] ];
		Print["H Heisenberg timing=", AbsoluteTiming[HJ=If[Norm[J]==0,0,N@AdatomHeisenberg[J,\[Lambda]n,Simp,Lx,Ly]];Dimensions@HJ] ];
		Print["H Zeeman timing=",     AbsoluteTiming[HZ=If[Norm[h]==0,0,N@AdatomZeeman[h,Simp,g,2 Lx Ly]];Dimensions@HZ] ];
		Print["H imp timing=",        AbsoluteTiming[HI=N@AdatomImp[1,Simp,2 Lx Ly]; Dimensions@HI] ];
		H0=HK+HJ+HZ; Clear[HK,HJ,HZ];
	]; Print[" "];Print["    Memory in use:  ",N[10^-9  MemoryInUse[]," GB" ]  ];
	
	Print["Starting JK Loop"];Print["Loop timing=",AbsoluteTiming[	
	Do[Module[{Himp,ev,JK },
		JK=KondoCouplings[[j]];
		Himp=N[  (H0+JK HI)];
		Print["    Eigenvalue timing=",AbsoluteTiming[(*ev=Sort@Eigenvalues[Himp,2 klevels];*)
		ev=Sort@(-Eigenvalues[-Himp, klevels,
		Method -> {"Arnoldi","Criteria"->"RealPart","MaxIterations"->1000,"Tolerance"->10^-9(*, "Shift"->-.7*)}]);  ]];
		
		dataAppend[datapath,{JK,ev}];
		(*AppendTo[eValues,{JK,ev}];*)
		Print["    j=",j,"/",Length@KondoCouplings];
		Print["    Memory in use:  ",N[10^-9  MemoryInUse[] ] ," GB" ];
],{j,1,Length@KondoCouplings}]] ];

(*AbsoluteTiming@dataWrite[datapath,eValues];*)

]


(* ::Subsection:: *)
(*Code -- Spin projection*)


Module[{Lx,Ly,J,\[Lambda]n,Simp,K,h,g,H0,HJ,HI,HK,HZ,eValues,path,info,datapath},
{{Lx,Ly},K,J,\[Lambda]n,h,Simp,g}=parameters[[1]]; 
{Lx,Ly}=Round@{Lx,Ly};
	
	datapath=dataPath[StringJoin[dataName,"_spin_components"],HamCoupling,Simp,{Lx,Ly},dataFolder];
	Print["Data path : ",datapath];

	HK=If[Norm[K]==0,0,N@AdatomKitaev[K,Simp,Lx,Ly]];
	HJ=If[Norm[J]==0,0,N@AdatomHeisenberg[J,\[Lambda]n,Simp,Lx,Ly]];
	HZ=If[Norm[h]==0,0,N@AdatomZeeman[h,Simp,g,2 Lx Ly]];
	HI=N@AdatomImp[1,Simp,2 Lx Ly]; 
	H0=HK+HJ+HZ;
	Clear[HK,HJ,HZ];
	(*info=StringReplace["simp=X_h=Y",{"X"->ToString@Simp,"Y"->ToString@N[Round[1000 Norm@h]/1000]}];
	path=dataPath[#,HamCoupling,Simp,{Lx,Ly},dataFolder]&@(StringJoin@{{"Hmatrices_"},{info}}); 
	Print["Uncompress time=", AbsoluteTiming[ {H0,HI}=dataZipImport[path]; ]];
	*)Print[" "];
	Print["    Memory in use:  ",N[10^-9  MemoryInUse[] ]  ];

	Print["Starting JK Loop"];
	Print["Loop timing=",N[AbsoluteTiming[
	
	Do[Module[{Himp,evec,JK,simp,s1,s2 },
		JK=KondoCouplings[[j]];
		Himp=N[  (H0+JK HI)];
		Print["    Eigenvector timing=",N[AbsoluteTiming[(*ev=Sort@Eigenvalues[Himp,2 klevels];*)
		evec=(-Eigensystem[-Himp, 1,
		Method -> {"Arnoldi","Criteria"->"RealPart","MaxIterations"->2000,"Tolerance"->10^-9}])[[2,1]];  ][[1]]/60], " min ;   j=",j,"/",Length@KondoCouplings];
		simp=Table[Conjugate[evec] . spinImpOp[Simp,2 Lx Ly ][[\[Gamma]]] . evec,{\[Gamma],1,3}]; 
		s1=Table[Conjugate[evec] . spinOp[Simp,1,2 Lx Ly ][[\[Gamma]]] . evec,{\[Gamma],1,3}];
		s2=Table[Conjugate[evec] . spinOp[Simp,2,2 Lx Ly ][[\[Gamma]]] . evec,{\[Gamma],1,3}];
		dataAppend[StringJoin[datapath,".txt"],
		Chop@{JK,{simp . avec,simp . bvec,simp . cvec},{s1 . avec,s1 . bvec,s1 . cvec},{s2 . avec,s2 . bvec,s2 . cvec}}    ];
		 
		(*Print["    Memory in use:  ",N[10^-9  MemoryInUse[] ]  ];*)
],{j,1,Length@KondoCouplings}]  ][[1]]/60], " min." ];

(*AbsoluteTiming@dataWrite[datapath,eValues];*)

]


(* ::Subsection:: *)
(*Closing Kernels*)


If[ ($FrontEnd===Null), 
	Print[];
	Print["Closing Kernels"];
	CloseKernels[];
];

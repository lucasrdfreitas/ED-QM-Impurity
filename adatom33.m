(* ::Package:: *)

(* ::Text:: *)
(*Adatom  -- varying Kondo coupling JK*)
(*	-compute eigenvalues for first "k" levels*)
(*	-compute eigenvector for the ground state*)
(*Coupling: Kitaev FM*)
(*dataName: KondoCoupling ranging from 0 to 1 with spacing .1*)


(* ::Subsection:: *)
(*preamble*)


If[ \[Not]($FrontEnd===Null), SetDirectory[NotebookDirectory[]] ];
$FileName=If[$FrontEnd === Null, $InputFileName, NotebookFileName[] ];
Get[ FileNameJoin[{Directory[],"definitions.wl" }] ]


systemDimensions = {{3,3}};
kitaev           = {{-1,-1,-1}};
heisenberg       = {0{1,1,1}};
anisotropy       = {0{1,1,1}};
hfields          = {0.4 cvec};
impuritySpin     = {1};
gs               = {1};
KondoCouplings   = {0,1,.05};
parameters       = N@Tuples[{systemDimensions,kitaev,heisenberg,anisotropy,hfields,impuritySpin,gs} ];
klevels          = 100;

HamCoupling="Kitaev_FM";
dataName=Module[{i,f,\[Delta],k}, 
				{i,f,\[Delta]}=KondoCouplings;  k=klevels;    {i,f,\[Delta],k}=ToString/@{i,f,\[Delta],k};
				StringReplace["JK=Range[i,f,d]_k=k0", {"i"->i,"f"->f,"d"->\[Delta],"k0"->k}]  ];
KondoCouplings=Range@@KondoCouplings;


(* klevels = 20 ~ 5 min per JK (Eigenvalue[]) *)


(* ::Subsection:: *)
(*Launching Kernel *)


If[ ($FrontEnd===Null),
	Print["Before Starting Kernels"];
	Needs["ClusterIntegration`"];
	(*kernels = LaunchKernels[SGE["micro4", 10]];*)
	Quiet[kernels = LaunchKernels[]];
	Print["Starting Kernels"];
];


(* ::Subsection:: *)
(*code -- save matrices*)


Print[];


Module[{Lx,Ly,J,\[Lambda]n,Simp,K,h,g,H0,HK,HJ,HZ,HI,eValues,info,path,Huncompressed}, 
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
	Print["Compressed data found - Skipping to compute eigenvalues"]
];

];


(* ::Subsection:: *)
(*Code -- eigenvalues*)


Print[];
Print["Computing Eigenvalues"];
Print[];


Module[{Lx,Ly,J,\[Lambda]n,Simp,K,h,g,H0,HJ,HI,HK,HZ,eValues,path,info,datapath},
{{Lx,Ly},K,J,\[Lambda]n,h,Simp,g}=parameters[[1]]; 
{Lx,Ly}=Round@{Lx,Ly};eValues={};
	
	datapath=dataPath[dataName,HamCoupling,Simp,{Lx,Ly},dataFolder];
	Print["Data path : ",datapath];

	info=StringReplace["simp=X_h=Y",{"X"->ToString@Simp,"Y"->ToString@N[Round[1000 Norm@h]/1000]}];
	path=dataPath[#,HamCoupling,Simp,{Lx,Ly},dataFolder]&@(StringJoin@{{"Hmatrices_"},{info}}); 
	Print["Uncompress time=", AbsoluteTiming[ {H0,HI}=dataZipImport[path]; ]];
	Print[" "];
	Print["    Memory in use:  ",N[10^-9  MemoryInUse[] ]  ];

	Print["Starting JK Loop"];
	Print["Loop timing=",AbsoluteTiming[
	
	Do[Module[{Himp,ev,JK },
		JK=KondoCouplings[[j]];
		Himp=N[  (H0+JK HI)];
		Print["    Eigenvalue timing=",AbsoluteTiming[(*ev=Sort@Eigenvalues[Himp,2 klevels];*)
		ev=Sort@(-Eigenvalues[-Himp, klevels,
		Method -> {"Arnoldi","Criteria"->"RealPart","MaxIterations"->600,"Tolerance"->10^-9(*, "Shift"->-.7*)}]);  ]];
		
		dataAppend[datapath,{JK,ev}];
		(*AppendTo[eValues,{JK,ev}];*)
		Print["    j=",j,"/",Length@KondoCouplings];
		Print["    Memory in use:  ",N[10^-9  MemoryInUse[] ]  ];
],{j,1,Length@KondoCouplings}]] ];

(*AbsoluteTiming@dataWrite[datapath,eValues];*)


]


N[5773/60]


5*21


(* ::Subsection:: *)
(*Closing Kernels*)


If[ ($FrontEnd===Null), 
	Print[];
	Print["Closing Kernels"];
	CloseKernels[];
];

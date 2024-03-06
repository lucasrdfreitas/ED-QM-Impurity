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
impuritySpin     = {1/2};
gs               = {1};
KondoCouplings   = Range[0,1,.1];
parameters       = N@Tuples[{systemDimensions,kitaev,heisenberg,anisotropy,hfields,impuritySpin,gs(*, KondoCouplings,*)} ];
klevels          = 50;

HamCoupling="Kitaev_FM";

dataName=Module[{i,f,\[Delta],k}, {i,f}=KondoCouplings[[{1,-1}]];\[Delta]=(f-i)/Length[KondoCouplings];
{i,f,\[Delta],k}=ToString/@{i,f,\[Delta],k};
StringReplace["JK=Range[i,f,d]_k=K",{"i"->i,"f"->f,"d"->\[Delta],"K"->k}]  ];


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
	Print["Memory in use:  ",     N[10^-9  MemoryInUse[] ]  ];
	Print["H Kitaev timing=",     AbsoluteTiming[HK=If[Norm[K]==0,0,N@AdatomKitaev[K,Simp,Lx,Ly]]; Dimensions@HK ] ];
	Print["H Heisenberg timing=", AbsoluteTiming[HJ=If[Norm[J]==0,0,N@AdatomHeisenberg[J,\[Lambda]n,Simp,Lx,Ly]];Dimensions@HJ] ];
	Print["H Zeeman timing=",     AbsoluteTiming[HZ=If[Norm[h]==0,0,N@AdatomZeeman[h,Simp,g,2 Lx Ly]];Dimensions@HZ] ];
	Print["H imp timing=",        AbsoluteTiming[HI=N@AdatomImp[1,Simp,2 Lx Ly]; Dimensions@HI] ];
	H0=HK+HJ+HZ;
	Clear[HK,HJ,HZ];
	Print["Memory in use:  ",     N[10^-9  MemoryInUse[] ]  ];
	Print["Compress time=",       AbsoluteTiming[ dataZipExport[path,{H0,HI}] ;]     ];
	,
	Print["Compressed data found - Skipping to compute eigenvalues"]
];(*Print["Loop timing=",AbsoluteTiming@Do[Module[{Himp,ev,JK },
	Print["Memory in use:  ",N[10^-9  MemoryInUse[] ]  ];
	JK=KondoCouplings[[j]];
	Himp=H0+JK HI;
	ev=Sort@Eigenvalues[N@ Himp,2 klevels];
	AppendTo[eValues,{JK,ev}];
],{j,1,Length@KondoCouplings}] ];

Print["Write timing=",
AbsoluteTiming@dataWrite[dataPath[dataName,HamCoupling,Simp,{Lx,Ly},dataFolder],eValues]]; *)
];


(* ::Subsection:: *)
(*Code -- eigenvalues*)


Print[];
Print["Computing Eigenvalues"];
Print[];


Module[{Lx,Ly,J,\[Lambda]n,Simp,K,h,g,H0,HJ,HI,HK,HZ,eValues,path,info},
{{Lx,Ly},K,J,\[Lambda]n,h,Simp,g}=parameters[[1]]; 
{Lx,Ly}=Round@{Lx,Ly};eValues={};

	info=StringReplace["simp=X_h=Y",{"X"->ToString@Simp,"Y"->ToString@N[Round[1000 Norm@h]/1000]}];
	path=dataPath[#,HamCoupling,Simp,{Lx,Ly},dataFolder]&@(StringJoin@{{"H0_"},{info}}); 
	Print["Uncompress time=", AbsoluteTiming[ {H0,HI}=dataZipImport[path]; ]];
	Print[" "];
	Print["Memory in use:  ",N[10^-9  MemoryInUse[] ]  ];


Print["Loop timing=",AbsoluteTiming@Do[Module[{Himp,ev,JK },
	Print["Memory in use:  ",N[10^-9  MemoryInUse[] ]  ];
	JK=KondoCouplings[[j]];
	Himp=H0+JK HI;
	ev=Sort@Eigenvalues[N@ Himp,2 klevels];
	AppendTo[eValues,{JK,ev}];
],{j,1,Length@KondoCouplings}] ];

Print["Write timing=",
AbsoluteTiming@dataWrite[dataPath[dataName,HamCoupling,Simp,{Lx,Ly},dataFolder],eValues]];

]


(* ::Subsection:: *)
(*Closing Kernels*)


If[ ($FrontEnd===Null), 
	Print[];
	Print["Closing Kernels"];
	CloseKernels[];
];

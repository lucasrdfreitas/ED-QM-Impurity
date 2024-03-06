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


(*Print["Before Starting Kernels"];
Needs["ClusterIntegration`"];
(*kernels = LaunchKernels[SGE["micro4", 10]];*)
Quiet[kernels = LaunchKernels[]];
Print["Starting Kernels"];*)


(* ::Subsection:: *)
(*code*)


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


(*Print[];
Print["Closing Kernels"];
CloseKernels[];*)

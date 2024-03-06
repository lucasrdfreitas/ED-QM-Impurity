(* ::Package:: *)

(* ::Text:: *)
(*Adatom  -- varying Kondo coupling JK*)
(*	-compute eigenvalues for first "k" levels*)
(*	-compute eigenvector for the ground state*)
(*Coupling: Kitaev FM  *)
(*dataName: KondoCoupling ranging from JKi to JKf with spacing \[Delta]J*)


HamCoupling="Kitaev_FM";
dataName=Module[{i=0,f=1,\[Delta]=.01},{i,f,\[Delta]}=ToString/@{i,f,\[Delta]};
StringReplace["JK_Range[i,f,d]",{"i"->i,"f"->f,"d"->\[Delta]}]];


(* ::Subsection:: *)
(*preamble*)


If[ $FrontEnd != Null, SetDirectory[NotebookDirectory[]] ];

Get[ FileNameJoin[{Directory[],"definitions.wl" }] ]


systemDimensions = {{2,2}};
kitaev           = {{-1,-1,-1}};
heisenberg       = {0{1,1,1}};
anisotropy       = {0{1,1,1}};
hfields          = {0.4 cvec};
impuritySpin     = {1/2};
gs               = {1};
KondoCouplings   = Range[0,1,.1];
parameters       = N@Tuples[{systemDimensions,kitaev,heisenberg,anisotropy,hfields,impuritySpin,gs(*, KondoCouplings,*)} ];
klevels          = 50;


(* ::Subsection:: *)
(*code*)


Module[{Lx,Ly,J,\[Lambda]n,Simp,K,h,g,H0,HJ,HI,HK,HZ,eValues},
{{Lx,Ly},K,J,\[Lambda]n,h,Simp,g}=parameters[[1]]; 
{Lx,Ly}=Round@{Lx,Ly};

Print["H Kitaev timing=",AbsoluteTiming[  
	HK=If[Norm[K]==0,0,N@AdatomKitaev[K,Simp,Lx,Ly]]; 
	Dimensions@HK
] ];
Print["H Heisenberg timing=", AbsoluteTiming[  
	HJ=If[Norm[J]==0,0,N@AdatomHeisenberg[J,\[Lambda]n,Simp,Lx,Ly]];
	Dimensions@HJ
] ];
Print["H Zeeman timing=", AbsoluteTiming[  
	HZ=If[Norm[h]==0,0,N@AdatomZeeman[h,Simp,g,2 Lx Ly]];
	Dimensions@HZ
] ];
Print["H imp timing=", AbsoluteTiming[  
	HI=N@AdatomImp[1,Simp,2 Lx Ly];
	Dimensions@HI]
];

Print["H0 sum timing=",AbsoluteTiming[H0=HK+HJ+HZ;
	Dimensions@H0] ]; 

Print[" "];

eValues={};
Print["Loop timing=",AbsoluteTiming@Do[Module[{Himp,ev,JK },
	JK=KondoCouplings[[j]];
	Himp=H0+JK HI;
	ev=Sort@Eigenvalues[N@ Himp,2klevels];
	AppendTo[eValues,{JK,ev[[1;;klevels]]-ev[[1]]}];
],{j,1,Length@KondoCouplings}] ];

Print["Write timing=",
AbsoluteTiming@dataWrite[dataName,HamCoupling,Simp,{Lx,Ly},eValues]];


];




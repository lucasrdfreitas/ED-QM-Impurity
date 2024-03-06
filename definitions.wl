(* ::Package:: *)

(* ::Title:: *)
(*Quantum Magnets + Impurity*)


(* ::Text:: *)
(*This script employs exact diagonalization to compute the magnon spectrum and spin expectation values in various Quantum Magnets featuring a magnetic impurity (adatom) coupled to a single site. These materials are effectively modeled by the JK\[CapitalGamma] extended Kitaev model.*)
(*	> we will focus mainly in two materials:  Subscript[RuCl, 3]  (e.g. 1706.06113 ) and Subscript[CrI, 3] (e.g. 1704.03849 )*)


(* ::Subsection:: *)
(*System files*)


(* ::Text:: *)
(*create data folder *)


If[ $FrontEnd != Null, SetDirectory[NotebookDirectory[]] ];
Print[Directory[] ];

NbPath =  If[$FrontEnd === Null, $InputFileName, NotebookFileName[] ];
Print[NbPath];

NbName = FileBaseName[NbPath];      (*StringSplit[ FileNameSplit[ NbPath ][[-1]]  ,"."][[1]];*)
Print[NbName];

NbPath =  If[$FrontEnd === Null, $FileName, NotebookFileName[] ];
Print[NbPath];

NbName = FileBaseName[NbPath];      
Print[NbName];


Abort[];

dataFolder=FileNameJoin[{Directory[],"data",NbName }];

If[Length@FileNames[dataFolder]==0,CreateDirectory@File@dataFolder];


(* ::Text:: *)
(*functions to load and write data*)


createDirectory[path_] :=Module[ {l=Length@FileNames[path]},
	If[ l==0, CreateDirectory@File@path; ,Null] ];


dataPath[dataName_,Ham_,Simp_,L_]:= Module[{folderName,folderPath,dataPath},
folderName=StringJoin[Ham,"_size=(",ToString@L[[1]],",",ToString@L[[2]],")_2Simp=",ToString@Round[2 Simp] ];
folderPath=FileNameJoin[{ dataFolder,folderName }];
dataPath=FileNameJoin[{ folderPath, dataName }];   
{folderPath,dataPath}
];


dataWrite[dataName_,Ham_,Simp_,L_,data_] :=
Module[ {folderPath,dataPath0,auxStream},	
	{folderPath,dataPath0}=dataPath[dataName,Ham,Simp,L];
	createDirectory@folderPath;
	auxStream = OpenWrite[dataPath0];
	Write[auxStream, data];
	Close[auxStream];                ];


dataAppend[dataName_,Ham_,Simp_,L_,data_] :=
Module[ {folderPath,dataPath0,auxStream},	
	{folderPath,dataPath0}=dataPath[dataName,Ham,Simp,L];
	createDirectory@folderPath;
	auxStream = OpenAppend[dataPath0];
	Write[auxStream, data];
	Close[auxStream];                ];


dataRead[datapath_]:=
Module[ {auxStream,data},
	auxStream = OpenRead[datapath];
	If[auxStream==$Failed, Print["Failed to OpenRead file at: "]; 
	Print[ datapath ]; Abort[] ];
	data=ReadList[auxStream];
	Close[auxStream];			
	data[[-1]]
];


(* ::Subsection::Closed:: *)
(*Basic definitions *)


(* ::Text:: *)
(*Base vectors*)


avec=1/Sqrt[6] {1,1,-2};
bvec=1/Sqrt[2] {-1,1,0};
cvec=1/Sqrt[3] {1,1,1};


(* ::Text:: *)
(*For impurity spin S we define spin components *)


spinmatrix[S_]:=Module[{n=Round[2S+1],S0,Sx,Sy,Sz}, If[FractionalPart[2 S]!=0.,Print["2(",S,")+1  is not an integer."];Abort[];];
S0=IdentityMatrix[n];
Sx=1/2 SparseArray[{ {i_,j_}/;i-j==1 :>Sqrt[(S+1)(i+j-1)-i j], {i_,j_}/;i-j==-1 :>Sqrt[(S+1)(i+j-1)-i j]},{n,n},0];
Sy=I/2 SparseArray[{ {i_,j_}/;i-j==1 :>Sqrt[(S+1)(i+j-1)-i j], {i_,j_}/;i-j==-1 :>-Sqrt[(S+1)(i+j-1)-i j] },{n,n},0];
Sz=SparseArray[{ {i_,i_}:>(S+1-i)},{n,n},0];
SparseArray/@N@{S0,Sx,Sy,Sz} ];


(* ::Subsection:: *)
(*Adatom Hamiltonian *)


(* ::Text:: *)
(*	Kitaev FM + AFM Kondo adatom impurity with N0=2LxLy bulk spins  *)


AdatomHamiltonian[J_,\[Lambda]n_,K_,h_,JK_,Simp_,g_,Lx_,Ly_]:=
AdatomKitaev[K,Simp,Lx,Ly]+AdatomHeisenberg[J,\[Lambda]n,Simp,Lx,Ly]+AdatomZeeman[h,Simp,g,2Lx Ly]+AdatomImp[JK,Simp,2Lx Ly];


(* ::Text:: *)
(*Note on optimization - if one want to vary some parameter (e.g. JK, h) don't use this function,  rather*)
(*first compute each Hamiltonian matrix and then sum them. To avoid compute the same matrix Hamiltonian over and over.*)


(* ::Text:: *)
(*Couplings -- K<0 = FM, K>0 AFM  (similarly for J, JK)*)
(*J=(Jx,Jy,Jz) Isotropic NN Heisenberg *)
(*\[Lambda]n=\[Lambda] (nx,ny,nz)  anisotropy \[Lambda] of Heisenberg coupling along the n direction*)
(* K=(Kx,Ky,Kz) Kitaev coupling  *)


(* ::Subsubsection:: *)
(*Zeeman Hamiltonian*)


(* ::Text:: *)
(*bulk: position 1, ..., 8  with spin 1/2,  and impurity spin Subscript[S, imp] at top of position 8 which is (chosen as) the 9th slot in the tensor product*)
(*		> spin/bond index \[Alpha] = x,y,z = 1,2,3 *)
(*		> bulk spins  bulk[[i,\[Alpha]]] for positions i=1, ..., 8 *)
(*			PadRight[ Subscript[s, \[Alpha]], 8, Subscript[s, 0], i-1 ] creates a list with spin Subscript[s, \[Alpha]] at position i and identity Subscript[s, 0] elsewhere*)
(*			  Use KroneckerProduct to perform the tensor product of the 8+1 spin matrices*)
(*		> impurity spin imp[[\[Alpha]]]*)
(* [Most of the computational time goes in the KroneckerProduct. Using SparseArrays (i.e. spinmatrix) makes the computation 2 or 3 orders of magnitude faster.]*)


AdatomZeeman[h_,Simp_,g_,N0_]:=Module[{s,S,bulk,imp}, 
s=spinmatrix[1/2]; 
S=spinmatrix[Simp];
bulk = Table[  KroneckerProduct@@Join[PadRight[  {s[[\[Alpha]+1]]}, N0, {s[[1]]}, i-1], {S[[1]]}     ]  ,{i,1,N0},{\[Alpha],1,3}];
imp  = Table[  KroneckerProduct@@Join[Table[     s[[1]],   N0],              {S[[\[Alpha]+1]]} ]  ,{\[Alpha],1,3}];

Sum[ -h[[\[Alpha]]] g imp[[\[Alpha]]] - h[[\[Alpha]]] Sum[bulk[[i,\[Alpha]]],{i,1,N0}]    ,{\[Alpha],1,3}]
];


(* ::Subsubsection::Closed:: *)
(*Impurity Hamiltonian*)


(* ::Text:: *)
(*Adatom impurity couples the 8th bulk spin to the impurity spin (the 9th slot in the tensor product) with AFM Heisenberg with coupling Subscript[J, K]*)


AdatomImp[JK_,Simp_,N0_]:=Module[{s,S,sS}, 
s=spinmatrix[1/2]; 
S=spinmatrix[Simp];

sS=Table[  KroneckerProduct@@Join[PadRight[ {s[[\[Alpha]+1]]} , N0, {s[[1]]}, N0-1] ,{ S[[\[Alpha]+1]] } ]  ,{\[Alpha],1,3}]; 
JK Sum[sS[[\[Alpha]]],{\[Alpha],1,3}]
]


(* measuring max memory allocated/ note that ByteCount overstimate the size (uncompressed value) *)
(*Do[Module[{hz}, Print@{N0+1,AbsoluteTiming[  hz=HZeeman[cvec,1/2,1,N0]; Dimensions@hz ], N[10^-9 ByteCount[hz]]  }  ] ,{N0,4,15}]*)


(*Do[Module[{himp}, Print@{N0+1,AbsoluteTiming[  himp=Himp[1,1/2,N0]; Dimensions@himp ], N[10^-9 ByteCount[himp]]  }  ] ,{N0,4,20}]*)


(* ::Subsubsection::Closed:: *)
(*Kitaev Hamiltonian*)


(* ::Text:: *)
(*Bulk Kitaev Hamiltonian  - *)
(*	> for the Honeycomb lattice consider Lx times Ly triangular lattice with N0=2LxLy bulk spins*)
(*	> unit cell index  r=m + n Lx+1 for m =0, ...,Lx-1  and n =0, ..., Ly-1 *)
(*	> position i = 2r-1 + \[Nu]   for sublattice label \[Nu]=0,1=A,B*)
(*The Kitaev interactions happens between the positions  bonds[[\[Alpha],r]]  *)
(*	> x-bonds:    (2r-1,2rx);         rx=Mod[ m+1,Lx]+ n Lx+1                                    e.g. for Lx=Ly=2:  (1,4), (3,2), (5,8), (7,6)*)
(*	> y-bonds:    (2r-1,2ry);         ry=m+Mod[ n+1,Ly] Lx+1                                                                        (1,6), (3,8), (5,2), (7,4)*)
(*	> z-bonds:    (2r-1,2r );                                                                                                                                      (1,2), (3,4), (5,6), (7,8)*)


Bonds[Lx_,Ly_]:=Transpose@Flatten[#,1]&@Table[ 
{   {2 (m+n Lx+1)-1, 2 (Mod[m+1,Lx]+n Lx+1) },{2 (m+n Lx+1)-1, 2 (m+Mod[n+1,Ly] Lx+1) },{2 (m+n Lx+1)-1, 2 (m+n Lx+1)}  }
,{n,0,Ly-1},{m,0,Lx-1}];


(* ::Text:: *)
(*      For each bond (i,j) we use Insert at the positions Min[i,j], and Max[i,j]  the corresponding spin s^\[Alpha] and s^ \[Beta]. *)
(*[Note that for Kitaev Hamiltonian it does not matter if i>j or i<j since in the Hamiltonian only appears Subscript[s, i]^\[Alpha] Subscript[s, j]^\[Beta] for \[Alpha]=\[Beta].*)
(*However, we must be careful with the ordering in the Kitaev-\[CapitalGamma] model ]*)
(*Let us denote by Subscript[H, \[Alpha]]  the terms in the Kitaev Hamiltonian for \[Alpha]-bonds.*)


AdatomKitaev[K_,Simp_,Lx_,Ly_]:=Module[{s,S,bonds,Hx,Hy,Hz,N0},
N0=2 Lx Ly;
s=spinmatrix[1/2]; 
S=spinmatrix[Simp];
bonds=Bonds[Lx,Ly];

{Hx,Hy,Hz} = Table[  K[[\[Alpha]]] Sum[ 
KroneckerProduct@@Join[ Insert[s[[\[Alpha]+1]],Max@bonds[[\[Alpha],r]] ]@Insert[s[[\[Alpha]+1]],Min@bonds[[\[Alpha],r]]]@Table[s[[1]],N0-2], {S[[1]]}  ]  
,{r,1,Lx Ly}],{\[Alpha],1,3}];

Hx+Hy+Hz
]


(*AbsoluteTiming@HKitaev[{1,1,1},1/2,3,3]*)


(* ::Subsubsection::Closed:: *)
(*Heisenberg Hamiltonian*)


(* ::Text:: *)
(*Bulk anisotropic Heisenberg Hamiltonian  - *)
(*	> follows the same structure as the Kitaev Hamiltonian, but each Subscript[H, \[Alpha]] has a additional sum over spin indices *)
(*	> we will also consider an anisotropy of strength \[Lambda] along the n axis *)


AdatomHeisenberg[J_,\[Lambda]n_,Simp_,Lx_,Ly_]:=Module[{s,S,bonds,HJ,H\[Lambda],N0,n,\[Lambda],sn},
N0=2 Lx Ly;
s=spinmatrix[1/2]; 
S=spinmatrix[Simp];
\[Lambda]=Norm[\[Lambda]n];
n=\[Lambda]n/\[Lambda];
sn=s[[2;;4]] . n;
bonds=Bonds[Lx,Ly];

HJ =Sum[J[[\[Alpha]]]KroneckerProduct@@Join[Insert[s[[\[Beta]+1]],Max@bonds[[\[Alpha],r]][[2]]]@Insert[s[[\[Beta]+1]],Min@bonds[[\[Alpha],r]]]@Table[s[[1]],N0-2], {S[[1]]} ]  
,{r,1,Lx Ly},{\[Beta],1,3},{\[Alpha],1,3}];
H\[Lambda]=Sum[\[Lambda] KroneckerProduct@@Join[Insert[sn,Max@bonds[[\[Alpha],r]][[2]]]@Insert[sn,Min@bonds[[\[Alpha],r]]]@Table[s[[1]],N0-2], {S[[1]]} ]  
,{r,1,Lx Ly},{\[Alpha],1,3}];

HJ+H\[Lambda]]

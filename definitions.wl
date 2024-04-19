(* ::Package:: *)

(* ::Title:: *)
(*Quantum Magnets + Impurity*)


(* ::Text:: *)
(*This script employs exact diagonalization to compute the magnon spectrum and spin expectation values in various Quantum Magnets featuring a magnetic impurity (adatom) coupled to a single site. These materials are effectively modeled by the JK\[CapitalGamma] extended Kitaev model.*)
(*	> we will focus mainly in two materials:  Subscript[RuCl, 3]  (e.g. 1706.06113 ) and Subscript[CrI, 3] (e.g. 1704.03849 )*)


(* ::Subsection::Closed:: *)
(*System files*)


(* ::Text:: *)
(*create data folder   (the variable $FileName should be defined externally when executing this file)*)


If[ $FrontEnd != Null, SetDirectory[NotebookDirectory[]] ];
NbPath =  If[$FrontEnd === Null, $FileName, NotebookFileName[] ];
NbName = FileBaseName[NbPath];      
Print[NbName];
dataFolder=FileNameJoin[{Directory[],"data",NbName }];
If[Length@FileNames[dataFolder]==0,CreateDirectory@File@dataFolder];


(* ::Text:: *)
(*functions to load and write data*)


createDirectory[path_] :=Module[ {l=Length@FileNames[path]},
	If[ l==0, CreateDirectory@File@path; ,Null] ];


dataPath[dataName_,Ham_,Simp_,L_,dataFolder_:dataFolder]:= Module[{folderName,folderPath,dataPath},
folderName=StringJoin[Ham,"_size=(",ToString@L[[1]],",",ToString@L[[2]],")_Simp=",ToString@NumberForm[N[Simp],{\[Infinity],1}]   ];
folderPath=FileNameJoin[{ dataFolder,folderName }];
dataPath=FileNameJoin[{ folderPath, dataName }];   
dataPath
];


dataPathTXT[dataName_,Ham_,Simp_,L_,dataFolder_:dataFolder]:= Module[{folderName,folderPath,dataPath,dataNametxt},
folderName=StringJoin[Ham,"_size=(",ToString@L[[1]],",",ToString@L[[2]],")_Simp=",ToString@NumberForm[N[Simp],{\[Infinity],1}]   ];
folderPath=FileNameJoin[{ dataFolder,folderName }];
dataNametxt=If[FileExtension["file"]!="txt",StringJoin[dataName,".txt"]  ];
dataPath=FileNameJoin[{ folderPath, dataNametxt }];   
dataPath
];


dataZipExport[datapath_,data_] :=Module[ {dirpath,auxStream,path},	
	dirpath=DirectoryName[datapath];
	createDirectory@dirpath;
	auxStream = OpenWrite[datapath];
	Write[auxStream, data];
	Close[auxStream];  
	CreateArchive[datapath,StringJoin[datapath,".zip"],
	OverwriteTarget->True,"CompressionLevel"->1];
	(* delete uncompressed data *) 
	DeleteFile[datapath];
	             ];


dataZipImport[datapath_] := Module[ {auxStream,data,dirpath},
	dirpath=DirectoryName[datapath];
	ExtractArchive[StringJoin[datapath,".zip"],dirpath,
	OverwriteTarget->True];
	
	auxStream = OpenRead[datapath];
	If[auxStream==$Failed, Print["Failed to OpenRead file at: "]; 
	Print[ datapath ]; Abort[] ];
	data=ReadList[auxStream];
	Close[auxStream];
	(* delete uncompressed data *)
	DeleteFile[datapath];
	data[[-1]]
	];


dataWrite[datapath_,data_] :=
Module[ {auxStream},	
	createDirectory@DirectoryName[datapath];
	auxStream = OpenWrite[datapath];
	Write[auxStream, data];
	Close[auxStream];                ];


dataAppend[datapath_,data_] :=
Module[ {auxStream},	
	createDirectory@DirectoryName[datapath];
	auxStream = OpenAppend[datapath];
	Write[auxStream, data];
	Close[auxStream];                ];


dataRead[datapath_]:=
Module[ {auxStream,data},
	auxStream = OpenRead[datapath];
	If[auxStream==$Failed, Print["Failed to OpenRead file at: "]; 
	Print[ datapath ]; Abort[] ];
	data=ReadList[auxStream];
	Close[auxStream];			
	data
];


(* ::Subsection::Closed:: *)
(*Basic definitions *)


round[x_,f_:10^4]:=N[Round[x f]/f]; 


(* ::Text:: *)
(*Base vectors*)


avec=1/Sqrt[6] {1,1,-2};
bvec=1/Sqrt[2] {-1,1,0};
cvec=1/Sqrt[3] {1,1,1};


(* ::Text:: *)
(*For impurity spin S we define spin components  *)


spinmatrix[S_]:=Module[{n=Round[2S+1],S0,Sx,Sy,Sz}, If[FractionalPart[2 S]!=0.,Print["2(",S,")+1  is not an integer."];Abort[];];
S0=IdentityMatrix[n];
Sx=1/2 SparseArray[{ {i_,j_}/;i-j==1 :>Sqrt[(S+1)(i+j-1)-i j], {i_,j_}/;i-j==-1 :>Sqrt[(S+1)(i+j-1)-i j]},{n,n},0];
Sy=I/2 SparseArray[{ {i_,j_}/;i-j==1 :>Sqrt[(S+1)(i+j-1)-i j], {i_,j_}/;i-j==-1 :>-Sqrt[(S+1)(i+j-1)-i j] },{n,n},0];
Sz=SparseArray[{ {i_,i_}:>(S+1-i)},{n,n},0];
SparseArray/@N@{Sx,Sy,Sz,S0} ];


(* ::Text:: *)
(*spin operators in Cartesian basis -- *)
(*	Adatom          N0 =  2 Lx Ly*)
(*	Substitution N0 =  2 Lx Ly -1*)


spinOp[Simp_,i_,N0_,Sbulk_:1/2]:= 
Module[{s,S},
	s=spinmatrix[Sbulk]; 
	S=spinmatrix[Simp];
	Table[  KroneckerProduct@@Join[{S[[4]]}, Insert[s[[\[Alpha]]],i]@Table[s[[4]],N0-1]  ]  ,{\[Alpha],1,3}]   ];
	
spinImpOp[Simp_,N0_,Sbulk_:1/2]:= 
Module[{s,S},
	s=spinmatrix[Sbulk]; 
	S=spinmatrix[Simp];
	Table[  KroneckerProduct@@Join[ {S[[\[Alpha]]]}, Table[s[[4]],N0]  ]  ,{\[Alpha],1,3}]   ];


(* ::Subsection::Closed:: *)
(*Adatom Hamiltonian *)


(* ::Text:: *)
(*	Kitaev FM + AFM Kondo adatom impurity with N0=2LxLy bulk spins  *)


AdatomHamiltonian[J_,\[Lambda]n_,K_,h_,JK_,Simp_,g_,Lx_,Ly_,Sbulk_:1/2]:=
AdatomKitaev[K,Simp,Lx,Ly,Sbulk]+AdatomHeisenberg[J,\[Lambda]n,Simp,Lx,Ly,Sbulk]+AdatomZeeman[h,Simp,g,2Lx Ly,Sbulk]+AdatomImp[JK,Simp,2Lx Ly,Sbulk];


(* ::Text:: *)
(*Note on optimization - if one want to vary some parameter (e.g. JK, h) don't use this function,  rather*)
(*first compute each Hamiltonian matrix and then sum them. To avoid compute the same matrix Hamiltonian over and over.*)


(* ::Text:: *)
(*Couplings -- K<0 = FM, K>0 AFM  (similarly for J, JK)*)
(*J=(Jx,Jy,Jz) Isotropic NN Heisenberg *)
(*\[Lambda]n=\[Lambda] (nx,ny,nz)  anisotropy \[Lambda] of Heisenberg coupling along the n direction*)
(* K=(Kx,Ky,Kz) Kitaev coupling  *)


(* ::Subsubsection::Closed:: *)
(*Zeeman Hamiltonian*)


(* ::Text:: *)
(*bulk: position 1, ..., N  with spin 1/2,  and impurity spin Subscript[S, imp] at top of position 1 which is (chosen as) the first slot in the tensor product*)
(*		> spin/bond index \[Alpha] = x,y,z = 1,2,3 *)
(*		> bulk spins  bulk[[i,\[Alpha]]] for positions i=1, ..., 8 *)
(*			Insert[ Subscript[s, \[Alpha]], i ] @Table[s[[4]],N0-1] creates a list with spin Subscript[s, \[Alpha]] at position i and identity Subscript[s, 0] elsewhere*)
(*			 Use KroneckerProduct to perform the tensor product of the N+1 spin matrices*)
(*		> impurity spin imp[[\[Alpha]]]*)
(* [Most of the computational time goes in the KroneckerProduct. Using SparseArrays (i.e. spinmatrix) makes the computation 2 or 3 orders of magnitude faster.]*)


AdatomZeeman[h_,Simp_,g_,N0_,Sbulk_:1/2]:=Module[{s,S,bulk,imp}, 
s=spinmatrix[Sbulk]; 
S=spinmatrix[Simp];
bulk = Table[  KroneckerProduct@@Join[ {S[[4]]}, Insert[s[[\[Alpha]]],i]@Table[s[[4]],N0-1]    ]  ,{i,1,N0},{\[Alpha],1,3}];
imp  = Table[  KroneckerProduct@@Join[ {S[[\[Alpha]]]}, Table[    s[[4]],   N0]               ]  ,{\[Alpha],1,3}];

Sum[ -h[[\[Alpha]]] g imp[[\[Alpha]]] - h[[\[Alpha]]] Sum[bulk[[i,\[Alpha]]],{i,1,N0}]    ,{\[Alpha],1,3}]
];


(* ::Subsubsection::Closed:: *)
(*Impurity Hamiltonian*)


(* ::Text:: *)
(*Adatom impurity couples the 8th bulk spin to the impurity spin (the 9th slot in the tensor product) with AFM Heisenberg with coupling Subscript[J, K]*)


AdatomImp[JK_,Simp_,N0_,Sbulk_:1/2]:=Module[{s,S,sS}, 
s=spinmatrix[Sbulk]; 
S=spinmatrix[Simp];

sS=Table[  KroneckerProduct@@Join[ {S[[\[Alpha]]]}, Insert[s[[\[Alpha]]],1]@Table[s[[4]],N0-1]  ]  ,{\[Alpha],1,3}]; 

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


AdatomKitaev[K_,Simp_,Lx_,Ly_,Sbulk_:1/2]:=Module[{s,S,bonds,Hx,Hy,Hz,N0},
N0=2 Lx Ly;
s=spinmatrix[Sbulk]; 
S=spinmatrix[Simp];
bonds=Bonds[Lx,Ly];

{Hx,Hy,Hz} = Table[  K[[\[Alpha]]] Sum[ 
KroneckerProduct@@Join[{S[[4]]},  Insert[s[[\[Alpha]]],Max@bonds[[\[Alpha],r]] ]@Insert[s[[\[Alpha]]],Min@bonds[[\[Alpha],r]]]@Table[s[[4]],N0-2]    ]  
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


AdatomHeisenberg[J_,\[Lambda]n_,Simp_,Lx_,Ly_,Sbulk_:1/2]:=Module[{s,S,bonds,HJ,H\[Lambda],N0,n,\[Lambda],sn},
N0=2 Lx Ly;
s=spinmatrix[Sbulk]; 
S=spinmatrix[Simp];
\[Lambda]=-Norm[\[Lambda]n];   (*negative for FM*)
n=\[Lambda]n/\[Lambda];
sn=n[[1]]s[[1]]+n[[2]]s[[2]]+n[[3]]s[[3]];
bonds=Bonds[Lx,Ly];

HJ = Sum[J[[\[Alpha]]]KroneckerProduct@@Join[{S[[4]]}, Insert[s[[\[Beta]]],Max@bonds[[\[Alpha],r]] ]@Insert[s[[\[Beta]]],Min@bonds[[\[Alpha],r]]]@Table[s[[4]],N0-2]  ]  
,{r,1,Lx Ly},{\[Beta],1,3},{\[Alpha],1,3}];
H\[Lambda] = Sum[\[Lambda] KroneckerProduct@@Join[   {S[[4]]}, Insert[sn,  Max@bonds[[\[Alpha],r]]  ]@Insert[sn,   Min@bonds[[\[Alpha],r]]]@Table[s[[4]],N0-2]  ]  
,{r,1,Lx Ly},{\[Alpha],1,3}];

HJ+H\[Lambda]]


(* ::Subsection:: *)
(*Substitution Hamiltonian *)


(* ::Text:: *)
(*	Kitaev FM + AFM Kondo substitutional impurity with N=2LxLy bulk spins  *)


substitutionHamiltonian[J_,\[Lambda]n_,K_,h_,JK_,Simp_,g_,Lx_,Ly_,Sbulk_:1/2]:=
substitutionKitaev[K,Simp,Lx,Ly]+substitutionHeisenberg[J,\[Lambda]n,Simp,Lx,Ly]+substitutionZeeman[h,Simp,g,2Lx Ly]+substitutionImp[JK,Simp,2Lx Ly];


(* ::Text:: *)
(*put the impurity where the bulk spin i=1 would  be*)


(* ::Subsubsection::Closed:: *)
(*Zeeman Hamiltonian*)


(* ::Text:: *)
(*	let N0=N-1 be the number of bulk spins *)


substitutionZeeman[h_,Simp_,g_,N1_,Sbulk_:1/2]:=Module[{s,S,bulk,imp,N0}, 
s=spinmatrix[Sbulk]; 
S=spinmatrix[Simp];
N0=N1-1;

bulk = Table[  KroneckerProduct@@Join[ {S[[4]]}, Insert[s[[\[Alpha]]],i]@Table[s[[4]],N0-1]    ]  ,{i,1,N0},{\[Alpha],1,3}];
imp  = Table[  KroneckerProduct@@Join[ {S[[\[Alpha]]]}, Table[ s[[4]],   N0]                  ]  ,{\[Alpha],1,3}]; 

Sum[ -h[[\[Alpha]]] g imp[[\[Alpha]]] - h[[\[Alpha]]] Sum[bulk[[i,\[Alpha]]],{i,1,N0}]    ,{\[Alpha],1,3}]
];


(* ::Subsubsection:: *)
(*Kitaev Hamiltonian*)


(* ::Text:: *)
(*Bulk Kitaev Hamiltonian  - *)


substitutionBonds[Lx_,Ly_]:=Module[{bulkBonds,impBonds,bonds},
bonds=Transpose@Flatten[#,1]&@Table[ 
{   {2 (m+n Lx+1)-1, 2 (Mod[m+1,Lx]+n Lx+1) },{2 (m+n Lx+1)-1, 2 (m+Mod[n+1,Ly] Lx+1) },{2 (m+n Lx+1)-1, 2 (m+n Lx+1)}  }
,{n,0,Ly-1},{m,0,Lx-1}];
impBonds= Table[ Select[ bonds[[\[Alpha]]] , (MemberQ[#,1])& ]  ,{\[Alpha],1,3}];
bulkBonds= Table[ Complement[  bonds[[\[Alpha]]] , impBonds[[\[Alpha]]]]  ,{\[Alpha],1,3}];
{impBonds,bulkBonds}
];


(* ::Text:: *)
(*>  we need to subtract 1 in bulkBonds=substitutionBonds[Lx,Ly][[2]]-1;  since the position will be shifted be +1 when add the impurity*)
(*	e.g. the position (0,0,B) has index r=2, but it is the first position in the bulk (i.e. before inserting the impurity) *)


substitutionKitaev[K_,Simp_,Lx_,Ly_,Sbulk_:1/2]:=Module[{s,S,bulkBonds,Hx,Hy,Hz,N0},
N0=2 Lx Ly-1;
s=spinmatrix[Sbulk]; 
S=spinmatrix[Simp];
bulkBonds=substitutionBonds[Lx,Ly][[2]]-1;

{Hx,Hy,Hz} = Table[  K[[\[Alpha]]] Sum[ 
KroneckerProduct@@Join[{S[[4]]},  Insert[s[[\[Alpha]]],Max@bulkBonds[[\[Alpha],r]] ]@Insert[s[[\[Alpha]]],Min@bulkBonds[[\[Alpha],r]]]@Table[s[[4]],N0-2]    ]  
,{r,1,Length@bulkBonds}],{\[Alpha],1,3}];

Hx+Hy+Hz
]


(* ::Subsubsection:: *)
(*Impurity Hamiltonian*)


(* ::Text:: *)
(*The impurity Hamiltonian for the substitutional is obtained by summing  Subscript[s, i]^\[Beta] Subscript[s, j]^\[Beta] for ij in an \[Alpha]-bond,*)
(*the positions for the bonds are in ImpBonds, where the impurity is at r=1 and each (one of the three) bulk spins are at Max@impBonds[[\[Alpha]]]*)
(**)
(*> impBonds has three entries, corresponding to x-, y-, and z- bonds,*)
(*> impBonds[[\[Alpha]]] = { {1,r\[Alpha]} } gives the position of the impurity at 1 and the bulk spin connected to the impurity in an \[Alpha]-bond, and*)
(*> Max@impBonds[[\[Alpha]]] gives r\[Alpha], which is the index of bulk spin (shifted by -1)*)


substitutionImp[JK_,Simp_,Lx_,Ly_,Sbulk_:1/2]:=Module[{s,S,sS,impBonds,N0}, 
N0=2 Lx Ly-1;
s=spinmatrix[Sbulk]; 
S=spinmatrix[Simp];
impBonds=substitutionBonds[Lx,Ly][[1]]-1; 

sS=Sum[  KroneckerProduct@@Join[ {S[[\[Beta]]]}, Insert[s[[\[Beta]]],Max@impBonds[[\[Alpha]]]]@Table[s[[4]],N0-1]  ]  ,{\[Alpha],1,3},{\[Beta],1,3}]; 
JK sS

]


(* measuring max memory allocated/ note that ByteCount overstimate the size (uncompressed value) *)
(*Do[Module[{hz}, Print@{N0+1,AbsoluteTiming[  hz=HZeeman[cvec,1/2,1,N0]; Dimensions@hz ], N[10^-9 ByteCount[hz]]  }  ] ,{N0,4,15}]*)


(*Do[Module[{himp}, Print@{N0+1,AbsoluteTiming[  himp=Himp[1,1/2,N0]; Dimensions@himp ], N[10^-9 ByteCount[himp]]  }  ] ,{N0,4,20}]*)


(* ::Subsubsection:: *)
(*Heisenberg Hamiltonian*)


substitutionHeisenberg[J_,\[Lambda]n_,Simp_,Lx_,Ly_,Sbulk_:1/2]:=Module[{s,S,bulkBonds,HJ,H\[Lambda],N0,n,\[Lambda],sn},
N0=2 Lx Ly-1;
s=spinmatrix[Sbulk]; 
S=spinmatrix[Simp];
bulkBonds=substitutionBonds[Lx,Ly][[2]]-1;
\[Lambda]=-Norm[\[Lambda]n];   (*negative for FM*)
n=\[Lambda]n/\[Lambda];
sn=n[[1]]s[[1]]+n[[2]]s[[2]]+n[[3]]s[[3]];

HJ =Sum[J[[\[Alpha]]]KroneckerProduct@@Join[{S[[4]]},Insert[s[[\[Beta]]],Max@bulkBonds[[\[Alpha],r]]]@Insert[s[[\[Beta]]],Min@bulkBonds[[\[Alpha],r]]]@Table[s[[4]],N0-2] ]  
,{r,1,Length@bulkBonds},{\[Beta],1,3},{\[Alpha],1,3}];
H\[Lambda]=Sum[\[Lambda] KroneckerProduct@@Join[{S[[4]]},Insert[sn,Max@bulkBonds[[\[Alpha],r]] ]@Insert[sn,Min@bulkBonds[[\[Alpha],r]]]@Table[s[[4]],N0-2]  ]  
,{r,1,Length@bulkBonds},{\[Alpha],1,3}];

HJ+H\[Lambda]]

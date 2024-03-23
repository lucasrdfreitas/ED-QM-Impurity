(* ::Package:: *)

(* ::Text:: *)
(*try to parallelize run multiples Mathematica files*)


(* ::Subsection:: *)
(*preamble*)


If[ \[Not]($FrontEnd===Null), SetDirectory[NotebookDirectory[]] ];
$FileName=If[$FrontEnd === Null, $InputFileName, NotebookFileName[] ];
Get[ FileNameJoin[{Directory[],"definitions.wl" }] ]


(* ::Subsection:: *)
(*Launching Kernel *)


If[ ($FrontEnd===Null),
	Print["Before Starting Kernels"];
	Needs["ClusterIntegration`"];
	(*kernels = LaunchKernels[SGE["micro4", 10]];*)
	Quiet[kernels = LaunchKernels[]];
	Print["Starting Kernels"];
];


Print[];
Print["Testing Mathematica Kernels "]
Print[];


Print["ProcessorCount=", $ProcessorCount ];

Print["$KernelID=",ParallelEvaluate[$KernelID] ];


(* ::Subsection:: *)
(*Running files*)


(* ::Text:: *)
(*a list with the name of the files to be executed *)


scripts= {
(*"adatom21.m" ,
"adatom21.m",
"substitution21.m",*)
"substitution21.m"
};


(* ::Text:: *)
(*loop executing the scripts in parallel*)


Print["        Executing the scripts in parallel"];Print[];
ParallelDo[

Print["----    Running script = ",scripts[[s]]," (",s,"/",Length@scripts,")    ----"];

Print["        Timing running ","(",s,"/",Length@scripts,") =",
				AbsoluteTiming[ Get[FileNameJoin[{Directory[],scripts[[s]] }]]; ]];

Print["----    Finishing script = ",scripts[[s]]," (",s,"/",Length@scripts,")    ----"];
Print[ ];Print[ ];
,{s,1,Length@scripts}]


(* ::Subsection:: *)
(*Closing Kernels*)


If[ ($FrontEnd===Null), 
	Print[];
	Print["Closing Kernels"];
	CloseKernels[];
];

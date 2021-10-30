
(* ::Package:: *)

(* ::Title:: *)
(*FTeV 0.13*)


(* ::Subtitle:: *)
(*Fast Tensors eValuation*)

(* ::Author:: *)
(*Davi C. Rodrigues*)
(*and Felipe Duarte dos Santos*)
(*Universiade Federal do Espirito Santo*)
(*This work was in part supported by  CNPq and FAPES.*)

(* ::Section:: *)
(*TO DO*)

(*
  - All FTeV functions should start with lower case letter
  - tensorPrint: incorporate option to set the indices positions. The default case should be all down. Exemplo: tensorPrint[tensor, "udu"]
  - tensorPrint: print indices positions (up or down) in the LHS
  - What could be the name of the FTeV function? tensorCompute or tC.
  - indices function should only require the target indices position (since the standard is all down)
  - T["ud", a,b] faz a conversão de "dd" para "ud". 
  - toTensor: should generate a function to yield the tensor components type T["dd",a,b] and reserve T["dd"] for the full tensor.
  - toTensor syntax:
    - toTensor[T[a], a] generates a List from the componets T[a]. By default, this tensor will be assumed to be covariant.
    - toTensor[T[a],a, "u"], the same of the previous case, but first it converts the current index status ("u") to covariant.
  - defTensor[List, T], from a given List, assumed to represent a covariant tensor, defines the symbol T as a tensor.
    This implies that T[] will be attributed to be the full covariant tensor.
    This defines T["udu"] to generate List with the new corresponding indices.
    This defines T["ddd", a,b] as the components of the fully covariant tensor
  - DCov syntax: 
    - DCov[T[], a] generates the component a of the covariant derivative of T (all indices down).
    - DCov[T[]] generates a List that fully represents the covariant derivative of T (all indices down).
  - Documentation
*)


(* ::Section:: *)
(*Begin package*)

BeginPackage["FTeV`"];

Print[Style["Starting FTeV 0.13 TEST", Bold]];
Print["Inside FTeV[] use: G for the Einstein Tensor; Riemann for Riemann tensor; Ricci for Ricci tensor; RicciS for Ricci Scalar; Chr for Christoffel Symbols; Weyl for Weyl tensor and kret for Kretschmann scalar."]
Print["The coordinates and the metric tensor must be named $Coordinates and $Metric, respectively."]


(* ::Subsection:: *)
(*Usage messages*)


FTeV::usage="FTeV[string] reads off the global $Coordinates and $Metric and takes in a string, which must be one of the following:
chr: For Christoffel symbols of the second kind (\!\(\*SubscriptBox[SuperscriptBox[\(gamma\), \(a\)], \(bc\)]\));
riemann: For Riemman tensor (\!\(\*SubscriptBox[SuperscriptBox[\(R\), \(d\)], \(abc\)]\));
ricci: For Covariant Ricci tensor (\!\(\*SubscriptBox[\(R\), \(ac\)]\));
riccis: For Ricci Scalar (R);
g: For Covariant Einstein Tensor (\!\(\*SubscriptBox[\(G\), \(ac\)]\));
weyl: For Weyl Tensor (\!\(\*SubscriptBox[\(C\), \(abcd\)]\));
kret: For Kretschmann scalar.";


indices::usage = "indices[\!\(\*
StyleBox[\"tensor\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"p_index\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"a_index\",\nFontSlant->\"Italic\"]\)] is used to raise or lower indices of tensors. Indices takes in a tensor as first argument, and two strings as the second and third arguments, where p_index indicates the indices prior to change and a_index indicates the desired change.";


TensorPrint::usage =
"print[\!\(\*
StyleBox[\"tensor\",\nFontSlant->\"Italic\"]\)] prints out the tensor in the form: [[i,j,..]]."


$Metric;
$Coordinates;
$PerturbationSymbol ;
$TensorPrintReplacements;
$TensorPrintReplacementsAlternative;

{g, gI, X, DCov, DDCov, spacetime, space, FirstOrder, FO, gamma, DimensionSpace, Contract, Metric, Coordinates, ToTensor, TensorBox, InsertTimeWeight};


(* ::Section:: *)
(*Begin Private*)


Begin["`Private`"];


(* ::Subsection:: *)
(*General purpose and "public" functions*)


g[a_,b_] := $Metric[[a+1,b+1]];
gInverse = Inverse@$Metric;
gI[a_,b_] := gInverse[[a+1,b+1]];
X[a_?NumberQ] := $Coordinates[[a + 1]];
spacetime := Sequence@@$Coordinates;
space := Sequence@@Drop[$Coordinates,1];
DimensionSpace := Length[$Coordinates]-1;


ClearAll[Contract];
SetAttributes[Contract,HoldAll];
SyntaxInformation[Contract] = {"ArgumentsPattern"->{_,__},"LocalVariables"->{"Integrate",{2,Infinity}}};
Contract[tensor_,indices__] := Module[
  {
    indList = {indices},
    indicesToSum
  },
  indicesToSum = Sequence @@ ({#,0,DimensionSpace}&) /@ indList;
  Sum[tensor,Evaluate[indicesToSum]]
];

ClearAll[ToTensor];
SetAttributes[ToTensor, HoldAll];
SyntaxInformation[ToTensor] = {"ArgumentsPattern"->{_,__},"LocalVariables"->{"Integrate",{2,Infinity}}};
ToTensor[tensorAbstractIndex_, indices__]:= Module[
  {
    indList={indices},
    indAux,
    indicesToRun
  },
  indicesToRun = Sequence @@ ({#,0,DimensionSpace}&)/@indList;
  Table[tensorAbstractIndex, Evaluate[indicesToRun]]
];


Clear[FirstOrder];
FO = FirstOrder;
FirstOrder := Series[#,{$PerturbationSymbol,0,1}] &;


InsertTimeWeight[expression_, weight_:Sqrt[$PerturbationSymbol]] := Module[
  {
    listOriginal,
    listFinal, 
    replace, 
    n,x,y,z,a,b,c,d,e,f 
  },
  listOriginal = Drop[FoldList[#1[#2] &, Derivative[n_,x_,y_,z_], {a__,b__,c__, d__, e__, f__}], 2];
  listFinal = weight^n Drop[FoldList[#1[#2]&, Derivative[n,x,y,z], {a,b,c, d,e,f}],2];
  AppendTo[listOriginal, Derivative[n_][x__][First@$Coordinates]];
  AppendTo[listFinal, weight^n Derivative[n][x][First@$Coordinates]];
  replace = Thread[listOriginal -> listFinal];
  expression /. replace
]; (*Derivatives on time receive a "weight", commonly useful for Newtonian or post-newtonian expansions. Works for time-only functions and 4D fields.*)


(* ::Subsection:: *)
(*FTeV function definition*)

Options[FTeV] = {Metric->$Metric, Coordinates->$Coordinates};
FTeV[string_?StringQ, OptionsPattern[]]:= Block[
  {
    coordvec = OptionValue[Coordinates],
    metric = OptionValue[Metric],
    d = Length[coordvec]-1,  (* number of space-time dimensions *) 
    Chr, 
    gI, 
    g,
    R,
    RS,
    G,
    x,
    i, j, k, l, m, n,
    Weyl,
    ruddd,
    ruuuu,
    rdddd,
    Kretsch
  },

  (* Defining coordinate componentes*)
  x[i_]:= coordvec[[i+1]];
  
  (* The metric and its inverse *)
  g[i_,j_] := metric[[i+1,j+1]];
  gI[i_,j_]:= Inverse[metric][[i+1,j+1]];

  (* Christoffel symbols *)
  Chr[k_,i_,j_]:= Chr[k,i,j] = 1/2 Sum[gI[k,l]  ( D[g[l,i],x[j]] + D[g[l,j],x[i]]- D[g[i,j],x[l]]),{l,0,d}];
  
  (* Riemman tensor *)
  R[l_,i_,k_,j_] := R[l,i,k,j] = D[Chr[l,i,j],x[k]] - D[Chr[l,i,k],x[j]] + Sum[Chr[l,m,k] Chr[m,i,j]- Chr[l,m,j]Chr[m,i, k], {m,0,d}];    
  
  (* Ricci tensor *)
  R[i_,j_] := R[i,j] = Sum[R[m,i,m,j], {m,0,d}];
  
  (* Ricci scalar*)
  RS := RS = Sum[R[i,m] gI[i, m], {m,0,d}, {i,0,d}];   
  
  (* Einstein tensor*)
  G[i_,j_] := G[i,j] = R[i,j] - 1/2  RS g[i,j];
    
  (* Weyl Tensor *)
  Weyl[l_,i_,k_,j_] := Weyl[l,i,k,j] = Sum[g[l,m]R[m,i,k,j],{m,0,d}]-1/(d-1) (R[l,k] g[i,j]+R[i,j] g[l,k]-R[l,j]g[i,k]-R[i,k]g[l,j])+1/(d (d-1)) RS(g[l,k] g[i,j] -g[i,k] g[l,j]);
  
  Kretsch:=(
    ruddd=Table[R[l,i,k,j],{l,0,d},{i,0,d},{j,0,d},{k,0,d}];
    ruuuu=indices[ruddd,"uddd","uuuu"];
    rdddd=indices[ruddd,"uddd","dddd"];
    Sum[ruuuu[[var1,var2,var3,var4]] rdddd[[var1,var2,var3,var4]] , {var1,d+1},{var2,d+1},{var3,d+1},{var4,d+1}]
  );

  Switch[ToLowerCase@string,
    "g", Table[G[i,j], {i,0,d},{j,0,d}],
    "chr",Table[Chr[k,i,j],{k,0,d},{i,0,d},{j,0,d}],
    "riemann",Table[R[l,i,k,j],{l,0,d},{i,0,d},{j,0,d},{k,0,d}],
    "ricci",Table[R[i,j],{i,0,d},{j,0,d}],
    "riccis", RS,
    "weyl", Table[Weyl[l,i,k,j],{l,0,d},{i,0,d},{j,0,d},{k,0,d}],
    "kret", Kretsch,
    _ , $Failed
  ]
];

(* ::Subsection:: *)
(*Covariant derivative (DCov)  & The Box (D'Alembertian) of a tensor (TensorBox)*)
(* DCov  is defined for each tensor rank independently. *)
(* The tensor to be derived is assumed to be covariant.*)


DCov::wrongRank ="The tensor rank should be one less than the number of indices used in dcov.";

(*For a scalar*)
DCov[scalar_][A_?NumberQ]:= D[scalar, X[A]] ;

(*For a vector*)
DCov[Vector_?ListQ][A_?NumberQ,B_?NumberQ] /; 
  If[Length@Dimensions@Vector == 1,True, Message[DCov::wrongRank]; False] := Module[
  {
    gammatensor = FTeV["Chr"],
    gamma,
    vec
  },
  vec[a_] := Vector[[a+1]];
  gamma[a_,b_,c_] := gammatensor[[a+1,b+1,c+1]];
  D[vec[B], X[A]] - Contract[gamma[a,A,B] vec[a], a]
];

(*For a rank 2 tensor*)
DCov[Tensor_?ListQ][A_?NumberQ,B_?NumberQ,C_?NumberQ] /;
  If[Length@Dimensions@Tensor == 2,True, Message[DCov::wrongRank]; False] := Module[
  {
    gammatensor = FTeV["Chr"],
    gamma,
    tensor
  },
  tensor[a_,b_]:= Tensor[[a+1,b+1]];
  gamma[a_,b_,c_]:= gammatensor[[a+1,b+1,c+1]];
  D[tensor[B,C],x[A]]- Contract[gamma[a,A,B]tensor[a, C],a] - Contract[ gamma[a,A,C]tensor[B,a], a]
];

(*For a rank 3 tensor*)
DCov[Tensor_?ListQ][A_?NumberQ,B_?NumberQ,C_?NumberQ, DD_?NumberQ] /; 
  If[Length@Dimensions@Tensor == 3,True, Message[DCov::wrongRank]; False] := Module[
  {
    gammatensor = FTeV["Chr"],
    gamma,
    tensor
  },
  tensor[a_,b_,c_]:= Tensor[[a+1,b+1, c+1]];
  gamma[a_,b_,c_]:= gammatensor[[a+1,b+1,c+1]];
  D[tensor[A, C, DD], x[A]]- Contract[gamma[a,A,B]tensor[a, C, DD],a] - Contract[gamma[a,A,C]tensor[B,a, DD],a]- Contract[gamma[a,A,DD]tensor[B,C, a] ,a]
];

(*This definition generates a tensor automatically, no indices are specified *)
DCov[tensor_][] := Module[
  {
    rank,
    GeneratedIndices,
    ind
  },
  If[ListQ@tensor,
    rank=TensorRank[tensor],
    (*else*)
    rank=0
  ];
  GeneratedIndices = Sequence@@(ind[#] & /@ Range[rank+1]);
  ToTensor[DCov[tensor][Evaluate[GeneratedIndices]], Evaluate[GeneratedIndices]]
];

DDCov[tensor_][] := DCov[DCov[tensor][]][]; (*Second covariant derivative*)

(*D'Alembertian for any tensor. The output is a tensor of the same rank of the input tensor.*)
TensorBox[tensor_]:= Block[
  {ddcov},
  ddcov[] = DCov[DCov[tensor][]][];
  ddcov[a_,b_]:=ddcov[][[a+1,b+1]];
  Contract[gI[a,b]ddcov[a,b],a,b]
];


(* ::Subsection:: *)
(*indices function: raises and lower indices.*)

Options[indices]={Metric->$Metric};
indices[tensor_?ListQ,string1_?StringQ,string2_?StringQ,OptionsPattern[]] := Block[
  {
    metric = OptionValue[Metric];,
    x,
    var,
    check, (* check defines which indices will be raised or lowered *)
    g,
    gI,
    Imetric = Inverse[metric],
    prodsum,
    varsum,
    tabsum
  },
  
  g[i_,j_]:= Indexed[metric,{i,j}];
  gI[i_,j_]:= Indexed[Imetric,{i,j}];

  Which[
    (*Option1*)
    Length[Characters[string1]]=!=Length[Characters[string2]],
    Return[$Failed],
    (*Option2*)
    Characters[string1]===Characters[string2],
    Return[tensor],
    (*Option3*)
    Characters[string1]=!=Characters[string2],
    check = If[
      Characters[string1][[#]] === Characters[string2][[#]], 
      Null, 
      (*else*) 
      If[Characters[string1][[#]] === "u" && Characters[string2][[#]] === "d", 
        g[Subscript[x, #], var[#]], 
        gI[Subscript[x, #], var[#]]
      ]
    ] & /@  Table[k, {k, Length @ Characters[string1]}];

    prodsum = Times @@ DeleteCases[check,Null] Indexed[tensor,
      {
        Apply[
          Sequence,
          If[check[[#]]=!= Null, var[#],Subscript[x, #]] &/@ Table[i,{i,1,Length @ Characters[string1]}]
        ]
      }
    ];
    
    varsum=Apply[Sequence,If[check[[#]]=!= Null, {var[#],1,Dimensions[tensor][[#]]},Nothing]&/@Table[i,{i,1,Length@Characters[string1]}]];

    tabsum=Apply[Sequence,{Subscript[x, #],1,Dimensions[tensor][[#]]}&/@Table[k,{k,1,Length@Characters[string1]}]];

    Table[Sum[prodsum,Evaluate@varsum],Evaluate@tabsum]
    (* evaluates the new tensor *)
  ]
];


(* ::Subsection:: *)
(*tensor print: a convenient display*)

$TensorPrintReplacements := {xx_[spacetime] ->  xx, xx_[space] -> xx };

$TensorPrintReplacementsAlternative := {xx_[spacetime] ->  xx[Style["\[Tau],x", Italic]], xx_[space] -> xx[Style["x",Italic]], xx_[numbers_?NumberQ] -> Subscript[xx, numbers] };

TensorPrint[tensor_] := Block[
  {
    variables,
    svariables,
    var,
    SimplifiedComponent
  },
  If[ListQ@tensor,
    Do[
      variables = Sequence@@Table[var[j], {j,1,TensorRank[tensor]}];
      svariables = StringDrop[StringDrop[ToString[{variables}-1],1],-1];
      SimplifiedComponent =Simplify@tensor[[variables]];
      If[SimplifiedComponent=!=0, (*Only prints the non-null components.*)
        Echo[ReplaceRepeated[SimplifiedComponent,  $TensorPrintReplacements], "\[LeftDoubleBracket]"<>svariables<>"\[RightDoubleBracket] : "] 
      ]; , 
      Evaluate[Sequence@@Table[{var[i],1,Dimensions[tensor][[i]]}, {i,1, TensorRank[tensor]}]]
    ],
    (*else*)
    Print[ReplaceRepeated[Simplify@tensor,  $TensorPrintReplacements]]
  ]
];


(* ::Section:: *)
(*End*)


End[];

EndPackage[];


(*Seems useful the idea below, but it is necessary to implement a definition of a function with arbitrary arguemnts...*)

(*It works, but only one time. needs revision.*)
(*DefineTensor[tensorName_,expressionWithAbstractIndices_, indices__]:=
Module[{rank=Length[{indices}], list, index},
ClearAll[Evaluate[tensorName]];
Evaluate[tensorName[]] = ToTensor[expressionWithAbstractIndices, indices];
Evaluate[kkk] = tensorName[][[Sequence@@(List[##] +1)]]&;
]*)


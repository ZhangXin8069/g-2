(* ::Package:: *)

(* Definetheparameters*)
X1=0;
X2=2;
Y1= X2 + 0.01;
Y2 = 30 + 0.01;
n = (X2 - X1)*1000 + 1;
m = (Y2 - Y1)*100 + 1;
hx = (X2 - X1)/(n - 1);
hy = (Y2 - Y1)/(m - 1);
x =Table[X1 + i*hx, {i, 0, n - 1}];
y =Table[Y1 + i*hy, {i, 0, m - 1}];


(* Define test functions *)
fV1 = 0.22; m1 = 0.78; gamma1 = 0.15;
fV2 = 0.19; m2 = 1.46; gamma2 = 0.4;
fV3 = 0.14; m3 = 1.72; gamma3 = 0.25;
fV4 = 0.14; m4 = 1.90; gamma4 = 0.1;


(* Define the F function *)
f[x_, M_, G_, F_] := F^2*((M^4 + M^2*G^2)/(24*Pi*((x - M^2)^2 + M^2*G^2)))
f1 = f[x, m1, gamma1, 1];
f2 = f[x, m2, gamma2, 1];
f3 = f[x, m3, gamma3, 1];
f4 = f[x, m4, gamma4, 1];
F = f1 + f2 + f3 + f4;


(* Define the trapz function *)
trapz[y, x] := Module[{},
  ExternalEvaluate[
    {"Python", "Executable" -> "/home/zhangxin/anaconda3/bin/python"},
    StringTemplate["
import numpy as np;
print(np.trapz(`MMAy`, `MMAx`));
"][
    <|
    "MMAy" -> ExportString[y, "PythonExpression"],
    "MMAx" -> ExportString[x, "PythonExpression"]|>
    ]
  ]
]
trapzG[y, x, F, m] := Module[{},
  ExternalEvaluate[
    {"Python", "Executable" -> "/home/zhangxin/anaconda3/bin/python"},
    StringTemplate["
import numpy as np;
y=np.array(`MMAy`);
x=np.array(`MMAx`);
F=np.array(`MMAF`);
m=int(`MMAm`);
G=[];
for i in range(m):
    G.append(np.trapz(1 / (y[i] - x) * F, x));
print(G);
"][
    <|
    "MMAy" -> ExportString[y, "PythonExpression"],
    "MMAx" -> ExportString[x, "PythonExpression"],
    "MMAF" -> ExportString[F, "PythonExpression"],
    "MMAm" -> ExportString[m, "PythonExpression"]|>
    ]
  ]
]
trapzA[Aphi, y, n] := Module[{},
  ExternalEvaluate[
    {"Python", "Executable" -> "/home/zhangxin/anaconda3/bin/python"},
    StringTemplate["
import numpy as np
Aphi=np.array(`MMAAphi`);
y=np.array(`MMAy`);
n=int(`MMAn`);
A = np.zeros((n, n))
for i in range(n):
    for j in range(n):
        A[i, j] = trapz(Aphi[:, i] * Aphi[:, j], y)
A
"][
    <|
    "MMAAphi" -> ExportString[Aphi "PythonExpression"],
    "MMAy" -> ExportString[y, "PythonExpression"],
    "MMAn" -> ExportString[n, "PythonExpression"]|>
    ]
  ]
]


(* Calculate G1 and G2 *)
F1 = f1;
G1 = trapzG[y, x, F1, m];
F2 = x - x;
G2 = trapzG[y, x, F2, m];
B1 = 1;
B2 = 2;
XD = 1.01;
b1 = XD*B1;
b2 = XD*B2;
fExact = B1*F1 + B2*F2;
fDelta = b1*F1 + b2*F2;


(* Calculate G3 *)
MM = fDelta[[1]];
NN = fDelta[[-1]];
VExact = (MM*(x - X2))/(X1 - X2) + (NN*(x - X1))/(X2 - X1);
G3 =Table[trapz[1/(y[[i]] - x)*VExact, x], {i, 1, m}];
UExact = fExact - VExact;


(* Calculate A and ba *)
A =Table[0, {n - 2}, {n - 2}];
Aphi =Table[0, {m}, {n - 2}];
ba =Table[0, {n - 2}];
For[i = 2, i <= n - 1, i++,
  Aphi[[All, i - 1]] = 
    1/hx*((y - x[[i - 1]])*Log[y - x[[i - 1]]] + (x[[i - 1]] - x[[i]]) + (y - x[[i + 1]])*Log[y - x[[i + 1]]] - 
       2*(y - x[[i]])*Log[y - x[[i]]] - (x[[i]] - x[[i + 1]]));
]
Aphi[[All, 1]] = Log[y - x[[1]]] + 
   1/hx*((y - x[[2]])*Log[y - x[[2]]] - (y - x[[1]])*Log[y - x[[1]]] - (x[[1]] - x[[2]]));
Aphi[[All, n - 2]] = -Log[y - x[[-1]]] + 
   1/hx*(-(y - x[[-1]])*Log[y - x[[-1]]] + (y - x[[-2]])*Log[y - x[[-2]]] + (x[[-2]] - x[[-1]]));
For[i = 1, i <= n - 2, i++,
  For[j = 1, j <= n - 2, j++,
    A[[i, j]] = trapz[Aphi[[All, i]]*Aphi[[All, j]], {y, Y1, Y2}];
  ];
  ba[[i]] = trapz[Aphi[[All, i]]*GDelta, {y, Y1, Y2}];
];


(* Create M matrix *)
M =Table[0, {n}, {n}];
M[[1, 1]] = 4;
M[[1, 2]] = 2;
For[i = 2, i <= n - 1, i++,
  M[[i, i]] = 8;
  M[[i, i - 1]] = 2;
  M[[i, i + 1]] = 2;
];
M[[-1, -1]] = 4;
M[[-1, -2]] = 2;
M = hx/12*M;



(* Create M1 matrix *)
M1 =Table[0, {n}, {n}];
M1[[1, 1]] = 1;
M1[[1, 2]] = -1;
For[i = 2, i <= n - 1, i++,
  M1[[i, i]] = 2;
  M1[[i, i - 1]] = -1;
  M1[[i, i + 1]] = -1;
];
M1[[-1, -1]] = 1;
M1[[-1, -2]] = -1;
M1 = (1/hx)*M1;


(* Create H1 matrix *)
H1 = M + M1;


(* Remove first and last rows and columns from A, ba, and H1 *)
A = A[[2 ;; -2, 2 ;; -2]];
ba = ba[[2 ;; -2]];
H1 = H1[[2 ;; -2, 2 ;; -2]];


(* Set the number of iterations *)
ite = 20;


(* Initialize arrays for results *)
alpha =Table[0, {ite}];
ua =Table[0, {n - 2}, {ite}];
Ua =Table[0, {n}, {ite}];
fa =Table[0, {n}, {ite}];
res =Table[0, {ite}];
err =Table[0, {ite}];
L =Table[0, {ite}];
LH =Table[0, {ite}];


(* Perform iterations *)
For[k = 1, k <= ite, k++,
  alpha[[k]] = 10^(-k);
  ua[[All, k]] = LinearSolve[A + alpha[[k]]*H1, ba];
  Ua[[All, k]] = Join[{0}, ua[[All, k]], {0}];
  fa[[All, k]] = Ua[[All, k]] + VExact;
  res[[k]] = Sqrt[trapz[(Aphi[[All, All, k]].fa[[All, k]] - gDelta)^2, {y, Y1, Y2}]];
  err[[k]] = Sqrt[trapz[(fa[[All, k]] - fExact)^2, x]];
  L[[k]] = Sqrt[trapz[fa[[All, k]]^2, x]]*res[[k]];
  LH[[k]] = Sqrt[trapz[fa[[All, k]]^2 + Grad[fa[[All, k]], {x}][[1]]^2, x]]*res[[k]];
  (* Plot the results *)
  Print[ListLinePlot[{x, fExact, fa[[All, k]]}, PlotStyle -> {Automatic, Red, Blue}, 
     PlotLegends -> {"Exact", "Approximated"}, 
     PlotLabel -> "alpha = " <> ToString[alpha[[k]]]]];
]



(* Find the indices of minimum L and LH *)
kl = FirstPosition[L, Min[L]][[1]];
klh = FirstPosition[LH, Min[LH]][[1]];


(* Calculate V and Rho *)
V =Table[0, {ite}];
For[i = 1, i <= ite, i++,
  VV = IdentityMatrix[m] - Aphi.Transpose[Aphi].Inverse[Transpose[Aphi].Aphi + alpha[[i]] IdentityMatrix[n - 2]].Transpose[Aphi];
  V1 = Norm[VV.gDelta]^2;
  V2 = Tr[VV]^2;
  V[[i]] = V1/V2;
]

Rho =Table[0, {ite}];
For[i = 1, i <= ite, i++,
  Rho[[i]] = Norm[alpha[[i]] Inverse[Transpose[Aphi].Aphi + alpha[[i]] IdentityMatrix[n - 2]].fa[[All, i]]];
]



(* Print the results *)
{{"L-curve+L^2", "L-curve+H^1", "GCV", "拟最优"}, {kl, klh, FirstPosition[V, Min[V]][[1]], FirstPosition[Rho, Min[Rho]][[1]]}}


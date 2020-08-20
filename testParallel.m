(* Wolfram Language package *)
1 + 1
ParallelMap[# + 1 &, Range[10]]
testfun[x_] := Module[{y} ,
  y = x^2;
  y + 1]
ParallelEvaluate[testfun[10]]
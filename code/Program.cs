using MathNet.Numerics.LinearAlgebra;

#region Inputs
var input = File.ReadAllLines(@"G:\My Drive\WN23\561 Core Des\HW\9\input.txt");

var h = double.Parse(input[0]);

var split = input[1].Split(' ', StringSplitOptions.RemoveEmptyEntries);
var alpha = split.Select(double.Parse).ToArray();

var k = int.Parse(input[2]);
var m = int.Parse(input[3]);
var n = k * m;
var dh = h / n;

var ncomp = int.Parse(input[4]);
var comp = new Composition[ncomp];
for (var i = 0; i < ncomp; ++i) {
  split = input[5 + i].Split(' ', StringSplitOptions.RemoveEmptyEntries);
  comp[i] = new Composition {
    d = double.Parse(split[0]),
    sigmaA = double.Parse(split[1]),
    nuSigmaF = double.Parse(split[2])
  };
}

var assignments = new int[k];
split = input[5 + ncomp].Split(' ', StringSplitOptions.RemoveEmptyEntries);
for (var i = 0; i < k; ++i) {
  assignments[i] = int.Parse(split[i]) - 1;
}

var normal = double.Parse(input[6 + ncomp]);

split = input[7 + ncomp].Split(' ', StringSplitOptions.RemoveEmptyEntries);
var kConverge = double.Parse(split[0]);
var fConverge = double.Parse(split[1]);
#endregion

// *** HELPER FUNCTIONS //
var getParam = (int i) => comp[assignments[i / m]];
var dTilde = (int i) => 1 / (dh / (2 * getParam(i).d) + dh / (2 * getParam(i + 1).d));
var dTildeBounds = (int i, Boundary a) => 1 / (dh / (2 * getParam(i).d) + 1 / alpha[(int)a]);
// *** //

var phi = Vector<double>.Build.DenseOfEnumerable(Enumerable.Repeat(1d, n));
var lambda = 1d;

#region Mesh Balance Equations
var a = new double[n];
var b = new double[n];
var c = new double[n];
var d = new double[n];

for (var i = 0; i < n; ++i) {
  d[i] = getParam(i).nuSigmaF * dh;

  if (i == 0) {
    // left boundary

    var c1 = -dTilde(i);
    b[i] = getParam(i).sigmaA * dh + dTildeBounds(i, Boundary.Left) - c1;
    c[i] = c1;

  } else if (i == n - 1) {
    // right boundary

    var an = -dTilde(i - 1);
    b[i] = getParam(i).sigmaA * dh + dTildeBounds(i, Boundary.Right) - an;
    a[i] = an;

  } else {
    // everything else

    var ai = -dTilde(i - 1);
    var ci = -dTilde(i);

    b[i] = getParam(i).sigmaA * dh - ai - ci;
    a[i] = ai;
    c[i] = ci;
  }
}
#endregion

var bigF = Matrix<double>.Build.SparseOfDiagonalArray(d);
var error = 1d;

#region Outer Iteration
for (var l = 0; ; ++l) {

  #region LU Factorization
  var aTilde = new double[n];
  var bTilde = new double[n];

  bTilde[0] = b[0];

  for (var i = 1; i < n; ++i) {
    aTilde[i] = a[i] / bTilde[i - 1];
    bTilde[i] = b[i] - aTilde[i] * c[i - 1];
  }
  
  var psi = (lambda * (bigF * phi));
  #endregion

  #region Forward Elimination
  var y = new double[n];

  for (var i = 0; i < n; ++i) {
    y[i] = psi[i];

    if (i != 0) {
      y[i] -= aTilde[i] * y[i - 1];
    }
  }
  #endregion

  #region Backward Substitution
  var newPhi = Vector<double>.Build.Dense(n);

  for (var i = n - 1; i >= 0; --i) {
    newPhi[i] = y[i] / bTilde[i];

    if (i != n - 1) {
      newPhi[i] -= (c[i] * newPhi[i + 1]) / bTilde[i];
    }
  }
  #endregion

  var newLambda = lambda * newPhi.DotProduct(phi) / newPhi.DotProduct(newPhi);
  var oldK = 1 / lambda;
  var newK = 1 / newLambda;

  var kDiff = Math.Abs(newK - oldK);
  var infNorm = newPhi.Zip(phi, (a, b) => Math.Abs(a - b)).Max();

  var converged = kDiff < kConverge && infNorm < fConverge;
  var ratio = infNorm / error;

  Console.WriteLine($"{kDiff,10:F8}\t{infNorm,10:F8}\t{ratio,10:F8}");
  
  if (converged) {
    #region Output
    Console.WriteLine($"\nConverged after {l} iterations.");
    Console.WriteLine($"k = {newK:f5}, dominance = {ratio:f8}");

    var i = 0;
    var phiNorm = newPhi.Sum();
    var phiString = newPhi.Aggregate("", (a, b) => a + $"{dh * i++}\t{b / phiNorm:F6}\n");
    phiString += $"{dh * i}\t0\n";
    File.WriteAllText(@"G:\My Drive\WN23\561 Core Des\HW\9\phi.txt", phiString);

    i = 0;
    var sourceNorm = d.Sum();
    var sourceString = d.Aggregate("", (a, b) => a + $"{dh * i++}\t{b / sourceNorm:F6}\n");
    sourceString += $"{dh * i}\t0\n";
    File.WriteAllText(@"G:\My Drive\WN23\561 Core Des\HW\9\source.txt", sourceString);
    #endregion

    break;
  }

  lambda = newLambda;
  phi = newPhi;
  error = infNorm;
}
#endregion

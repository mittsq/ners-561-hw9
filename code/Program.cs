using MathNet.Numerics.LinearAlgebra;

// Wielandt Shift Skip Value
//    set to `int.MaxValue` to disable
const int SKIP = int.MaxValue;

#region Inputs
var input = File.ReadAllLines(@"G:\My Drive\WN23\561 Core Des\HW\9\input-360.txt");

var h = double.Parse(input[0]);

var split = input[1].Split(' ', StringSplitOptions.RemoveEmptyEntries);
var alpha = split.Select(double.Parse).ToArray();

var r = int.Parse(input[2]);
var m = int.Parse(input[3]);
var n = r * m;
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

var assignments = new int[r];
split = input[5 + ncomp].Split(' ', StringSplitOptions.RemoveEmptyEntries);
for (var i = 0; i < r; ++i) {
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
var shift = 0d;
var k = 1d;

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
  // perform wielandt shift
  var kTilde = 1 / (1 / k - shift);
  
  var bShift = (int i) => b[i] - shift * d[i];

  #region LU Factorization
  var aTilde = new double[n];
  var bTilde = new double[n];

  bTilde[0] = bShift(0);

  for (var i = 1; i < n; ++i) {
    aTilde[i] = a[i] / bTilde[i - 1];
    bTilde[i] = bShift(i) - aTilde[i] * c[i - 1];
  }
  
  var psi = (1 / kTilde) * (bigF * phi);
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

  var source = bigF * newPhi;
  var gamma = newPhi.DotProduct(phi) / newPhi.DotProduct(newPhi);
  var newK = 1 / (gamma / k + shift * (1 - gamma));
  var newKTilde = kTilde / gamma;

  var newShift = 0d;
  if (l >= SKIP) {
    newShift = shift + 1 / newKTilde - 1e-5;
  }

  var kDiff = Math.Abs(newK - k);
  var infNorm = newPhi.Zip(phi, (a, b) => Math.Abs(a - b)).Max();

  var converged = kDiff < kConverge && infNorm < fConverge;
  var ratio = infNorm / error;

  Console.WriteLine($"delta k   = {kDiff,10:F6}");
  Console.WriteLine($"infNorm   = {infNorm,10:F6}");
  Console.WriteLine($"ratio     = {ratio,10:F6}");
  Console.WriteLine($"k         = {newK,10:F6}");
  Console.WriteLine($"shift     = {newShift,10:F6}");
  Console.WriteLine();

  if (double.IsNaN(kDiff)) {
    break;
  }
  
  if (converged) {
    #region Output
    Console.WriteLine($"\nConverged after {l} iterations.");
    Console.WriteLine($"k = {newK:f5}, dominance = {ratio:f8}");

    var i = 0;
    var phiNorm = newPhi.Sum();
    var phiString = newPhi.Aggregate("", (a, b) => a + $"{dh * i++}\t{b / phiNorm:F6}\n");
    phiString += $"{dh * i}\t{0:F6}\n";
    File.WriteAllText(@$"G:\My Drive\WN23\561 Core Des\HW\9\phi-{h}-{SKIP}.txt", phiString);

    i = 0;
    var sourceNorm = d.Sum();
    var sourceString = d.Aggregate("", (a, b) => a + $"{dh * i++}\t{b / sourceNorm:F6}\n");
    sourceString += $"{dh * i}\t{0:F6}\n";
    File.WriteAllText(@$"G:\My Drive\WN23\561 Core Des\HW\9\source-{h}-{SKIP}.txt", sourceString);
    #endregion

    break;
  }

  k = newK;
  phi = newPhi;
  error = infNorm;
  shift = newShift;
}
#endregion

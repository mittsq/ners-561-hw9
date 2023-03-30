#region Constants
using MathNet.Numerics.LinearAlgebra;
#endregion

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

#region Getters
var getParam = (int i) => comp[assignments[i / m]];
var dTilde = (int i) => 1 / (dh / (2 * getParam(i).d) + dh / (2 * getParam(i + 1).d));
var dTildeBounds = (int i, Boundary a) => 1 / (dh / (2 * getParam(i).d) + 1 / alpha[(int)a]);
#endregion

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
var bigM = Matrix<double>.Build.Sparse(n, n);

for (var i = 0; i < n; ++i) {
  bigM[i, i] = b[i];

  if (i != 0 && i != n - 1) {
    bigM[i, i - 1] = a[i];
    bigM[i, i + 1] = c[i];
  }
}

for (var l = 0; ; ++l) {
  var newPhi = lambda * (bigM.Inverse() * bigF) * phi;
  var newLambda = lambda * newPhi.DotProduct(phi) / newPhi.DotProduct(newPhi);

  var kDiff = Math.Abs(newLambda.Flip() - lambda.Flip());
  var infNorm = newPhi.Zip(phi, (a, b) => Math.Abs(a - b)).Max();

  Console.WriteLine($"kDiff = {kDiff,10:F7}, infNorm = {infNorm,10:F6}");

  if (kDiff < kConverge && infNorm < fConverge) {
    Console.WriteLine($"\nConverged after {l} iterations.");
    Console.WriteLine($"k = {newLambda.Flip():F5}");

    var i = 0;
    var phiString = phi.Aggregate(string.Empty, (a, b) => a + $"{i++}\t{b:F6}\n");
    File.WriteAllText(@"G:\My Drive\WN23\561 Core Des\HW\9\phi.txt", phiString);
    break;
  }

  phi = newPhi;
  lambda = newLambda;
}

#region Output
#endregion

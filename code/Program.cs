#region Constants
using MathNet.Numerics.LinearAlgebra;

const int LEFT = 0;
const int RIGHT = 1;
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
var dTildeBounds = (int i, int a) => 1 / (dh / (2 * getParam(i).d) + 1 / alpha[a]);
#endregion

var phi = Enumerable.Repeat(1d, n).ToArray();
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
    b[i] = getParam(i).sigmaA * dh + dTildeBounds(i, LEFT) - c1;
    c[i + 1] = c1;

  } else if (i == n - 1) {
    // right boundary

    var an = -dTilde(i - 1);
    b[i] = getParam(i).sigmaA * dh + dTildeBounds(i, RIGHT) - an;
    a[i - 1] = an;

  } else {
    // everything else

    var ai = -dTilde(i - 1);
    var ci = -dTilde(i);

    b[i] = getParam(i).sigmaA * dh - ai - ci;
    a[i - 1] = ai;
    c[i + 1] = ci;
  }
}
#endregion

var bigF = Matrix<double>.Build.SparseOfDiagonalArray(d);

while (true) {
  #region Inner Iteration
  var phiV = Vector<double>.Build.DenseOfArray(phi);

  #region LU Factorization
  var aTilde = new double[n];
  var bTilde = new double[n];
  var cTilde = new double[n];

  bTilde[0] = b[0];

  for (var i = 1; i < n; ++i) {
    aTilde[i - 1] = a[i - 1] / bTilde[i - 1];
    bTilde[i] = b[i] - aTilde[i] * c[i];
  }

  var psi = (lambda * (bigF * phiV)).ToArray();
  #endregion

  #region Forward Elimination
  var y = new double[n];

  for (var i = 0; i < n; ++i) {
    y[i] = psi[i];

    if (i != 0) {
      y[i] += aTilde[i - 1] * y[i - 1];
    }
  }
  #endregion

  #region Backward Substitution
  var newPhi = new double[n];

  for (var i = n - 1; i >= 0; --i) {
    var invB = 1 / bTilde[i];
    newPhi[i] = y[i] * invB;

    if (i != n - 1) {
      newPhi[i] -= cTilde[i + 1] * newPhi[i + 1] * invB;
    }
  }
  #endregion

  var newPhiV = Vector<double>.Build.DenseOfArray(newPhi);
  var newLambda = lambda * newPhiV.DotProduct(phiV) / newPhiV.DotProduct(newPhiV);
  var newK = 1 / newLambda;

  // loop, terminate if converged
  var converged = Math.Abs(newK - 1 / lambda) < kConverge;
  var infNorm = newPhi.Zip(phi, (a, b) => Math.Abs(a - b)).Max();
  converged &= infNorm < fConverge;

  Console.WriteLine($"{newK:f6}\t{infNorm:f6}");
  if (converged) {
    Console.WriteLine();
    
    for (var i = 0; i < n; ++i) {
      Console.WriteLine($"{newPhi[i]:f6}");
    }

    break;
  }

  lambda = newLambda;
  phi = newPhi;

  #endregion
}

#region Output
#endregion

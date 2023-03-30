struct Composition {
  public double d;
  public double sigmaA;
  public double nuSigmaF;
}

enum Boundary {
  Left, Right
}

static class Extensions {
  public static double Flip(this double d) => 1 / d;
}

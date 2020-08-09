#ifndef STIFF_DECSOL_HPP
#define STIFF_DECSOL_HPP

#include "stiff/vector_matrix.hpp"

namespace stiff {

template <typename T, int SizeN = Eigen::Dynamic, int SizeM = Eigen::Dynamic>
int dec(Matrix<T, SizeN, SizeM> &a, Vector<int, SizeN> &ip) {
  /* System generated locals */
  int a_dim1, a_offset, i1, i2, i3;
  T d1, d2;

  /* Local variables */
  int i, j, k, m;
  T t;
  int nm1, kp1, ier;

  /* VERSION REAL DOUBLE PRECISION */
  /* ----------------------------------------------------------------------- */
  /*  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION. */
  /*  INPUT.. */
  /*     N = ORDER OF MATRIX. */
  /*     NDIM = DECLARED DIMENSION OF ARRAY  A . */
  /*     A = MATRIX TO BE TRIANGULARIZED. */
  /*  OUTPUT.. */
  /*     A(I,J), I.LE.J = UPPER TRIANGULAR FACTOR, U . */
  /*     A(I,J), I.GT.J = MULTIPLIERS = LOWER TRIANGULAR FACTOR, I - L. */
  /*     IP(K), K.LT.N = INDEX OF K-TH PIVOT ROW. */
  /*     IP(N) = (-1)**(NUMBER OF INTERCHANGES) OR O . */
  /*     IER = 0 IF MATRIX A IS NONSINGULAR, OR K IF FOUND TO BE */
  /*           SINGULAR AT STAGE K. */
  /*  USE  SOL  TO OBTAIN SOLUTION OF LINEAR SYSTEM. */
  /*  DETERM(A) = IP(N)*A(1,1)*A(2,2)*...*A(N,N). */
  /*  IF IP(N)=O, A IS SINGULAR, SOL WILL DIVIDE BY ZERO. */

  /*  REFERENCE.. */
  /*     C. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER, */
  /*     C.A.C.M. 15 (1972), P. 274. */
  /* ----------------------------------------------------------------------- */

  const int n = a.rows();
  const int ndim = a.cols();
  /* Function Body */
  ier = 0;
  ip(n - 1) = 1;
  if (n == 1) {
    goto L70;
  }
  nm1 = n - 1;
  i1 = nm1;
  for (k = 1; k <= i1; ++k) {
    kp1 = k + 1;
    m = k;
    i2 = n;
    for (i = kp1; i <= i2; ++i) {
      if ((d1 = a(i - 1, k - 1), abs(d1)) > (d2 = a(m - 1, k - 1), abs(d2))) {
        m = i;
      }
      /* L10: */
    }
    ip(k - 1) = m;
    t = a(m - 1, k - 1);
    if (m == k) {
      goto L20;
    }
    ip(n - 1) = -ip(n - 1);
    a(m - 1, k - 1) = a(k - 1, k - 1);
    a(k - 1, k - 1) = t;
  L20:
    if (t == 0.) {
      goto L80;
    }
    t = 1. / t;
    i2 = n;
    for (i = kp1; i <= i2; ++i) {
      /* L30: */
      a(i - 1, k - 1) = -a(i - 1, k - 1) * t;
    }
    i2 = n;
    for (j = kp1; j <= i2; ++j) {
      t = a(m - 1, j - 1);
      a(m - 1, j - 1) = a(k - 1, j - 1);
      a(k - 1, j - 1) = t;
      if (t == 0.) {
        goto L45;
      }
      i3 = n;
      for (i = kp1; i <= i3; ++i) {
        /* L40: */
        a(i - 1, j - 1) += a(i - 1, k - 1) * t;
      }
    L45:
        /* L50: */
        ;
    }
    /* L60: */
  }
L70:
  k = n;
  if (a(n - 1, n - 1) == 0.) {
    goto L80;
  }
  return ier;
L80:
  ier = k;
  ip(n - 1) = 0;
  return ier;
  /* ----------------------- END OF SUBROUTINE DEC ------------------------- */
}

template <typename T, int SizeN = Eigen::Dynamic, int SizeM = Eigen::Dynamic>
void sol(const Matrix<T, SizeN, SizeM> &a, Vector<T, SizeN> &b,
         const Vector<int, SizeN> &ip) {
  /* System generated locals */
  int a_dim1, a_offset, i1, i2;

  /* Local variables */
  int i, k, m;
  T t;
  int kb, km1, nm1, kp1;

  /* VERSION REAL DOUBLE PRECISION */
  /* ----------------------------------------------------------------------- */
  /*  SOLUTION OF LINEAR SYSTEM, A*X = B . */
  /*  INPUT.. */
  /*    N = ORDER OF MATRIX. */
  /*    NDIM = DECLARED DIMENSION OF ARRAY  A . */
  /*    A = TRIANGULARIZED MATRIX OBTAINED FROM DEC. */
  /*    B = RIGHT HAND SIDE VECTOR. */
  /*    IP = PIVOT VECTOR OBTAINED FROM DEC. */
  /*  DO NOT USE IF DEC HAS SET IER .NE. 0. */
  /*  OUTPUT.. */
  /*    B = SOLUTION VECTOR, X . */
  /* ----------------------------------------------------------------------- */
  /* Parameter adjustments */
  const int n = a.rows();
  const int ndim = a.cols();

  /* Function Body */
  if (n == 1) {
    goto L50;
  }
  nm1 = n - 1;
  i1 = nm1;
  for (k = 1; k <= i1; ++k) {
    kp1 = k + 1;
    m = ip(k - 1);
    t = b(m - 1);
    b(m - 1) = b(k - 1);
    b(k - 1) = t;
    i2 = n;
    for (i = kp1; i <= i2; ++i) {
      /* L10: */
      b(i - 1) += a(i - 1, k - 1) * t;
    }
    /* L20: */
  }
  i1 = nm1;
  for (kb = 1; kb <= i1; ++kb) {
    km1 = n - kb;
    k = km1 + 1;
    b(k - 1) /= a(k - 1, k - 1);
    t = -b(k - 1);
    i2 = km1;
    for (i = 1; i <= i2; ++i) {
      /* L30: */
      b(i - 1) += a(i - 1, k - 1) * t;
    }
    /* L40: */
  }
L50:
  b(1 - 1) /= a(1 - 1, 1 - 1);
}

template <typename T, int SizeN = Eigen::Dynamic, int SizeM = Eigen::Dynamic>
int dech(Matrix<T, SizeN, SizeM> &a, const int lb, Vector<int, SizeN> &ip) {
  /* System generated locals */
  int a_dim1, a_offset, i1, i2, i3;
  T d1, d2;

  /* Local variables */
  int i, j, k, m;
  T t;
  int na, nm1, kp1, ier;

  /* VERSION REAL DOUBLE PRECISION */
  /* ----------------------------------------------------------------------- */
  /*  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION OF A HESSENBERG */
  /*  MATRIX WITH LOWER BANDWIDTH LB */
  /*  INPUT.. */
  /*     N = ORDER OF MATRIX A. */
  /*     NDIM = DECLARED DIMENSION OF ARRAY  A . */
  /*     A = MATRIX TO BE TRIANGULARIZED. */
  /*     LB = LOWER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED, LB.GE.1). */
  /*  OUTPUT.. */
  /*     A(I,J), I.LE.J = UPPER TRIANGULAR FACTOR, U . */
  /*     A(I,J), I.GT.J = MULTIPLIERS = LOWER TRIANGULAR FACTOR, I - L. */
  /*     IP(K), K.LT.N = INDEX OF K-TH PIVOT ROW. */
  /*     IP(N) = (-1)**(NUMBER OF INTERCHANGES) OR O . */
  /*     IER = 0 IF MATRIX A IS NONSINGULAR, OR K IF FOUND TO BE */
  /*           SINGULAR AT STAGE K. */
  /*  USE  SOLH  TO OBTAIN SOLUTION OF LINEAR SYSTEM. */
  /*  DETERM(A) = IP(N)*A(1,1)*A(2,2)*...*A(N,N). */
  /*  IF IP(N)=O, A IS SINGULAR, SOL WILL DIVIDE BY ZERO. */

  /*  REFERENCE.. */
  /*     THIS IS A SLIGHT MODIFICATION OF */
  /*     C. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER, */
  /*     C.A.C.M. 15 (1972), P. 274. */
  /* ----------------------------------------------------------------------- */
  /* Parameter adjustments */
  const int n = a.rows();
  const int ndim = a.cols();

  /* Function Body */
  ier = 0;
  ip(n - 1) = 1;
  if (n == 1) {
    goto L70;
  }
  nm1 = n - 1;
  i1 = nm1;
  for (k = 1; k <= i1; ++k) {
    kp1 = k + 1;
    m = k;
    /* Computing MIN */
    i2 = n, i3 = lb + k;
    na = std::min(i2, i3);
    i2 = na;
    for (i = kp1; i <= i2; ++i) {
      if ((d1 = a(i - 1, k - 1), abs(d1)) > (d2 = a(m - 1, k - 1), abs(d2))) {
        m = i;
      }
      /* L10: */
    }
    ip(k - 1) = m;
    t = a(m - 1, k - 1);
    if (m == k) {
      goto L20;
    }
    ip(n - 1) = -ip(n - 1);
    a(m - 1, k - 1) = a(k - 1, k - 1);
    a(k - 1, k - 1) = t;
  L20:
    if (t == 0.) {
      goto L80;
    }
    t = 1. / t;
    i2 = na;
    for (i = kp1; i <= i2; ++i) {
      /* L30: */
      a(i - 1, k - 1) = -a(i - 1, k - 1) * t;
    }
    i2 = n;
    for (j = kp1; j <= i2; ++j) {
      t = a(m - 1, j - 1);
      a(m - 1, j - 1) = a(k - 1, j - 1);
      a(k - 1, j - 1) = t;
      if (t == 0.) {
        goto L45;
      }
      i3 = na;
      for (i = kp1; i <= i3; ++i) {
        /* L40: */
        a(i - 1, j - 1) += a(i - 1, k - 1) * t;
      }
    L45:
        /* L50: */
        ;
    }
    /* L60: */
  }
L70:
  k = n;
  if (a(n - 1, n - 1) == 0.) {
    goto L80;
  }
  return ier;
L80:
  ier = k;
  ip(n - 1) = 0;
  return ier;
}

template <typename T, int SizeN = Eigen::Dynamic, int SizeM = Eigen::Dynamic>
void solh(const Matrix<T, SizeN, SizeM> &a, const int lb, Vector<T, SizeN> &b,
          const Vector<int, SizeN> &ip) {
  /* System generated locals */
  int a_dim1, a_offset, i1, i2, i3;

  /* Local variables */
  int i, k, m;
  T t;
  int kb, na, km1, nm1, kp1;

  /* VERSION REAL DOUBLE PRECISION */
  /* ----------------------------------------------------------------------- */
  /*  SOLUTION OF LINEAR SYSTEM, A*X = B . */
  /*  INPUT.. */
  /*    N = ORDER OF MATRIX A. */
  /*    NDIM = DECLARED DIMENSION OF ARRAY  A . */
  /*    A = TRIANGULARIZED MATRIX OBTAINED FROM DECH. */
  /*    LB = LOWER BANDWIDTH OF A. */
  /*    B = RIGHT HAND SIDE VECTOR. */
  /*    IP = PIVOT VECTOR OBTAINED FROM DEC. */
  /*  DO NOT USE IF DECH HAS SET IER .NE. 0. */
  /*  OUTPUT.. */
  /*    B = SOLUTION VECTOR, X . */
  /* ----------------------------------------------------------------------- */
  /* Parameter adjustments */
  const int n = a.rows();
  const int ndim = a.cols();

  /* Function Body */
  if (n == 1) {
    goto L50;
  }
  nm1 = n - 1;
  i1 = nm1;
  for (k = 1; k <= i1; ++k) {
    kp1 = k + 1;
    m = ip(k - 1);
    t = b(m - 1);
    b(m - 1) = b(k - 1);
    b(k - 1) = t;
    /* Computing MIN */
    i2 = n, i3 = lb + k;
    na = std::min(i2, i3);
    i2 = na;
    for (i = kp1; i <= i2; ++i) {
      /* L10: */
      b(i - 1) += a(i - 1, k - 1) * t;
    }
    /* L20: */
  }
  i1 = nm1;
  for (kb = 1; kb <= i1; ++kb) {
    km1 = n - kb;
    k = km1 + 1;
    b(k - 1) /= a(k - 1, k - 1);
    t = -b(k - 1);
    i2 = km1;
    for (i = 1; i <= i2; ++i) {
      /* L30: */
      b(i - 1) += a(i - 1, k - 1) * t;
    }
    /* L40: */
  }
L50:
  b(1 - 1) /= a[a_dim1 + 1];
}

template <typename T, int SizeN = Eigen::Dynamic, int SizeM = Eigen::Dynamic>
int decc(Matrix<T, SizeN, SizeM> &ar, Matrix<T, SizeN, SizeM> &ai,
         Vector<int, SizeN> &ip) {
  /* System generated locals */
  int ar_dim1, ar_offset, ai_dim1, ai_offset, i1, i2, i3;
  T d1, d2, d__3, d__4;

  /* Local variables */
  int i, j, k, m;
  T ti, tr;
  int nm1, kp1, ier;
  T den, prodi, prodr;

  /* VERSION COMPLEX DOUBLE PRECISION */
  /* ----------------------------------------------------------------------- */
  /*  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION */
  /*  ------ MODIFICATION FOR COMPLEX MATRICES -------- */
  /*  INPUT.. */
  /*     N = ORDER OF MATRIX. */
  /*     NDIM = DECLARED DIMENSION OF ARRAYS  AR AND AI . */
  /*     (AR, AI) = MATRIX TO BE TRIANGULARIZED. */
  /*  OUTPUT.. */
  /*     AR(I,J), I.LE.J = UPPER TRIANGULAR FACTOR, U ; REAL PART. */
  /*     AI(I,J), I.LE.J = UPPER TRIANGULAR FACTOR, U ; IMAGINARY PART. */
  /*     AR(I,J), I.GT.J = MULTIPLIERS = LOWER TRIANGULAR FACTOR, I - L. */
  /*                                                    REAL PART. */
  /*     AI(I,J), I.GT.J = MULTIPLIERS = LOWER TRIANGULAR FACTOR, I - L. */
  /*                                                    IMAGINARY PART. */
  /*     IP(K), K.LT.N = INDEX OF K-TH PIVOT ROW. */
  /*     IP(N) = (-1)**(NUMBER OF INTERCHANGES) OR O . */
  /*     IER = 0 IF MATRIX A IS NONSINGULAR, OR K IF FOUND TO BE */
  /*           SINGULAR AT STAGE K. */
  /*  USE  SOL  TO OBTAIN SOLUTION OF LINEAR SYSTEM. */
  /*  IF IP(N)=O, A IS SINGULAR, SOL WILL DIVIDE BY ZERO. */

  /*  REFERENCE.. */
  /*     C. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER, */
  /*     C.A.C.M. 15 (1972), P. 274. */
  /* ----------------------------------------------------------------------- */
  /* Parameter adjustments */
  const int n = ar.rows();
  const int ndim = ar.cols();

  /* Function Body */
  ier = 0;
  ip(n - 1) = 1;
  if (n == 1) {
    goto L70;
  }
  nm1 = n - 1;
  i1 = nm1;
  for (k = 1; k <= i1; ++k) {
    kp1 = k + 1;
    m = k;
    i2 = n;
    for (i = kp1; i <= i2; ++i) {
      if ((d1 = ar(i - 1, k - 1), abs(d1)) + (d2 = ai(i - 1, k - 1), abs(d2)) >
          (d__3 = ar(m - 1, k - 1), abs(d__3)) +
              (d__4 = ai(m - 1, k - 1), abs(d__4))) {
        m = i;
      }
      /* L10: */
    }
    ip(k - 1) = m;
    tr = ar(m - 1, k - 1);
    ti = ai(m - 1, k - 1);
    if (m == k) {
      goto L20;
    }
    ip(n - 1) = -ip(n - 1);
    ar(m - 1, k - 1) = ar(k - 1, k - 1);
    ai(m - 1, k - 1) = ai(k - 1, k - 1);
    ar(k - 1, k - 1) = tr;
    ai(k - 1, k - 1) = ti;
  L20:
    if (abs(tr) + abs(ti) == 0.) {
      goto L80;
    }
    den = tr * tr + ti * ti;
    tr /= den;
    ti = -ti / den;
    i2 = n;
    for (i = kp1; i <= i2; ++i) {
      prodr = ar(i - 1, k - 1) * tr - ai(i - 1, k - 1) * ti;
      prodi = ai(i - 1, k - 1) * tr + ar(i - 1, k - 1) * ti;
      ar(i - 1, k - 1) = -prodr;
      ai(i - 1, k - 1) = -prodi;
      /* L30: */
    }
    i2 = n;
    for (j = kp1; j <= i2; ++j) {
      tr = ar(m - 1, j - 1);
      ti = ai(m - 1, j - 1);
      ar(m - 1, j - 1) = ar(k - 1, j - 1);
      ai(m - 1, j - 1) = ai(k - 1, j - 1);
      ar(k - 1, j - 1) = tr;
      ai(k - 1, j - 1) = ti;
      if (abs(tr) + abs(ti) == 0.) {
        goto L48;
      }
      if (ti == 0.) {
        i3 = n;
        for (i = kp1; i <= i3; ++i) {
          prodr = ar(i - 1, k - 1) * tr;
          prodi = ai(i - 1, k - 1) * tr;
          ar(i - 1, j - 1) += prodr;
          ai(i - 1, j - 1) += prodi;
          /* L40: */
        }
        goto L48;
      }
      if (tr == 0.) {
        i3 = n;
        for (i = kp1; i <= i3; ++i) {
          prodr = -ai(i - 1, k - 1) * ti;
          prodi = ar(i - 1, k - 1) * ti;
          ar(i - 1, j - 1) += prodr;
          ai(i - 1, j - 1) += prodi;
          /* L45: */
        }
        goto L48;
      }
      i3 = n;
      for (i = kp1; i <= i3; ++i) {
        prodr = ar(i - 1, k - 1) * tr - ai(i - 1, k - 1) * ti;
        prodi = ai(i - 1, k - 1) * tr + ar(i - 1, k - 1) * ti;
        ar(i - 1, j - 1) += prodr;
        ai(i - 1, j - 1) += prodi;
        /* L47: */
      }
    L48:
        /* L50: */
        ;
    }
    /* L60: */
  }
L70:
  k = n;
  if ((d1 = ar(n - 1, n - 1), abs(d1)) + (d2 = ai(n - 1, n - 1), abs(d2)) ==
      0.) {
    goto L80;
  }
  return ier;
L80:
  ier = k;
  ip(n - 1) = 0;
  return ier;
}

template <typename T, int SizeN = Eigen::Dynamic, int SizeM = Eigen::Dynamic>
void solc(const Matrix<T, SizeN, SizeM> &ar, const Matrix<T, SizeN, SizeM> &ai,
          Vector<T, SizeN> &br, Vector<T, SizeN> &bi,
          const Vector<int, SizeN> &ip) {
  /* System generated locals */
  int ar_dim1, ar_offset, ai_dim1, ai_offset, i1, i2;

  /* Local variables */
  int i, k, m, kb;
  T ti, tr;
  int km1, nm1, kp1;
  T den, prodi, prodr;

  /* VERSION COMPLEX DOUBLE PRECISION */
  /* ----------------------------------------------------------------------- */
  /*  SOLUTION OF LINEAR SYSTEM, A*X = B . */
  /*  INPUT.. */
  /*    N = ORDER OF MATRIX. */
  /*    NDIM = DECLARED DIMENSION OF ARRAYS  AR AND AI. */
  /*    (AR,AI) = TRIANGULARIZED MATRIX OBTAINED FROM DEC. */
  /*    (BR,BI) = RIGHT HAND SIDE VECTOR. */
  /*    IP = PIVOT VECTOR OBTAINED FROM DEC. */
  /*  DO NOT USE IF DEC HAS SET IER .NE. 0. */
  /*  OUTPUT.. */
  /*    (BR,BI) = SOLUTION VECTOR, X . */
  /* ----------------------------------------------------------------------- */
  /* Parameter adjustments */
  const int n = ar.rows();

  /* Function Body */
  if (n == 1) {
    goto L50;
  }
  nm1 = n - 1;
  i1 = nm1;
  for (k = 1; k <= i1; ++k) {
    kp1 = k + 1;
    m = ip(k - 1);
    tr = br(m - 1);
    ti = bi(m - 1);
    br(m - 1) = br(k - 1);
    bi(m - 1) = bi(k - 1);
    br(k - 1) = tr;
    bi(k - 1) = ti;
    i2 = n;
    for (i = kp1; i <= i2; ++i) {
      prodr = ar(i - 1, k - 1) * tr - ai(i - 1, k - 1) * ti;
      prodi = ai(i - 1, k - 1) * tr + ar(i - 1, k - 1) * ti;
      br(i - 1) += prodr;
      bi(i - 1) += prodi;
      /* L10: */
    }
    /* L20: */
  }
  i1 = nm1;
  for (kb = 1; kb <= i1; ++kb) {
    km1 = n - kb;
    k = km1 + 1;
    den = ar(k - 1, k - 1) * ar(k - 1, k - 1) +
          ai(k - 1, k - 1) * ai(k - 1, k - 1);
    prodr = br(k - 1) * ar(k - 1, k - 1) + bi(k - 1) * ai(k - 1, k - 1);
    prodi = bi(k - 1) * ar(k - 1, k - 1) - br(k - 1) * ai(k - 1, k - 1);
    br(k - 1) = prodr / den;
    bi(k - 1) = prodi / den;
    tr = -br(k - 1);
    ti = -bi(k - 1);
    i2 = km1;
    for (i = 1; i <= i2; ++i) {
      prodr = ar(i - 1, k - 1) * tr - ai(i - 1, k - 1) * ti;
      prodi = ai(i - 1, k - 1) * tr + ar(i - 1, k - 1) * ti;
      br(i - 1) += prodr;
      bi(i - 1) += prodi;
      /* L30: */
    }
    /* L40: */
  }
L50:
  den = ar(0, 0) * ar(0, 0) + ai(0, 0) * ai(0, 0);
  prodr = br(0) * ar(0, 0) + bi(0) * ai(0, 0);
  prodi = bi(0) * ar(0, 0) - br(0) * ai(0, 0);
  br(0) = prodr / den;
  bi(0) = prodi / den;
}

template <typename T, int SizeN = Eigen::Dynamic, int SizeM = Eigen::Dynamic>
int dechc(Matrix<T, SizeN, SizeM> &ar, Matrix<T, SizeN, SizeM> &ai,
          const int lb, Vector<int, SizeN> &ip) {
  /* System generated locals */
  int ar_dim1, ar_offset, ai_dim1, ai_offset, i1, i2, i3;
  T d1, d2, d__3, d__4;

  /* Local variables */
  int i, j, k, m, na;
  T ti, tr;
  int nm1, kp1, ier;
  T den, prodi, prodr;

  /* VERSION COMPLEX DOUBLE PRECISION */
  /* ----------------------------------------------------------------------- */
  /*  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION */
  /*  ------ MODIFICATION FOR COMPLEX MATRICES -------- */
  /*  INPUT.. */
  /*     N = ORDER OF MATRIX. */
  /*     NDIM = DECLARED DIMENSION OF ARRAYS  AR AND AI . */
  /*     (AR, AI) = MATRIX TO BE TRIANGULARIZED. */
  /*  OUTPUT.. */
  /*     AR(I,J), I.LE.J = UPPER TRIANGULAR FACTOR, U ; REAL PART. */
  /*     AI(I,J), I.LE.J = UPPER TRIANGULAR FACTOR, U ; IMAGINARY PART. */
  /*     AR(I,J), I.GT.J = MULTIPLIERS = LOWER TRIANGULAR FACTOR, I - L. */
  /*                                                    REAL PART. */
  /*     AI(I,J), I.GT.J = MULTIPLIERS = LOWER TRIANGULAR FACTOR, I - L. */
  /*                                                    IMAGINARY PART. */
  /*     LB = LOWER BANDWIDTH OF A (DIAGONAL NOT COUNTED), LB.GE.1. */
  /*     IP(K), K.LT.N = INDEX OF K-TH PIVOT ROW. */
  /*     IP(N) = (-1)**(NUMBER OF INTERCHANGES) OR O . */
  /*     IER = 0 IF MATRIX A IS NONSINGULAR, OR K IF FOUND TO BE */
  /*           SINGULAR AT STAGE K. */
  /*  USE  SOL  TO OBTAIN SOLUTION OF LINEAR SYSTEM. */
  /*  IF IP(N)=O, A IS SINGULAR, SOL WILL DIVIDE BY ZERO. */

  /*  REFERENCE.. */
  /*     C. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER, */
  /*     C.A.C.M. 15 (1972), P. 274. */
  /* ----------------------------------------------------------------------- */
  const int n = ar.row();

  /* Function Body */
  ier = 0;
  ip(n - 1) = 1;
  if (lb == 0) {
    goto L70;
  }
  if (n == 1) {
    goto L70;
  }
  nm1 = n - 1;
  i1 = nm1;
  for (k = 1; k <= i1; ++k) {
    kp1 = k + 1;
    m = k;
    /* Computing MIN */
    i2 = n, i3 = lb + k;
    na = std::min(i2, i3);
    i2 = na;
    for (i = kp1; i <= i2; ++i) {
      if ((d1 = ar(i - 1, k - 1), abs(d1)) + (d2 = ai(i - 1, k - 1), abs(d2)) >
          (d__3 = ar(m - 1, k - 1), abs(d__3)) +
              (d__4 = ai(m - 1, k - 1), abs(d__4))) {
        m = i;
      }
      /* L10: */
    }
    ip(k - 1) = m;
    tr = ar(m - 1, k - 1);
    ti = ai(m - 1, k - 1);
    if (m == k) {
      goto L20;
    }
    ip(n - 1) = -ip(n - 1);
    ar(m - 1, k - 1) = ar(k - 1, k - 1);
    ai(m - 1, k - 1) = ai(k - 1, k - 1);
    ar(k - 1, k - 1) = tr;
    ai(k - 1, k - 1) = ti;
  L20:
    if (abs(tr) + abs(ti) == 0.) {
      goto L80;
    }
    den = tr * tr + ti * ti;
    tr /= den;
    ti = -ti / den;
    i2 = na;
    for (i = kp1; i <= i2; ++i) {
      prodr = ar(i - 1, k - 1) * tr - ai(i - 1, k - 1) * ti;
      prodi = ai(i - 1, k - 1) * tr + ar(i - 1, k - 1) * ti;
      ar(i - 1, k - 1) = -prodr;
      ai(i - 1, k - 1) = -prodi;
      /* L30: */
    }
    i2 = n;
    for (j = kp1; j <= i2; ++j) {
      tr = ar(m - 1, j - 1);
      ti = ai(m - 1, j - 1);
      ar(m - 1, j - 1) = ar(k - 1, j - 1);
      ai(m - 1, j - 1) = ai(k - 1, j - 1);
      ar(k - 1, j - 1) = tr;
      ai(k - 1, j - 1) = ti;
      if (abs(tr) + abs(ti) == 0.) {
        goto L48;
      }
      if (ti == 0.) {
        i3 = na;
        for (i = kp1; i <= i3; ++i) {
          prodr = ar(i - 1, k - 1) * tr;
          prodi = ai(i - 1, k - 1) * tr;
          ar(i - 1, j - 1) += prodr;
          ai(i - 1, j - 1) += prodi;
          /* L40: */
        }
        goto L48;
      }
      if (tr == 0.) {
        i3 = na;
        for (i = kp1; i <= i3; ++i) {
          prodr = -ai(i - 1, k - 1) * ti;
          prodi = ar(i - 1, k - 1) * ti;
          ar(i - 1, j - 1) += prodr;
          ai(i - 1, j - 1) += prodi;
          /* L45: */
        }
        goto L48;
      }
      i3 = na;
      for (i = kp1; i <= i3; ++i) {
        prodr = ar(i - 1, k - 1) * tr - ai(i - 1, k - 1) * ti;
        prodi = ai(i - 1, k - 1) * tr + ar(i - 1, k - 1) * ti;
        ar(i - 1, j - 1) += prodr;
        ai(i - 1, j - 1) += prodi;
        /* L47: */
      }
    L48:;
    }
  }
L70:
  k = n;
  if ((d1 = ar(n - 1, n - 1), abs(d1)) + (d2 = ai(n - 1, n - 1), abs(d2)) ==
      0.0) {
    goto L80;
  }
  return ier;
L80:
  ier = k;
  ip(n - 1) = 0;
  return ier;
}

template <typename T, int SizeN = Eigen::Dynamic, int SizeM = Eigen::Dynamic>
void solhc(const Matrix<T, SizeN, SizeM> &ar, const Matrix<T, SizeN, SizeM> &ai,
           const int lb, Vector<T, SizeN> &br, Vector<T, SizeN> &bi,
           const Vector<int, SizeN> &ip) {
  /* System generated locals */
  int ar_dim1, ar_offset, ai_dim1, ai_offset, i1, i2, i3, i4;

  /* Local variables */
  int i, k, m, kb;
  T ti, tr;
  int km1, nm1, kp1;
  T den, prodi, prodr;

  /* VERSION COMPLEX DOUBLE PRECISION */
  /* ----------------------------------------------------------------------- */
  /*  SOLUTION OF LINEAR SYSTEM, A*X = B . */
  /*  INPUT.. */
  /*    N = ORDER OF MATRIX. */
  /*    NDIM = DECLARED DIMENSION OF ARRAYS  AR AND AI. */
  /*    (AR,AI) = TRIANGULARIZED MATRIX OBTAINED FROM DEC. */
  /*    (BR,BI) = RIGHT HAND SIDE VECTOR. */
  /*    LB = LOWER BANDWIDTH OF A. */
  /*    IP = PIVOT VECTOR OBTAINED FROM DEC. */
  /*  DO NOT USE IF DEC HAS SET IER .NE. 0. */
  /*  OUTPUT.. */
  /*    (BR,BI) = SOLUTION VECTOR, X . */
  /* ----------------------------------------------------------------------- */

  const int n = ar.row();

  /* Function Body */
  if (n == 1) {
    goto L50;
  }
  nm1 = n - 1;
  if (lb == 0) {
    goto L25;
  }
  i1 = nm1;
  for (k = 1; k <= i1; ++k) {
    kp1 = k + 1;
    m = ip(k - 1);
    tr = br(m - 1);
    ti = bi(m - 1);
    br(m - 1) = br(k - 1);
    bi(m - 1) = bi(k - 1);
    br(k - 1) = tr;
    bi(k - 1) = ti;
    /* Computing MIN */
    i3 = n, i4 = lb + k;
    i2 = std::min(i3, i4);
    for (i = kp1; i <= i2; ++i) {
      prodr = ar(i - 1, k - 1) * tr - ai(i - 1, k - 1) * ti;
      prodi = ai(i - 1, k - 1) * tr + ar(i - 1, k - 1) * ti;
      br(i - 1) += prodr;
      bi(i - 1) += prodi;
      /* L10: */
    }
    /* L20: */
  }
L25:
  i1 = nm1;
  for (kb = 1; kb <= i1; ++kb) {
    km1 = n - kb;
    k = km1 + 1;
    den = ar(k - 1, k - 1) * ar(k - 1, k - 1) +
          ai(k - 1, k - 1) * ai(k - 1, k - 1);
    prodr = br(k - 1) * ar(k - 1, k - 1) + bi(k - 1) * ai(k - 1, k - 1);
    prodi = bi(k - 1) * ar(k - 1, k - 1) - br(k - 1) * ai(k - 1, k - 1);
    br(k - 1) = prodr / den;
    bi(k - 1) = prodi / den;
    tr = -br(k - 1);
    ti = -bi(k - 1);
    i2 = km1;
    for (i = 1; i <= i2; ++i) {
      prodr = ar(i - 1, k - 1) * tr - ai(i - 1, k - 1) * ti;
      prodi = ai(i - 1, k - 1) * tr + ar(i - 1, k - 1) * ti;
      br(i - 1) += prodr;
      bi(i - 1) += prodi;
      /* L30: */
    }
    /* L40: */
  }
L50:
  den = ar(0, 0) * ar(0, 0) + ai(0, 0) * ai(0, 0);
  prodr = br(0) * ar(0, 0) + bi(0) * ai(0, 0);
  prodi = bi(0) * ar(0, 0) - br(0) * ai(0, 0);
  br(0) = prodr / den;
  bi(0) = prodi / den;
}

template <typename T, int SizeN = Eigen::Dynamic, int SizeM = Eigen::Dynamic>
int decb(Matrix<T, SizeN, SizeM> &a, const int ml, const int mu,
         Vector<int, SizeN> &ip) {
  /* System generated locals */
  int a_dim1, a_offset, i1, i2, i3, i4;
  T d1, d2;

  /* Local variables */
  int i, j, k, m;
  T t;
  int md, jk, mm, ju, md1, nm1, kp1, mdl, ijk, ier;

  /* ----------------------------------------------------------------------- */
  /*  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION OF A BANDED */
  /*  MATRIX WITH LOWER BANDWIDTH ML AND UPPER BANDWIDTH MU */
  /*  INPUT.. */
  /*     N       ORDER OF THE ORIGINAL MATRIX A. */
  /*     NDIM    DECLARED DIMENSION OF ARRAY  A. */
  /*     A       CONTAINS THE MATRIX IN BAND STORAGE.   THE COLUMNS */
  /*                OF THE MATRIX ARE STORED IN THE COLUMNS OF  A  AND */
  /*                THE DIAGONALS OF THE MATRIX ARE STORED IN ROWS */
  /*                ML+1 THROUGH 2*ML+MU+1 OF  A. */
  /*     ML      LOWER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED). */
  /*     MU      UPPER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED). */
  /*  OUTPUT.. */
  /*     A       AN UPPER TRIANGULAR MATRIX IN BAND STORAGE AND */
  /*                THE MULTIPLIERS WHICH WERE USED TO OBTAIN IT. */
  /*     IP      INDEX VECTOR OF PIVOT INDICES. */
  /*     IP(N)   (-1)**(NUMBER OF INTERCHANGES) OR O . */
  /*     IER     = 0 IF MATRIX A IS NONSINGULAR, OR  = K IF FOUND TO BE */
  /*                SINGULAR AT STAGE K. */
  /*  USE  SOLB  TO OBTAIN SOLUTION OF LINEAR SYSTEM. */
  /*  DETERM(A) = IP(N)*A(MD,1)*A(MD,2)*...*A(MD,N)  WITH MD=ML+MU+1. */
  /*  IF IP(N)=O, A IS SINGULAR, SOLB WILL DIVIDE BY ZERO. */

  /*  REFERENCE.. */
  /*     THIS IS A MODIFICATION OF */
  /*     C. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER, */
  /*     C.A.C.M. 15 (1972), P. 274. */
  /* ----------------------------------------------------------------------- */
  const int n = a.rows();

  /* Function Body */
  ier = 0;
  ip(n - 1) = 1;
  md = ml + mu + 1;
  md1 = md + 1;
  ju = 0;
  if (ml == 0) {
    goto L70;
  }
  if (n == 1) {
    goto L70;
  }
  if (n < mu + 2) {
    goto L7;
  }
  i1 = n;
  for (j = mu + 2; j <= i1; ++j) {
    i2 = ml;
    for (i = 1; i <= i2; ++i) {
      a(i - 1, j - 1) = 0.;
    }
  }
L7:
  nm1 = n - 1;
  i2 = nm1;
  for (k = 1; k <= i2; ++k) {
    kp1 = k + 1;
    m = md;
    /* Computing MIN */
    i1 = ml, i3 = n - k;
    mdl = std::min(i1, i3) + md;
    i1 = mdl;
    for (i = md1; i <= i1; ++i) {
      if ((d1 = a(i - 1, k - 1), abs(d1)) > (d2 = a(m - 1, k - 1), abs(d2))) {
        m = i;
      }
      /* L10: */
    }
    ip(k - 1) = m + k - md;
    t = a(m - 1, k - 1);
    if (m == md) {
      goto L20;
    }
    ip(n - 1) = -ip(n - 1);
    a(m - 1, k - 1) = a(md - 1, k - 1);
    a(md - 1, k - 1) = t;
  L20:
    if (t == 0.) {
      goto L80;
    }
    t = 1. / t;
    i1 = mdl;
    for (i = md1; i <= i1; ++i) {
      /* L30: */
      a(i - 1, k - 1) = -a(i - 1, k - 1) * t;
    }
    /* Computing MIN */
    /* Computing MAX */
    i3 = ju, i4 = mu + ip(k - 1);
    i1 = std::max(i3, i4);
    ju = std::min(i1, n);
    mm = md;
    if (ju < kp1) {
      goto L55;
    }
    i1 = ju;
    for (j = kp1; j <= i1; ++j) {
      --m;
      --mm;
      t = a(m - 1, j - 1);
      if (m == mm) {
        goto L35;
      }
      a(m - 1, j - 1) = a(mm - 1, j - 1);
      a(mm - 1, j - 1) = t;
    L35:
      if (t == 0.) {
        goto L45;
      }
      jk = j - k;
      i3 = mdl;
      for (i = md1; i <= i3; ++i) {
        ijk = i - jk;
        /* L40: */
        a(ijk - 1, j - 1) += a(i - 1, k - 1) * t;
      }
    L45:;
    }
  L55:;
  }
L70:
  k = n;
  if (a(md - 1, n - 1) == 0.0) {
    goto L80;
  }
  return ier;
L80:
  ier = k;
  ip(n - 1) = 0;
  return ier;
}

template <typename T, int SizeN = Eigen::Dynamic, int SizeM = Eigen::Dynamic>
void solb(const Matrix<T, SizeN, SizeM> &a, const int ml, const int mu,
          Vector<T, SizeN> &b, const Vector<int, SizeN> &ip) {
  /* System generated locals */
  int a_dim1, a_offset, i1, i2, i3;

  /* Local variables */
  int i, k, m;
  T t;
  int kb, md, lm, md1, nm1, imd, kmd, mdl, mdm, ier;

  /* ----------------------------------------------------------------------- */
  /*  SOLUTION OF LINEAR SYSTEM, A*X = B . */
  /*  INPUT.. */
  /*    N      ORDER OF MATRIX A. */
  /*    NDIM   DECLARED DIMENSION OF ARRAY  A . */
  /*    A      TRIANGULARIZED MATRIX OBTAINED FROM DECB. */
  /*    ML     LOWER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED). */
  /*    MU     UPPER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED). */
  /*    B      RIGHT HAND SIDE VECTOR. */
  /*    IP     PIVOT VECTOR OBTAINED FROM DECB. */
  /*  DO NOT USE IF DECB HAS SET IER .NE. 0. */
  /*  OUTPUT.. */
  /*    B      SOLUTION VECTOR, X . */
  /* ----------------------------------------------------------------------- */
  const int n = a.rows();

  /* Function Body */
  md = ml + mu + 1;
  md1 = md + 1;
  mdm = md - 1;
  nm1 = n - 1;
  if (ml == 0) {
    goto L25;
  }
  if (n == 1) {
    goto L50;
  }
  i1 = nm1;
  for (k = 1; k <= i1; ++k) {
    m = ip(k - 1);
    t = b(m - 1);
    b(m - 1) = b(k - 1);
    b(k - 1) = t;
    /* Computing MIN */
    i2 = ml, i3 = n - k;
    mdl = std::min(i2, i3) + md;
    i2 = mdl;
    for (i = md1; i <= i2; ++i) {
      imd = i + k - md;
      b(imd - 1) += a(i - 1, k - 1) * t;
    }
  }
L25:
  i1 = nm1;
  for (kb = 1; kb <= i1; ++kb) {
    k = n + 1 - kb;
    b(k - 1) /= a(md - 1, k - 1);
    t = -b(k - 1);
    kmd = md - k;
    /* Computing MAX */
    i2 = 1, i3 = kmd + 1;
    lm = std::max(i2, i3);
    i2 = mdm;
    for (i = lm; i <= i2; ++i) {
      imd = i - kmd;
      /* L30: */
      b(imd - 1) += a(i - 1, k - 1) * t;
    }
    /* L40: */
  }
L50:
  b(0) /= a(md - 1, 0);
}

template <typename T, int SizeN = Eigen::Dynamic, int SizeM = Eigen::Dynamic>
int decbc(Matrix<T, SizeN, SizeM> &ar, Matrix<T, SizeN, SizeM> &ai,
          const int ml, const int mu, Vector<int, SizeN> &ip) {
  /* System generated locals */
  int ar_dim1, ar_offset, ai_dim1, ai_offset, i1, i2, i3, i4;
  T d1, d2, d__3, d__4;

  /* Local variables */
  int i, j, k, m, md, jk, mm;
  T ti;
  int ju;
  T tr;
  int md1, nm1, kp1;
  T den;
  int mdl, ijk, ier;
  T prodi, prodr;

  /* ----------------------------------------------------------------------- */
  /*  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION OF A BANDED COMPLEX */
  /*  MATRIX WITH LOWER BANDWIDTH ML AND UPPER BANDWIDTH MU */
  /*  INPUT.. */
  /*     N       ORDER OF THE ORIGINAL MATRIX A. */
  /*     NDIM    DECLARED DIMENSION OF ARRAY  A. */
  /*     AR, AI     CONTAINS THE MATRIX IN BAND STORAGE.   THE COLUMNS */
  /*                OF THE MATRIX ARE STORED IN THE COLUMNS OF  AR (REAL */
  /*                PART) AND AI (IMAGINARY PART)  AND */
  /*                THE DIAGONALS OF THE MATRIX ARE STORED IN ROWS */
  /*                ML+1 THROUGH 2*ML+MU+1 OF  AR AND AI. */
  /*     ML      LOWER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED). */
  /*     MU      UPPER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED). */
  /*  OUTPUT.. */
  /*     AR, AI  AN UPPER TRIANGULAR MATRIX IN BAND STORAGE AND */
  /*                THE MULTIPLIERS WHICH WERE USED TO OBTAIN IT. */
  /*     IP      INDEX VECTOR OF PIVOT INDICES. */
  /*     IP(N)   (-1)**(NUMBER OF INTERCHANGES) OR O . */
  /*     IER     = 0 IF MATRIX A IS NONSINGULAR, OR  = K IF FOUND TO BE */
  /*                SINGULAR AT STAGE K. */
  /*  USE  SOLBC  TO OBTAIN SOLUTION OF LINEAR SYSTEM. */
  /*  DETERM(A) = IP(N)*A(MD,1)*A(MD,2)*...*A(MD,N)  WITH MD=ML+MU+1. */
  /*  IF IP(N)=O, A IS SINGULAR, SOLBC WILL DIVIDE BY ZERO. */

  /*  REFERENCE.. */
  /*     THIS IS A MODIFICATION OF */
  /*     C. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER, */
  /*     C.A.C.M. 15 (1972), P. 274. */
  /* ----------------------------------------------------------------------- */
  const int n = ar.rows();

  /* Function Body */
  ier = 0;
  ip(n - 1) = 1;
  md = ml + mu + 1;
  md1 = md + 1;
  ju = 0;
  if (ml == 0) {
    goto L70;
  }
  if (n == 1) {
    goto L70;
  }
  if (n < mu + 2) {
    goto L7;
  }
  i1 = n;
  for (j = mu + 2; j <= i1; ++j) {
    i2 = ml;
    for (i = 1; i <= i2; ++i) {
      ar(i - 1, j - 1) = 0.;
      ai(i - 1, j - 1) = 0.;
      /* L5: */
    }
  }
L7:
  nm1 = n - 1;
  i2 = nm1;
  for (k = 1; k <= i2; ++k) {
    kp1 = k + 1;
    m = md;
    /* Computing MIN */
    i1 = ml, i3 = n - k;
    mdl = std::min(i1, i3) + md;
    i1 = mdl;
    for (i = md1; i <= i1; ++i) {
      if ((d1 = ar(i - 1, k - 1), abs(d1)) + (d2 = ai(i - 1, k - 1), abs(d2)) >
          (d__3 = ar(m - 1, k - 1), abs(d__3)) +
              (d__4 = ai(m - 1, k - 1), abs(d__4))) {
        m = i;
      }
      /* L10: */
    }
    ip(k - 1) = m + k - md;
    tr = ar(m - 1, k - 1);
    ti = ai(m - 1, k - 1);
    if (m == md) {
      goto L20;
    }
    ip(n - 1) = -ip(n - 1);
    ar(m - 1, k - 1) = ar(md - 1, k - 1);
    ai(m - 1, k - 1) = ai(md - 1, k - 1);
    ar(md - 1, k - 1) = tr;
    ai(md - 1, k - 1) = ti;
  L20:
    if (abs(tr) + abs(ti) == 0.) {
      goto L80;
    }
    den = tr * tr + ti * ti;
    tr /= den;
    ti = -ti / den;
    i1 = mdl;
    for (i = md1; i <= i1; ++i) {
      prodr = ar(i - 1, k - 1) * tr - ai(i - 1, k - 1) * ti;
      prodi = ai(i - 1, k - 1) * tr + ar(i - 1, k - 1) * ti;
      ar(i - 1, k - 1) = -prodr;
      ai(i - 1, k - 1) = -prodi;
      /* L30: */
    }
    /* Computing MIN */
    /* Computing MAX */
    i3 = ju, i4 = mu + ip(k - 1);
    i1 = std::max(i3, i4);
    ju = std::min(i1, n);
    mm = md;
    if (ju < kp1) {
      goto L55;
    }
    i1 = ju;
    for (j = kp1; j <= i1; ++j) {
      --m;
      --mm;
      tr = ar(m - 1, j - 1);
      ti = ai(m - 1, j - 1);
      if (m == mm) {
        goto L35;
      }
      ar(m - 1, j - 1) = ar(mm - 1, j - 1);
      ai(m - 1, j - 1) = ai(mm - 1, j - 1);
      ar(mm - 1, j - 1) = tr;
      ai(mm - 1, j - 1) = ti;
    L35:
      if (abs(tr) + abs(ti) == 0.) {
        goto L48;
      }
      jk = j - k;
      if (ti == 0.) {
        i3 = mdl;
        for (i = md1; i <= i3; ++i) {
          ijk = i - jk;
          prodr = ar(i - 1, k - 1) * tr;
          prodi = ai(i - 1, k - 1) * tr;
          ar(ijk - 1, j - 1) += prodr;
          ai(ijk - 1, j - 1) += prodi;
          /* L40: */
        }
        goto L48;
      }
      if (tr == 0.) {
        i3 = mdl;
        for (i = md1; i <= i3; ++i) {
          ijk = i - jk;
          prodr = -ai(i - 1, k - 1) * ti;
          prodi = ar(i - 1, k - 1) * ti;
          ar(ijk - 1, j - 1) += prodr;
          ai(ijk - 1, j - 1) += prodi;
          /* L45: */
        }
        goto L48;
      }
      i3 = mdl;
      for (i = md1; i <= i3; ++i) {
        ijk = i - jk;
        prodr = ar(i - 1, k - 1) * tr - ai(i - 1, k - 1) * ti;
        prodi = ai(i - 1, k - 1) * tr + ar(i - 1, k - 1) * ti;
        ar(ijk - 1, j - 1) += prodr;
        ai(ijk - 1, j - 1) += prodi;
      }
    L48:;
    }
  L55:;
  }
L70:
  k = n;
  if (abs(ar(md - 1, n - 1)) + abs(ai(md - 1, n - 1)) == 0.0) {
    goto L80;
  }
  return ier;
L80:
  ier = k;
  ip(n - 1) = 0;
  return ier;
}

template <typename T, int SizeN = Eigen::Dynamic, int SizeM = Eigen::Dynamic>
void solbc(const Matrix<T, SizeN, SizeM> &ar, const Matrix<T, SizeN, SizeM> &ai,
           const int ml, const int mu, Vector<T, SizeN> &br,
           Vector<T, SizeN> &bi, const Vector<int, SizeN> &ip) {
  /* System generated locals */
  int ar_dim1, ar_offset, ai_dim1, ai_offset, i1, i2, i3;

  /* Local variables */
  int i, k, m, kb, md, lm;
  T ti, tr;
  int md1, nm1;
  T den;
  int imd, kmd, mdl, mdm;
  T prodi, prodr;

  /* ----------------------------------------------------------------------- */
  /*  SOLUTION OF LINEAR SYSTEM, A*X = B , */
  /*                  VERSION BANDED AND COMPLEX-DOUBLE PRECISION. */
  /*  INPUT.. */
  /*    N      ORDER OF MATRIX A. */
  /*    NDIM   DECLARED DIMENSION OF ARRAY  A . */
  /*    AR, AI TRIANGULARIZED MATRIX OBTAINED FROM DECB (REAL AND IMAG. PART).
   */
  /*    ML     LOWER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED). */
  /*    MU     UPPER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED). */
  /*    BR, BI RIGHT HAND SIDE VECTOR (REAL AND IMAG. PART). */
  /*    IP     PIVOT VECTOR OBTAINED FROM DECBC. */
  /*  DO NOT USE IF DECB HAS SET IER .NE. 0. */
  /*  OUTPUT.. */
  /*    BR, BI SOLUTION VECTOR, X (REAL AND IMAG. PART). */
  /* ----------------------------------------------------------------------- */
  /* Parameter adjustments */
  const int n = ar.rows();

  /* Function Body */
  md = ml + mu + 1;
  md1 = md + 1;
  mdm = md - 1;
  nm1 = n - 1;
  if (ml == 0) {
    goto L25;
  }
  if (n == 1) {
    goto L50;
  }
  i1 = nm1;
  for (k = 1; k <= i1; ++k) {
    m = ip(k - 1);
    tr = br(m - 1);
    ti = bi(m - 1);
    br(m - 1) = br(k - 1);
    bi(m - 1) = bi(k - 1);
    br(k - 1) = tr;
    bi(k - 1) = ti;
    /* Computing MIN */
    i2 = ml, i3 = n - k;
    mdl = std::min(i2, i3) + md;
    i2 = mdl;
    for (i = md1; i <= i2; ++i) {
      imd = i + k - md;
      prodr = ar(i - 1, k - 1) * tr - ai(i - 1, k - 1) * ti;
      prodi = ai(i - 1, k - 1) * tr + ar(i - 1, k - 1) * ti;
      br(imd - 1) += prodr;
      bi(imd - 1) += prodi;
      /* L10: */
    }
    /* L20: */
  }
L25:
  i1 = nm1;
  for (kb = 1; kb <= i1; ++kb) {
    k = n + 1 - kb;
    den = ar(md - 1, k - 1) * ar(md - 1, k - 1) +
          ai(md - 1, k - 1) * ai(md - 1, k - 1);
    prodr = br(k - 1) * ar(md - 1, k - 1) + bi(k - 1) * ai(md - 1, k - 1);
    prodi = bi(k - 1) * ar(md - 1, k - 1) - br(k - 1) * ai(md - 1, k - 1);
    br(k - 1) = prodr / den;
    bi(k - 1) = prodi / den;
    tr = -br(k - 1);
    ti = -bi(k - 1);
    kmd = md - k;
    /* Computing MAX */
    i2 = 1, i3 = kmd + 1;
    lm = std::max(i2, i3);
    i2 = mdm;
    for (i = lm; i <= i2; ++i) {
      imd = i - kmd;
      prodr = ar(i - 1, k - 1) * tr - ai(i - 1, k - 1) * ti;
      prodi = ai(i - 1, k - 1) * tr + ar(i - 1, k - 1) * ti;
      br(imd - 1) += prodr;
      bi(imd - 1) += prodi;
    }
  }
  den = ar(md - 1, 0) * ar(md - 1, 0) + ai(md - 1, 0) * ai(md - 1, 0);
  prodr = br(0) * ar(md - 1, 0) + bi(0) * ai(md - 1, 0);
  prodi = bi(0) * ar(md - 1, 0) - br(0) * ai(md - 1, 0);
  br(0) = prodr / den;
  bi(0) = prodi / den;
L50:;
}

template <typename T, int SizeN = Eigen::Dynamic, int SizeM = Eigen::Dynamic>
void elmhes(const int low, const int igh, Matrix<T, SizeN, SizeM> &a,
            Vector<int, SizeN> &inter) {
  /* System generated locals */
  int a_dim1, a_offset, i1, i2, i3;
  T d1;

  /* Local variables */
  int i, j, m;
  T x, y;
  int la, mm1, kp1, mp1;

  /*     this subroutine is a translation of the algol procedure elmhes, */
  /*     num. math. 12, 349-368(1968) by martin and wilkinson. */
  /*     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971). */

  /*     given a real general matrix, this subroutine */
  /*     reduces a submatrix situated in rows and columns */
  /*     low through igh to upper hessenberg form by */
  /*     stabilized elementary similarity transformations. */

  /*     on input: */

  /*      nm must be set to the row dimension of two-dimensional */
  /*        array parameters as declared in the calling program */
  /*        dimension statement; */

  /*      n is the order of the matrix; */

  /*      low and igh are ints determined by the balancing */
  /*        subroutine  balanc.      if  balanc  has not been used, */
  /*        set low=1, igh=n; */

  /*      a contains the input matrix. */

  /*     on output: */

  /*      a contains the hessenberg matrix.  the multipliers */
  /*        which were used in the reduction are stored in the */
  /*        remaining triangle under the hessenberg matrix; */
  /*      int contains information on the rows and columns */
  /*        interchanged in the reduction. */
  /*        only elements low through igh are used. */
  /*     questions and comments should be directed to b. s. garbow, */
  /*     applied mathematics division, argonne national laboratory */
  /*     ------------------------------------------------------------------ */

  const int nm = a.rows();
  const int n = a.cols();

  la = igh - 1;
  kp1 = low + 1;
  if (la < kp1) {
    goto L200;
  }

  i1 = la;
  for (m = kp1; m <= i1; ++m) {
    mm1 = m - 1;
    x = 0.;
    i = m;

    i2 = igh;
    for (j = m; j <= i2; ++j) {
      if ((d1 = a(j - 1, mm1 - 1), abs(d1)) <= abs(x)) {
        goto L100;
      }
      x = a(j - 1, mm1 - 1);
      i = j;
    L100:;
    }

    inter[m] = i;
    if (i == m) {
      goto L130;
    }
    i2 = n;
    for (j = mm1; j <= i2; ++j) {
      y = a(i - 1, j - 1);
      a(i - 1, j - 1) = a(m - 1, j - 1);
      a(m - 1, j - 1) = y;
    }

    i2 = igh;
    for (j = 1; j <= i2; ++j) {
      y = a(j - 1, i - 1);
      a(j - 1, i - 1) = a(j - 1, m - 1);
      a(j - 1, m - 1) = y;
    }
  L130:
    if (x == 0.) {
      goto L180;
    }
    mp1 = m + 1;

    i2 = igh;
    for (i = mp1; i <= i2; ++i) {
      y = a(i - 1, mm1 - 1);
      if (y == 0.) {
        goto L160;
      }
      y /= x;
      a(i - 1, mm1 - 1) = y;
      i3 = n;
      for (j = m; j <= i3; ++j) {
        a(i - 1, j - 1) -= y * a(m - 1, j - 1);
      }
      i3 = igh;
      for (j = 1; j <= i3; ++j) {
        a(j - 1, m - 1) += y * a(j - 1, i - 1);
      }
    L160:;
    }
  L180:;
  }
L200:;
}
} // namespace stiff

#endif // STIFF_DECSOL_HPP

/*************************************************************************
 *  Compilation:  javac FFT.java
 *  Execution:    java FFT N
 *  Dependencies: Complex.java
 *
 *  Compute the FFT and inverse FFT of a length N complex sequence.
 *  Bare bones implementation that runs in O(N log N) time. Our goal
 *  is to optimize the clarity of the code, rather than performance.
 *
 *  Limitations
 *  -----------
 *   -  assumes N is a power of 2
 *
 *   -  not the most memory efficient algorithm (because it uses
 *      an object type for representing complex numbers and because
 *      it re-allocates memory for the subarray, instead of doing
 *      in-place or reusing a single temporary array)
 *
 *************************************************************************/

public class RFFT extends FFT {

    // compute the FFT of x[], assuming its length is a power of 2
    public static Complex[] fft(Complex[] x) {
        int N = x.length;

        // base case
        if (N == 1) return new Complex[] { x[0] };

        // radix 2 Cooley-Tukey FFT
        if (N % 2 != 0) { throw new RuntimeException("N is not a power of 2"); }

        // fft of even terms
        Complex[] even = new Complex[N/2];
        for (int k = 0; k < N/2; k++) {
            even[k] = x[2*k];
        }
        Complex[] q = fft(even);

        // fft of odd terms
        Complex[] odd  = even;  // reuse the array
        for (int k = 0; k < N/2; k++) {
            odd[k] = x[2*k + 1];
        }
        Complex[] r = fft(odd);

        // combine
        Complex[] y = new Complex[N];
        for (int k = 0; k < N/2; k++) {
            double kth = -2 * k * Math.PI / N;
            Complex wk = new Complex(Math.cos(kth), Math.sin(kth));
            y[k]       = q[k].plus(wk.times(r[k]));
            y[k + N/2] = q[k].minus(wk.times(r[k]));
        }
        return y;
    }


    // compute the inverse FFT of x[], assuming its length is a power of 2
    public static Complex[] ifft(Complex[] x) {
        int N = x.length;
        Complex[] y = new Complex[N];

        // take conjugate
        for (int i = 0; i < N; i++) {
            y[i] = x[i].conjugate();
        }

        // compute forward FFT
        y = fft(y);

        // take conjugate again
        for (int i = 0; i < N; i++) {
            y[i] = y[i].conjugate();
        }

        // divide by N
        for (int i = 0; i < N; i++) {
            y[i] = y[i].times(1.0 / N);
        }

        return y;

    }

    // compute the circular convolution of x and y
    public static Complex[] cconvolve(Complex[] x, Complex[] y) {

        // should probably pad x and y with 0s so that they have same length
        // and are powers of 2
        if (x.length != y.length) { throw new RuntimeException("Dimensions don't agree"); }

        int N = x.length;

        // compute FFT of each sequence
        Complex[] a = fft(x);
        Complex[] b = fft(y);

        // point-wise multiply
        Complex[] c = new Complex[N];
        for (int i = 0; i < N; i++) {
            c[i] = a[i].times(b[i]);
        }

        // compute inverse FFT
        return ifft(c);
    }


    // compute the linear convolution of x and y
    public static Complex[] convolve(Complex[] x, Complex[] y) {
        Complex ZERO = new Complex(0, 0);

        Complex[] a = new Complex[2*x.length];
        for (int i = 0;        i <   x.length; i++) a[i] = x[i];
        for (int i = x.length; i < 2*x.length; i++) a[i] = ZERO;

        Complex[] b = new Complex[2*y.length];
        for (int i = 0;        i <   y.length; i++) b[i] = y[i];
        for (int i = y.length; i < 2*y.length; i++) b[i] = ZERO;

        return cconvolve(a, b);
    }

    // display an array of Complex numbers to standard output
    public static void show(Complex[] x, String title) {
        System.out.println(title);
        System.out.println("-------------------");
        for (int i = 0; i < x.length; i++) {
            System.out.println(x[i]);
        }
        System.out.println();
    }


    public static Complex[] YW_to_ComplexYW(double[] SZ)
    {
        int count = SZ.length;
        Complex[] C_SZ = new Complex[count];
        for (int i = 0; i < count; i++)
        {
            Complex d = new Complex(SZ[i],0);
            C_SZ[i] = d;
        }
        return C_SZ;
    }

    public static Complex[] rfft(float[] r) {
        int N = r.length;
        Complex[] x = toComplex(r);
        Complex[] y = fft(x);
        return y;
    }

    public static Complex[] toComplex(float[] r){
        int N = r.length;
        Complex[] x = new Complex[N];
        for (int i = 0; i < N; i++) {
            x[i] = new Complex(r[i], 0);
        }
        return x;
    }
    public static float[] toFloat_bak(Complex[] c){
        int N = c.length*2;
        float[] x = new float[N];
        int k=0;
        for (int i = 0; i < c.length; i++) {
            x[k++] = (float)c[i].re();
            x[k++] = (float)c[i].im();
        }
        return x;
    }
    public static float[] toFloat(Complex[] c){
        int N = c.length;
        float[] x = new float[N];
        int k=0;
        for (int i = 0; i < c.length; i++) {
            x[i] = (float)c[i].re();
        }
        return x;
    }

    public static float[] conv(float[] x, float[] y) {
        int Ly = x.length + y.length - 1;
        int Ly2 = pow2(nextpow2(Ly));
        x = patch(x, Ly2);
        y = patch(y, Ly2);
        Complex[] dat = toComplex(x);
        Complex[] h = toComplex(y);
        Complex[] c = convolve(dat,h);
        float[] rc = toFloat(c);
        rc = cut(rc, Ly);
        return rc;
    }

    public static float[] conv_bak(float[] x, float[] y) {
        int Ly = x.length + y.length-1;
        int Ly2 = pow2(nextpow2(Ly));
        x = patch(x, Ly2);
        y = patch(y, Ly2);
        Complex[] X = rfft(x);
        Complex[] H = rfft(y);
        Complex[] Y = dotMult(X, H);
        Complex[] yy = ifft(Y);
        float[] z = real(yy);
        //z = wkeep(z, Ly);
        //z = cut(z, Ly);
        return z;
    }
    //返回序列中间部分--见Matlab wkeep

    public static float[] wkeep(float[] d, int len) {
        if (len >= d.length) {
            return d;
        }
        int dif = d.length - len,
                start = dif / 2;
        if ((start == 0) && (dif > 0)) {
            start = 1;
        }
        float[] keep = new float[len];
        for (int i = 0; i < len; i++) {
            keep[i] = d[i + start];
        }
        return keep;
    }

    public static int pow2(int p) {
        return (int) Math.pow(2, p);
    }

    public static int nextpow2(int n) {
        double log2 = (Math.log(n) / Math.log(2));
        if ((int) log2 < log2) {
            log2 += 1;
        }
        return (int) log2;
    }

    public static float[] real(Complex[] x) {
        int N = x.length;
        float[] r = new float[N];
        for (int i = 0; i < r.length; i++) {
            r[i] = (float) x[i].re();
        }
        return r;
    }

    public static Complex[] dotMult(Complex[] x, Complex[] h) {
        Complex[] xh = new Complex[x.length];
        for (int i = 0; i < xh.length; i++) {
            xh[i] = x[i].times(h[i]);
        }
        return xh;
    }

    public static float[] patch(float[] x, int len) {
        float[] xx = new float[len];
        for (int i = 0; i < x.length; i++) {
            xx[i] = x[i];
        }
        return xx;
    }

    public static float[] cut(float[] x, int len) {
        if (x.length > len) {
            float[] xx = new float[len];
            for (int i = 0; i < len; i++) {
                xx[i] = x[i];
            }
            return xx;
        } else {
            return x;
        }
    }

    // display an array of Complex numbers to standard output
    public static void show(float[] x, String title) {
        System.out.println(title);
        System.out.println("-------------------");
        for (int i = 0; i < x.length; i++) {
            System.out.println(x[i]);
        }
        System.out.println();
    }

    public static void main(String[] args) {
//        int N = 16;
//        float[] x = new float[N];
//        for (int i = 0; i < N; i++) {
//            x[i] = (float) (-2 * Math.random() + 1.0);
//        }
		float[] x ={94, 94, 112, 112, 112, 94, 94, 94, 112, 112, 112, 112, 112, 112, 112, 112, 112, 112, 130, 112, 112, 112, 130, 130, 130, 130, 130, 130, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 166, 148, 148, 148, 166, 166, 166, 148, 148, 148, 166, 166, 166, 166, 166, 166, 166, 166, 184, 166, 166, 166, 166, 166, 166, 166, 184, 166, 166, 166, 184, 184, 184, 184, 184, 166, 166, 166, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 166, 166, 184, 220, 184, 166, 184, 202, 184, 166, 166, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 166, 166, 166, 184, 166, 166, 166, 166, 166, 166, 166, 166, 148, 148, 148, 166, 166, 166, 166, 166, 148, 148, 148, 166, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 130, 130, 130, 148, 148, 148, 148, 148, 130, 130, 130, 148, 130, 130, 130, 130, 130, 130, 112, 112, 112, 112, 112, 112, 112, 112, 112, 130, 112, 112, 112, 112, 94, 94, 94, 112, 112, 112, 112, 112, 112, 112, 112, 112, 94, 94, 94, 112, 94, 94, 94, 94, 94, 112, 94, 94, 94, 94, 94, 112, 94, 94, 94, 112, 94, 94, 94, 94, 94, 94, 94, 94, 76, 76, 76, 76, 76, 76, 58, 58, 58, 76, 58, 58, 58, 58, 58, 58, 40, 40, 40, 40, 40, 40, 40, 58, 40, 40, 40, 58, 58, 76, 40, 22, 40, 76, 58, 40, 40, 40, 40, 58, 40, 40, 40, 58, 40, 40, 40, 58, 40, 40, 40, 58, 58, 58, 40, 40, 40, 58, 40, 40, 40, 58, 58, 58, 58, 58, 58, 58, 58, 76, 58, 58, 58, 76, 58, 58, 58, 76, 76, 76, 76, 76, 58, 58, 58, 76, 58, 58, 58, 76, 76, 76, 76, 76, 76, 76, 76, 76, 76, 76, 76, 76, 76, 76, 76, 76, 76, 76, 76, 76, 76, 76, 76, 76, 76, 76, 76, 94, 76, 76, 76, 76, 76, 76, 76, 94, 94, 94, 94, 94, 76, 76, 76, 94, 94, 94, 94, 94, 94, 94, 94, 94, 94, 94, 94, 112, 112, 112, 112, 112, 112, 112, 112, 112, 112, 130, 112, 112, 112, 130, 112, 112, 112, 130, 130, 130, 130, 148, 130, 130, 130, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 166, 202, 166, 148, 166, 184, 166, 166, 166, 184, 166, 166, 166, 166, 166, 184, 166, 166, 166, 184, 184, 184, 166, 166, 166, 166, 166, 184, 184, 184, 166, 166, 166, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 166, 166, 166, 184, 166, 166, 166, 166, 166, 184, 166, 166, 166, 184, 166, 166, 166, 184, 166, 166, 166, 184, 184, 184, 184, 184, 184, 184, 184, 184, 166, 166, 166, 184, 184, 184, 184, 184, 166, 166, 166, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 166, 166, 166, 184, 166, 166, 166, 166, 166, 166, 166, 166, 166, 166, 148, 148, 148, 166, 166, 166, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 130, 130, 130, 148, 130, 112, 112, 130, 130, 130, 130, 130, 112, 112, 112, 112, 112, 112, 112, 112, 112, 112, 112, 112, 94, 94, 94, 94, 94, 94, 94, 94, 94, 94, 94, 94, 76, 76, 76, 94, 94, 94, 94, 94, 94, 94, 76, 76, 76, 94, 76, 76, 76, 76, 76, 76, 76, 76, 76, 76, 76, 76, 58, 58, 58, 76, 76, 76, 76, 76, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 76, 58, 58, 58, 76, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 40, 40, 40, 58, 58, 58, 58, 58, 58, 58, 40, 40, 40, 58, 40, 40, 58, 76, 58, 58, 58, 58, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 76, 76, 76, 58, 58, 58, 76, 58, 58, 58, 76, 58, 58, 58, 58, 58, 76, 76, 94, 76, 58, 76, 94, 76, 76, 94, 112, 112, 112, 94, 94, 94, 112, 94, 94, 94, 112, 94, 94, 112, 130, 112, 112, 94, 94, 94, 94, 94, 94, 94, 112, 112, 112, 112, 112, 112, 130, 112, 112, 112, 130, 112, 112, 112, 130, 112, 112, 112, 112, 112, 112, 112, 112, 112, 112, 112, 112, 112, 112, 112, 112, 112, 112, 112, 112, 112, 112, 112, 112, 112, 112, 112, 130, 112, 112, 112, 112, 112, 112, 112, 112, 112, 112, 112, 130, 130, 130, 130, 130, 112, 112, 112, 130, 112, 112, 112, 130, 130, 130, 130, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 166, 166, 166, 148, 148, 148, 166, 166, 166, 166, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 202, 184, 184, 184, 184, 184, 202, 202, 202, 184, 184, 184, 202, 184, 184, 184, 202, 202, 202, 202, 220, 202, 184, 184, 202, 202, 202, 202, 202, 202, 202, 202, 202, 202, 202, 202, 202, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 166, 166, 166, 166, 166, 166, 166, 166, 166, 166, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 130, 130, 130, 148, 130, 130, 130, 148, 130, 130, 130, 148, 130, 130, 130, 148, 130, 130, 130};
    	//Complex[] x = YW_to_ComplexYW(realInA);
        show(x, "x");
        // FFT of original data
        Complex[] y = rfft(x);
        show(y, "y = fft(x)");
        // take inverse FFT
        Complex[] z = ifft(y);
        show(z, "z = irfft(y)");

        float[] a = {4, 5, 6, 7};
        float[] b = {7, 6, 5, 4};
        float[] c = conv(a, b);
        show(c, "c = conv(a,b)");

    }
}

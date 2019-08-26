
import std;

enum M = 1000;
enum N = 1000;

void main()
{
    size_t total_bits = 0;
    size_t error_bits = 0;
    foreach(_; 0 .. 100) {
        writefln!"iter: %s"(_);
        byte[] bits;
        foreach(i; 0 .. M) {
            auto p = uniform01();
            bits ~= p > 0.5 ? 1 : -1;
        }

        Complex!double[][] h = new Complex!double[][](M, N);
        foreach(i; 0 .. M)
            foreach(j; 0 .. N) {
                h[i][j] = complexGaussian01!double(rndGen);
            }


        Complex!double[] sum = new Complex!double[](N);
        sum[] = Complex!double(0);
        foreach(i; 0 .. M) {
            foreach(j; 0 .. N)
                sum[j] += h[i][j] * bits[i];
        }

        auto bpresult = gaussian_bp(M, N, 0, sum, h);
        foreach(i; 0 .. M) {
            total_bits += 1;
            error_bits += bits[i] != bpresult[i] ? 1 : 0;
        }

    }

    writefln("Total: %s bits", total_bits);
    writefln("Error: %s bits", error_bits);
    writefln("BER: %s", error_bits * 1.0L / total_bits);
}


byte[] gaussian_bp(size_t m, size_t n, double N0, Complex!double[] y, Complex!double[][] h)
{
    double[] softbit = new double[m];
    Complex!double[] yic = new Complex!double[n];
    double[] sigma = new double[n];

    softbit[] = 0;

    foreach(_; 0 .. 10) {
        yic[] = y[];
        sigma[] = N0;
        foreach(j; 0 .. n) {
            foreach(k; 0 .. m) {
                yic[j] -= h[k][j] * softbit[k];
                sigma[j] += h[k][j].sqAbs * (1 - softbit[k])^^2;
            }
        }

        foreach(i; 0 .. m) {
            double new_llr = 0;
            foreach(j; 0 .. n) {
                Complex!double y_ji = yic[j] + h[i][j] * softbit[i];
                double sig_ji = sigma[j] - h[i][j].sqAbs * (1 - softbit[i])^^2;
                new_llr += 4 * (y_ji * h[i][j].conj).re / sig_ji;
            }

            new_llr = min(max(new_llr, -4), +4);
            softbit[i] = tanh(new_llr/2);
        }
    }

    byte[] dst = new byte[m];
    foreach(i; 0 .. m)
        dst[i] = softbit[i] > 0 ? 1 : -1;

    return dst;
}


Complex!F complexGaussian01(F = real, Rnd)(ref Rnd rnd)
{
    import std.complex : expi;

    F x = uniform01(rnd),
      y = uniform01(rnd);

    typeof(return) dst = sqrt(-1 * log(x)) * expi(2 * PI * y);

    return dst;
}

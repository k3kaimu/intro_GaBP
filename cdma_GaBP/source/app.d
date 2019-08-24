import dffdd;

import std;

enum M = 200;

void main()
{
    size_t total_bits = 0;
    size_t error_bits = 0;
    foreach(_; 0 .. 1024) {
        writefln!"iter: %s"(_);
        byte[] bits;
        foreach(i; 0 .. M) {
            auto p = uniform01();
            bits ~= p > 0.5 ? 1 : -1;
        }

        double[][] h = new double[][](M, L1CACode.codeLength);
        double[] sum = new double[](L1CACode.codeLength);
        sum[] = 0;
        foreach(prn; 0 .. M) {
            auto s = gen_l1ca(prn+1);
            sum[] += s[] * double(bits[prn]);
            h[prn] = s;
        }

        Complex!double[] recv = sum.map!(a => Complex!double(a)).array();
        Complex!double[][] hmat = h.map!(as => as.map!(a => Complex!double(a)).array()).array();

        auto bpresult = gaussian_bp(h.length, sum.length, 0, recv, hmat);
        foreach(i; 0 .. M) {
            total_bits += 1;
            error_bits += bits[i] != bpresult[i] ? 1 : 0;
        }

    }

    writefln("Total: %s bits", total_bits);
    writefln("Error: %s bits", error_bits);
    writefln("BER: %s", error_bits * 1.0L / total_bits);
}


double[] gen_l1ca(uint prn)
{
    auto code = L1CACode(prn);
    double[] signal = new double[code.codeLength];

    code.read(signal);

    return signal;
}


byte[] gaussian_bp(size_t m, size_t n, double N0, Complex!double[] y, Complex!double[][] h)
{
    double[] softbit = new double[m];
    Complex!double[] yic = new Complex!double[n];
    double[] sigma = new double[n];

    softbit[] = 0;

    foreach(_; 0 .. 20) {
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
                new_llr += (-(y_ji - h[i][j]).sqAbs + (y_ji + h[i][j]).sqAbs) / sig_ji;
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

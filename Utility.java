
import java.util.Random;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author hoangcuong2011
 */
public class Utility {    
    public static double sampleGamma(double k, double theta) {
        Random rng = new Random();
        boolean accept = false;
        if (k < 1) {
            // Weibull algorithm
            double c = (1 / k);
            double d = ((1 - k) * Math.pow(k, (k / (1 - k))));
            double u, v, z, e, x;
            do {
                u = rng.nextDouble();
                v = rng.nextDouble();
                z = -Math.log(u);
                e = -Math.log(v);
                x = Math.pow(z, c);
                if ((z + e) >= (d + x)) {
                    accept = true;
                }
            } while (!accept);
            return (x * theta);
        } else {
            // Cheng's algorithm
            double b = (k - Math.log(4));
            double c = (k + Math.sqrt(2 * k - 1));
            double lam = Math.sqrt(2 * k - 1);
            double cheng = (1 + Math.log(4.5));
            double u, v, x, y, z, r;
            do {
                u = rng.nextDouble();
                v = rng.nextDouble();
                y = ((1 / lam) * Math.log(v / (1 - v)));
                x = (k * Math.exp(y));
                z = (u * v * v);
                r = (b + (c * y) - x);
                if ((r >= ((4.5 * z) - cheng))
                        || (r >= Math.log(z))) {
                    accept = true;
                }
            } while (!accept);
            return (x * theta);
        }
    }
}

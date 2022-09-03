#import "Basic";
#import,file "module.jai"(2,1,1);

main :: () {
    println("% : %", fmt(e12), formatInt(e12.basis, 2));
    println("% : %", fmt(e21), formatInt(e21.basis, 2));

    println("e0e0 %", fmt(e0*e0));
    println("e1e1 %", fmt(e1*e1));
    println("e1e2 %", fmt(e1*e2));
    println("e2e1 %", fmt(e2*e1));

    println("hamming_weight(5) = %", hamming_weight(5));

    println("42.0 * e2 = %", fmt(42.0 * e2));
    println("e3 * 313.0 = %", fmt(e3 * 313.0));

    x := mv(1.0, -42.0*e2, e3*313.0, e12 / 10.0);
    println("x = % = %", fmt(x), x);

    x = x + x;
    println("x = %", fmt(x));

    x = x + 100.0*e2 - 1.0*e23;
    println("x = %", fmt(x));
    println("-x = %", fmt(-x));

    y := 3.0 - 4.0 * e1 + 0.123 * e13 + 0.321 * e234;
    println("e1   = % = %", e1  , fmt(e1  ));
    println("e13  = % = %", e13 , fmt(e13 ));
    println("e234 = % = %", e234, fmt(e234));
    println("y = % = %", fmt(y), y);
    println("x*y = %", fmt(x*y));
    println("y*x = %", fmt(y*x));
    println("x*x = %", fmt(x*x));
    println("y*y = %", fmt(y*y));


    println("N_ELEMENTS = %", N_ELEMENTS);
}
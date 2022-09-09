#import "Basic";
#import,file "module.jai"(3,0,0);

main :: () {
    println("% : %", fmt(e12), formatInt(e12.basis, 2));
    println("% : %", fmt(e21), formatInt(e21.basis, 2));

    println("e00 %", fmt(e0*e0));
    println("e11 %", fmt(e1*e1));
    println("e12 %", fmt(e1*e2));
    println("e21 %", fmt(e2*e1));

    // // println("hamming_weight(5) = %", hamming_weight(5));

    println("42.0 * e2 = %", fmt(42.0 * e2));
    println("e3 * 313.0 = %", fmt(e3 * 313.0));


    x := MultiVector(e1.basis, e2.basis).{.[1.0, 2.0]};
    println("x : % = %", type_of(x), x);

    y := e1 + e1;
    println("y : % = %", type_of(y), y);

    z := e1 + e2;
    println("z : % = %", type_of(z), z);

    w := 1.0 + e1;
    println("w : % = %", type_of(w), w);

    v := y + w + y;
    println("v : % = %", type_of(v), v);

    v2 := v + z;
    println("v2 : % = %", type_of(v2), v2);

    v3 := v2 - x - 42.0 + e123;
    println("v3 : % = %", type_of(v3), v3);




    // println("% + % = %", e1, e1, e1+e1);
    // x := mv(1.0, -42.0*e2, e3*313.0, e12 / 10.0);
    // println("x = % = %", fmt(x), x);

    // x = x + x;
    // println("x = %", fmt(x));

    // x = x + 100.0*e2 - 1.0*e23;
    // println("x = %", fmt(x));
    // println("-x = %", fmt(-x));

    // y := 3.0 - 4.0 * e1 + 0.123 * e13 + 0.321 * e234;
    // println("e1   = % = %", e1  , fmt(e1  ));
    // println("e13  = % = %", e13 , fmt(e13 ));
    // println("e234 = % = %", e234, fmt(e234));
    // println("y = % = %", fmt(y), y);
    // println("x*y = %", fmt(x*y));
    // println("y*x = %", fmt(y*x));
    // println("x*x = %", fmt(x*x));
    // println("y*y = %", fmt(y*y));


    // println("N_ELEMENTS = %", N_ELEMENTS);
}
#import "Basic";
#import,file "module.jai"(1,2,0);

main :: () {
    println("e00 %", fmt(e0*e0));
    println("e11 %", fmt(e1*e1));
    println("e12 %", fmt(e1*e2));
    println("e21 %", fmt(e2*e1));

    // println("hamming_weight(5) = %", hamming_weight(5));

    println("42.0 * e2 = %", fmt(42.0 * e2));
    println("e3 * 313.0 = %", fmt(e3 * 313.0));


    x := MultiVector(1,2).{.[1.0, 2.0]};


    mv_print :: (name: string, mv: $V/MultiVector) {
        println("% : % = %", name, type_of(mv), fmt(mv));
    }

    println("x : % = %", type_of(x), fmt(x));

    y := e1 + e1;
    println("y : % = %", type_of(y), fmt(y));

    z := e1 + e2;
    println("z : % = %", type_of(z), fmt(z));

    w := 1.0 + e1;
    println("w : % = %", type_of(w), fmt(w));

    v := y + w + y;
    println("v : % = %", type_of(v), fmt(v));

    v2 := v + z;
    println("v2 : % = %", type_of(v2), fmt(v2));

    v3 := v2 - x - 42.0 + e123;
    println("v3 : % = %", type_of(v3), fmt(v3));

    v4 := (1.0 + e123) + v3;
    println("v4 : % = %", type_of(v4), fmt(v4));

    k0 := v4 * v3;
    mv_print("k0", k0);

    k1 := v3 * w;
    mv_print("k1", k1);

    k2 := MV_ZERO;
    k2 += k1;
    for i: 0..3 {
        mv_print("k2", k2);
        k2 *= k1;
    }

    v3_dual := dual(v3);
    mv_print("v3_dual", v3_dual);

    v3_dual_reverse := reverse(v3_dual);
    mv_print("v3_dual_reverse", v3_dual_reverse);

    println("select 7 of v3_dual_reverse: %", select(v3_dual_reverse, 7));
    println("select 2 of v3_dual_reverse: %", select(v3_dual_reverse, 2));

    // mv_print("k2", k1 * k1);
    // mv_print("k2", k1 * k1 * k1);
    // mv_print("k2", k1 * k1 * k1 * k1);
    // mv_print("k2", k1 * k1 * k1 * k1 * k1);
}
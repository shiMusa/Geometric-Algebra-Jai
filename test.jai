#import "Basic";
#import,file "module.jai"(2,0,1,float32,"hodge");

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

    v3_dual := dual_polar(v3);
    mv_print("v3_dual", v3_dual);

    v3_dual_reverse := reverse(v3_dual);
    mv_print("v3_dual_reverse", v3_dual_reverse);

    println("select 0 of v3_dual_reverse: %", fmt(select(v3_dual_reverse, 0)));
    println("select 1 of v3_dual_reverse: %", fmt(select(v3_dual_reverse, 1)));
    println("select 2 of v3_dual_reverse: %", fmt(select(v3_dual_reverse, 2)));
    println("select 3 of v3_dual_reverse: %", fmt(select(v3_dual_reverse, 3)));

    println("polar of v3: %", fmt(dual_polar(v3)));

    println("v3 * k0 = %", fmt(v3*k0));
    println("meet of v3 and k0: %", fmt(meet(v3, k0)));
    println("inner of v3 and k0: %", fmt(inner(v3, k0)));
    println("comm of v3 and k0: %", fmt(comm(v3, k0)));

    println("x * w = %", fmt(x*w));
    println("meet of x and w: %", fmt(meet(x, w)));
    println("join of x and w: %", fmt(join(x, w)));
    println("inner of x and w: %", fmt(inner(x, w)));
    println("comm of x and w: %", fmt(comm(x, w)));

    println("x * e3 = %", fmt(x*e3));
    println("meet of x and e3: %", fmt(meet(x, e3)));
    println("join of x and e3: %", fmt(join(x, e3)));
    println("inner of x and e3: %", fmt(inner(x, e3)));
    println("comm of x and e3: %", fmt(comm(x, e3)));

    println("v : % = %", type_of(v), fmt(v));
    println("w : % = %", type_of(w), fmt(w));
    println("v * w = %", fmt(v*w));
    println("meet of w and v: %", fmt(meet(w, v)));
    println("join of w and v: %", fmt(join(w, v)));
    println("inner of w and v: %", fmt(inner(w, v)));
    println("comm of w and v: %", fmt(comm(w, v)));

    println("e1 * e2 : % = %, e2 * e1 : % = %", type_of(e1*e2), fmt(e1*e2), type_of(e2*e1), fmt(e2*e1));
    println("e1 * e2 * e1 : % = %, e2 * e1 * e1 : % = %", type_of(e1*e2*e1), fmt(e1*e2*e1), type_of(e2*e1*e1), fmt(e2*e1*e1));

    #if POS == 2 && NEG == 0 && ZER == 1 {
        println("2D PGA");

        l := join(point(-1.0, 0.0), point(1.0, 0.0), .HODGE);
        p := point(0.5, 1.0);
        perp := orthogonal_line(l, p);
        intersection := meet(l, perp);
        println("orthogonal line to line % through point % = %", fmt_euclid(l), fmt_euclid(p), fmt_euclid(perp));
        println("both lines intersect at %", fmt_euclid(intersection));

        p2 := point(10.0, -20.0);
        println("project % onto % : %", fmt_euclid(p2), fmt_euclid(l), fmt_euclid(project(p2, l)));
        println("project % onto % : %", fmt_euclid(p2), fmt_euclid(perp), fmt_euclid(project(p2, perp)));
        
        println("project % onto % : % < %", fmt_euclid(l), fmt_euclid(p2), fmt_euclid(project(l, p2)), type_of(project(l, p2)));
        println("orthogonal direction: %", fmt_euclid(orthogonal_direction(project(l, p2))));
        println("project % onto % : % < %", fmt_euclid(perp), fmt_euclid(p2), fmt_euclid(project(perp, p2)), type_of(project(perp, p2)));
        println("orthogonal direction: %", fmt_euclid(orthogonal_direction(project(perp, p2))));


        println("distance between % and %: %", fmt_euclid(p), fmt_euclid(p2), distance(p, p2));
        l2 := join(point(3.0, 42.0), point(-10.1, 42.0), .HODGE);
        println("distance between % and % = %", fmt_euclid(l), fmt_euclid(l2), distance(l,l2));
    }
}
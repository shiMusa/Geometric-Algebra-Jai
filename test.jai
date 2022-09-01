#import "Basic";
#import,file "module.jai"(1,3,0);

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

    println("mv = %", mv(1.0, -42.0*e2, e3*313.0, e01 / 10.0));
    println("mv = %", fmt(mv(1.0, -42.0*e2, e3*313.0, e01 / 10.0)));
}
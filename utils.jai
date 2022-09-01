

powi :: (base: int, exponent: int) -> int {
    res := 1;
    for 0..exponent-1
        res *= base;
    return res;
}

// * Comparing to cheat sheet, the following notation changes need to be considered:
// * here  <>  cheat sheet
// * e3    <>  e0

Point     :: #run type_of(e23 + e31 + e12);
Direction :: #run type_of(e23 + e31);
Line      :: #run type_of(e1 + e2 + e3);

#run {
    println("2D-PGA:");
    println("Point     = %", Point);
    println("Direction = %", Direction);
    println("Line      = %", Line);
}

line :: (x: T, y: T, c: T = zero) -> Line {
    return x * e1 + y * e2 + c * e3;
}
x :: inline (l: Line) -> T {
    return l.values[0];
}
y :: inline (l: Line) -> T {
    return l.values[1];
}
c :: inline (l: Line) -> T {
    return l.values[2];
}
fmt_euclid :: (l: Line) -> string {
    return sprint("% x + % y + % = 0", x(l), y(l), c(l));
}

point :: (x: T, y: T) -> Point {
    return x * e23 + y * e13 + e12;
}
x :: inline (p: Point) -> T {
    return p.values[2]/p.values[0];
}
y :: inline (p: Point) -> T {
    return p.values[1]/p.values[0];
}
fmt_euclid :: (p: Point) -> string {
    return sprint("(%, %)", x(p), y(p));
}

direction :: (x: T, y: T) -> Direction {
    return x * e23 + y * e13;
}
x :: inline (d: Direction) -> T {
    return d.values[1];
}
y :: inline (d: Direction) -> T {
    return d.values[0];
}
fmt_euclid :: (d: Direction) -> string {
    return sprint("|%,%)", x(d), y(d));
}

orthogonal_line :: (to_line: Line, through_point: Point) -> Line {
    return inner(to_line, through_point);
}
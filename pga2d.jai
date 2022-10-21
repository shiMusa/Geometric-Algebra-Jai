
Point     :: #run type_of(e31 + e12 + e23);
Direction :: #run type_of(e31 + e12);
Line      :: #run type_of(e2 + e3 + e1);

#run {
    println("2D-PGA:");
    println("Point     = %", Point);
    println("Direction = %", Direction);
    println("Line      = %", Line);
}

line :: (x: T, y: T, c: T = zero) -> Line {
    return x * e2 + y * e3 + c * e1;
}

point :: (x: T, y: T) -> Point {
    return x * e31 + y * e12 + e23;
}

direction :: (x: T, y: T) -> Direction {
    return x * e31 + y * e12;
}

orthogonal_line :: (to: Line, through: Point) -> Line {
    return inner(to, through);
}
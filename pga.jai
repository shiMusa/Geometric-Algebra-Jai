
// * Comparing to cheat sheet, the following notation changes need to be considered:
// * here  <>  cheat sheet
// * e3    <>  e0

#run {
    println("%D-PGA:", POS);
}

// * 2D-PGA
#if POS == 2 {
    Point     :: #run type_of(e23 + e31 + e12);
    Direction :: #run type_of(e23 + e31);
    Line      :: #run type_of(join(point(0.0,1.0), point(1.0,0.0)));
    
    #run {
        println("Point     = %", Point);
        println("Direction = %", Direction);
        println("Line      = %", Line);
    }
}


line :: (x: T, y: T, c: T = zero) -> Line {
    res : Line  = ---;
    res.values[0] = x;
    res.values[1] = y;
    res.values[2] = c;
    return res;
    // return x * e1 + y * e2 + c * e3;
}
line :: inline (a: $_A/MultiVector) -> Line {
    return select(a, Line);
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
    res : Point = ---;
    res.values[0] = one;
    res.values[1] = y;
    res.values[2] = x;
    return res;
    // return x * e23 + y * e13 + e12;
}
point :: (a: $_A/MultiVector) -> Point {
    res := select(a, Point);
    // * normalize
    res.values[1] /= res.values[0];
    res.values[2] /= res.values[0];
    res.values[0] = one;
    return res;
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
direction :: inline (a: $_A/MultiVector) -> Direction {
    return select(a, Direction);
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


// * ##     ## ######## ######## ##     ##  #######  ########   ######
// * ###   ### ##          ##    ##     ## ##     ## ##     ## ##    ##
// * #### #### ##          ##    ##     ## ##     ## ##     ## ##
// * ## ### ## ######      ##    ######### ##     ## ##     ##  ######
// * ##     ## ##          ##    ##     ## ##     ## ##     ##       ##
// * ##     ## ##          ##    ##     ## ##     ## ##     ## ##    ##
// * ##     ## ########    ##    ##     ##  #######  ########   ######

orthogonal_line :: (to_line: Line, through_point: Point) -> Line {
    return dot(to_line, through_point);
}

project :: (pt: Point, onto_line: Line) -> Point {
    DEBUG :: false;

    #if DEBUG {
        println(      "project -------------------");
        defer println("^^^^^^^^^^^^^^^^^^^^^^^^^^^");
        tmp := dot(onto_line, pt) * onto_line;
        println("point = %", fmt(pt));
        println("onto_line = %", fmt(onto_line));
        println("dot = %", fmt(dot(onto_line, pt)));
        println("tmp = %", fmt(tmp));

        tmp2 := point(tmp);
        println("tmp2 = %", fmt(tmp2));
    }

    return point(dot(pt, onto_line) * onto_line);
}

project :: (ln: Line, onto_point: Point) -> Line {
    return line(dot(ln, onto_point) * onto_point);
}

orthogonal_direction :: (to_line: Line) -> Direction {
    return direction(dual_polar(to_line));
}

norm :: (l: Line) -> T {
    return sqrt(l.values[0]*l.values[0] + l.values[1]*l.values[1]);
}
norm :: (p: Point) -> T {
    return abs(p.values[0]);
}
ideal :: (p: Point) -> T {
    return abs(p.values[1]*p.values[1] + p.values[2]*p.values[2]);
}
ideal :: inline (l: Line) -> T {
    return l.values[2];
}
#if Line != MV_I {
    ideal :: inline (pseudo_scalar: MV_I) -> T {
        return pseudo_scalar.values[0];
    }
}


distance :: (a: Point, b: Point) -> T {
    return norm(join(a,b));
}
angle :: (a: Line, b: Line) -> T {
    an := a/norm(a);
    bn := b/norm(b);
    return arccos(dot(an,bn));
}
// * for parallel lines
distance :: (a: Line, b: Line) -> T {
    an := a/norm(a);
    bn := b/norm(b);
    return ideal(meet(an,bn));
}
#import "Basic";
#import "Sort";
#load "utils.jai";



Foo :: struct(T: Type) {
    val: T;
}
foo :: (val: $T) -> Foo(T) {
    res : Foo(T) = ---;
    res.val = val;
    return res;
}

operator + :: (a: Foo($T), b: Foo(T)) -> Foo(T) {
    return foo(a.val + b.val);
}
operator * :: (a: Foo($T), b: Foo(T)) -> Foo(T) {
    return foo(a.val * b.val);
}


// *  ######  ##    ## ##     ## ########   #######  ##       ####  ######
// * ##    ##  ##  ##  ###   ### ##     ## ##     ## ##        ##  ##    ##
// * ##         ####   #### #### ##     ## ##     ## ##        ##  ##
// *  ######     ##    ## ### ## ########  ##     ## ##        ##  ##
// *       ##    ##    ##     ## ##     ## ##     ## ##        ##  ##
// * ##    ##    ##    ##     ## ##     ## ##     ## ##        ##  ##    ##
// *  ######     ##    ##     ## ########   #######  ######## ####  ######

Types :: enum {
    VAR;
    OP;
}


Symbolic :: struct {
    types : [..]Types;
    extra : [..]int;

    // * Operators
    op_back : [..]int;
    op_type : [..]OpTypes;
    op_args : [..][2]ID;
}
restart :: () {
    using context.sym;

    array_reset(*types);
    array_reset(*extra);

    array_reset(*op_type);
    array_reset(*op_args);
    array_reset(*op_back);
}

#add_context sym: Symbolic;


ID :: struct {
    val: int;
}
id :: (val: int) -> ID {
    res : ID = ---;
    res.val = val;
    return res;
}
operator== :: inline (a: ID, b: ID) -> bool {
    return a.val == b.val;
}

var :: () -> ID {
    using context.sym;

    i := types.count;
    array_add(*types, .VAR);
    array_add(*extra, -1);
    return id(i);
}


OpTypes :: enum {
    ADD;
    MUL;
}


op :: (a: ID, b: ID, type: OpTypes) -> ID {
    using context.sym;
    
    // * check, if this operation already exists
    for oi, oii: op_args {
        if oi[0] == a && oi[1] == b && op_type[oii] == type {
            return id(op_back[oii]);
        }
    }

    // * new operation, so create
    tmp : [2]ID = ---;
    tmp[0] = a;
    tmp[1] = b;
    
    i := types.count;
    j := op_type.count;
    
    array_add(*types, .OP);
    array_add(*extra, j);

    array_add(*op_back, i);
    array_add(*op_type, type);
    array_add(*op_args, tmp);

    res := id(i);
    // println("op -> %", res);
    return res;
}


operator * :: (a: ID, b: ID) -> ID {
    return op(a, b, .MUL);
}
operator + :: (a: ID, b: ID) -> ID {
    return op(a, b, .ADD);
}



// *  ######  ######## ########
// * ##    ##    ##    ##     ##
// * ##          ##    ##     ##
// *  ######     ##    ########
// *       ##    ##    ##   ##
// * ##    ##    ##    ##    ##
// *  ######     ##    ##     ##


indent :: (num: int) -> string {
    builder : String_Builder;

    for i: 0..num-1 {
        append(*builder, " ");
    }

    return builder_to_string(*builder);
}

str :: (id: ID, shift := 0) -> string {
    using context.sym;
    
    ind := indent(shift);

    builder : String_Builder;

    print_to_builder(*builder, "%1%2(%3)", ind, types[id.val], id.val);

    if types[id.val] == .VAR
        return builder_to_string(*builder);

    // * operation then
    op_indx := extra[id.val];
    print_to_builder(*builder, " - %\n", op_type[op_indx]);
    
    arg_ids := op_args[op_indx];
    print_to_builder(*builder, "%\n", str(arg_ids[0], shift + 4));
    print_to_builder(*builder, "%"  , str(arg_ids[1], shift + 4));

    return builder_to_string(*builder);
}



// *  ######  #### ##     ## ########  ##       #### ######## ##    ##
// * ##    ##  ##  ###   ### ##     ## ##        ##  ##        ##  ##
// * ##        ##  #### #### ##     ## ##        ##  ##         ####
// *  ######   ##  ## ### ## ########  ##        ##  ######      ##
// *       ##  ##  ##     ## ##        ##        ##  ##          ##
// * ##    ##  ##  ##     ## ##        ##        ##  ##          ##
// *  ######  #### ##     ## ##        ######## #### ##          ##

simplify :: (id: ID) -> ID {
    using context.sym;
    
    if types[id.val] == .OP {
        i := extra[id.val];
        ids := op_args[i];

        // println("ids = %", ids);

        ids[0] = simplify(ids[0]);
        ids[1] = simplify(ids[1]);

        if op_type[i] == .ADD {

            if types[ids[0].val] == .OP && types[ids[1].val] == .OP {
                j0 := extra[ids[0].val];
                j1 := extra[ids[1].val];

                // println("j0 = %, j1 = %", j0, j1);

                if op_type[j0] == .MUL && op_type[j1] == .MUL {
                    ids0 := op_args[j0];
                    ids1 := op_args[j1];

                    // println("ids0 = %", ids0);
                    // println("ids1 = %", ids1);

                    // * a*b + a*c = a*(b+c)
                    if ids0[0] == ids1[0] return ids0[0] * (ids0[1] + ids1[1]);
                    if ids0[0] == ids1[1] return ids0[0] * (ids0[1] + ids1[0]);
                    if ids0[1] == ids1[0] return ids0[1] * (ids0[0] + ids1[1]);
                    if ids0[1] == ids1[1] return ids0[1] * (ids0[0] + ids1[0]);
                }
            }
        }
    }

    // * no simplification found
    return id;
}




// *  ######   #######  ##     ## ########  #### ##       ########
// * ##    ## ##     ## ###   ### ##     ##  ##  ##       ##
// * ##       ##     ## #### #### ##     ##  ##  ##       ##
// * ##       ##     ## ## ### ## ########   ##  ##       ######
// * ##       ##     ## ##     ## ##         ##  ##       ##
// * ##    ## ##     ## ##     ## ##         ##  ##       ##
// *  ######   #######  ##     ## ##        #### ######## ########

get_vars :: (id: ID) -> [..]int {
    using context.sym;
    
    res : [..]int;

    if types[id.val] == .VAR {
        array_add_if_unique(*res, id.val);
    }
    if types[id.val] == .OP {
        i := extra[id.val];
        ids := op_args[i];

        for get_vars(ids[0])
            array_add_if_unique(*res, it);

        for get_vars(ids[1])
            array_add_if_unique(*res, it);
    }


    cmp :: inline (a: int, b: int) -> s64 {
        return a - b;
    }
    bubble_sort(res, cmp);

    return res;
}

compile :: (id: ID) -> string {
    using context.sym;
    
    builder : String_Builder;

    if types[id.val] == .VAR
        print_to_builder(*builder, "v%", id.val);
    else if types[id.val] == .OP {
        i := extra[id.val];
        type := op_type[i];
        args := op_args[i];
        s0 := compile(args[0]);
        s1 := compile(args[1]);

        if type == .ADD
            print_to_builder(*builder, "% + %", s0, s1);
        else if type == .MUL
            print_to_builder(*builder, "(%) * (%)", s0, s1);
    }

    return builder_to_string(*builder);
}

make_function :: (name: string, id: ID) -> string {
    using context.sym;
    
    builder : String_Builder;

    print_to_builder(*builder, "% :: (", name);

    vars := get_vars(id);
    for vars {
        print_to_builder(*builder, "v%: float", it);
        if it_index < vars.count-1 {
            append(*builder, ", ");
        }
    }

    append(*builder, ") -> float {\n");

    print_to_builder(*builder, "    return %;\n}\n", compile(id));


    return builder_to_string(*builder);
}

// * ##     ##    ###    #### ##    ##
// * ###   ###   ## ##    ##  ###   ##
// * #### ####  ##   ##   ##  ####  ##
// * ## ### ## ##     ##  ##  ## ## ##
// * ##     ## #########  ##  ##  ####
// * ##     ## ##     ##  ##  ##   ###
// * ##     ## ##     ## #### ##    ##

main :: () {
    {
        x := foo(42.0);
        y := foo(0.313);
        z := foo(-1.0);

        res := x * z + y * z + y*z;

        println("res = %", res);
    }
    {
        restart();

        // * since `Foo` takes any `T: Type`, we can just give it the symbolic `ID`
        // * instead, and all the math will be constructed symbolically
        x : Foo(ID) = .{val = var()};
        // x := foo(var());
        y := foo(var());
        z := foo(var());

        res := x * z + y * z + y*z;

        println("res = %\n%", res, str(res.val));
    }


    // * Do the simplification during compile-time and generate the optimized
    // * code; then call the optimized function `foo`.
    #insert -> string {
        restart();

        x := var();
        y := var();
        z := var();
        println("x = %", str(x));
        println("y = %", str(y));
        println("z = %", str(z));

        println("x+y = %", str(x+y));
        println("x*y = %", str(x*y));
        println("x*y = %, reoccuring operator recognized?", str(x*y));

        res := x*z + y*z + y*z;
        
        println("res = %", res);
        println("%", str(res));

        res2 := simplify(res);
        println("res2 = %", res2);
        println("%", str(res2));

        func_src := compile(res2);
        println("compiled code:\n%", func_src);

        func := make_function("foo", res2);
        println("generated function code:\n%", func);

        return func;
    };

    println("foo(42.0, 0.313, -1.0) = %", foo(42.0, 0.313, -1.0));
}
#import "Basic";
#import "Math";
#import "Sort";
#load "utils.jai";
#module_parameters(_POS: int, _NEG: int, _ZER: int, T := float32, _DEFAULT_DUAL:="poincare");
POS :: _POS;
NEG :: _NEG;
ZER :: _ZER;
#if _DEFAULT_DUAL == "poincare" then DEFAULT_DUAL :: DualType.POINCARE;
#if _DEFAULT_DUAL == "hodge"    then DEFAULT_DUAL :: DualType.HODGE;
#if _DEFAULT_DUAL == "polar"    then DEFAULT_DUAL :: DualType.POLAR;

#if POS == 2 && NEG == 0 && ZER == 1 #load "pga.jai";

one  :: #run cast(T)1;
zero :: #run cast(T)0;

N :: POS + NEG + ZER;
N_ELEMENTS :: #run powi(2, N+1);

#if N <= 8 {
    B :: u8;
} else #if N <= 16 {
    B :: u16;
} else #if N <= 32 {
    B :: u32; 
} else {
    B :: u64;
    #run assert(N <= 64, "Only algebras of size up to 64 are supported, but chosen %.", N);
}

#run print("Generating algebra with basis (%,%,%). Basis elements represented as %\n", POS, NEG, ZER, B);


fmt :: inline (x: T) -> string {
    return sprint("%", x);
}

fmt :: (n: B, $dual:= false) -> string {
    builder: String_Builder;

    if n == 0
        return "";

    indx := 1;
    x := n;
    while x > 0 {
        if x & 1 #if dual {
            print_to_builder(*builder, "E%", indx);
        } else {
            print_to_builder(*builder, "e%", indx);
        }
        x = x >> 1;
        indx += 1;
    }

    return builder_to_string(*builder);
}

// *  ##     ## ##     ## ##       ######## #### ##     ## ########  ######  ########  #######  ######## 
// *  ###   ### ##     ## ##          ##     ##  ##     ## ##       ##    ##    ##    ##     ## ##     ## 
// *  #### #### ##     ## ##          ##     ##  ##     ## ##       ##          ##    ##     ## ##     ## 
// *  ## ### ## ##     ## ##          ##     ##  ##     ## ######   ##          ##    ##     ## ######## 
// *  ##     ## ##     ## ##          ##     ##   ##   ##  ##       ##          ##    ##     ## ##   ## 
// *  ##     ## ##     ## ##          ##     ##    ## ##   ##       ##    ##    ##    ##     ## ##    ## 
// *  ##     ##  #######  ########    ##    ####    ###    ########  ######     ##     #######  ##     ## 

assert_duality :: inline (a: $_A/MultiVector, b: $_B/MultiVector) #expand {
    #run assert(
        a.poincare == b.poincare, 
        "Cannot add multivectors of different duality: \na Poincare dual: %\nb Poincare dual: %", 
        a.poincare, b.poincare);
}

MultiVector :: struct(bases: ..B, poincare: bool = false) {
    values: [bases.count]T;
}

e0 :: MultiVector(0).{.[one]};

MV_I     :: MultiVector(1<<N - 1);
// MV_I_inv :: MultiVector(1<<N - 1, poincare=true);


eI :: MV_I.{.[one]};
eI_polar :: #run ifx (eI*eI).values[0] < 0 then -1*eI else eI;
EI :: eI_polar;

#run {
    println("eI       : % = %", type_of(eI), fmt(eI));
    println("eI_polar : % = %", type_of(eI_polar), fmt(eI_polar));
}

#insert -> string {
    builder : String_Builder;

    append(*builder, "FULL_MV :: MultiVector(0");
    for i: 1..1<<N-1 {
        print_to_builder(*builder, ", %", i);
    }
    append(*builder, ");\n");

    // * full mv zero
    append(*builder, "MV_ZERO :: FULL_MV.{.[zero");
    for i: 1..1<<N-1 {
        append(*builder, ", zero");
    }
    append(*builder, "]};\n");

    return builder_to_string(*builder);
};

#insert -> string {
    DEBUG :: false;

    builder : String_Builder;

    for grade: 1..N {
        #if DEBUG print("grade %:\n", grade);

        NN := powi(N, grade);
        #if DEBUG print("    NN = %\n", NN);

        indices := NewArray(grade, int); defer array_free(indices);

        // * iterate over all possible index combinations
        for i: 0..NN-1 {
            #if DEBUG print("    %: ", i);

            // * increase indices
            indices[0] = i;
            for g: 0..grade-1-1 {
                if indices[g] != 0 && indices[g] % N == 0 then
                    indices[g+1] += 1;
            }
            for *indices { <<it %= N; }
            #if DEBUG print("%\n", indices);

            // * skip if the same element appears twice
            for i0: 0..grade-1 {
                for j0: i0+1..grade-1 {
                    if indices[i0] == indices[j0] {
                        #if DEBUG print("        same elements % & %, continue\n", i0, j0);
                        continue i;
                    }
                }
            }

            if grade == 1 {
                for indices
                    print_to_builder(*builder, "e% :: MultiVector(%).{.[one]};\n", it+1, 1 << it);
            } else {
                append(*builder, "e");
                for indices {
                    print_to_builder(*builder, "%", it+1);
                }
                append(*builder, " :: #run ");
                for indices {
                    print_to_builder(*builder, "e%", it+1);
                    if it_index != indices.count-1 then
                        append(*builder, " * ");
                }
                append(*builder, ";\n");
            }
        }
    }

    return builder_to_string(*builder);
}

fmt :: (mv: $V/MultiVector) -> string {
    builder : String_Builder;

    for b, i: V.bases {
        if i > 0 && mv.values[i] >= 0
            append(*builder, " +");
        print_to_builder(*builder, " %'%", mv.values[i], fmt(b, V.poincare));
    }

    return builder_to_string(*builder);
}


// *  ######  ######## ##       ########  ######  ########
// * ##    ## ##       ##       ##       ##    ##    ##
// * ##       ##       ##       ##       ##          ##
// *  ######  ######   ##       ######   ##          ##
// *       ## ##       ##       ##       ##          ##
// * ##    ## ##       ##       ##       ##    ##    ##
// *  ######  ######## ######## ########  ######     ##


SELECT :: struct(a: Type, grade: int) {
    MV :: #run select(a,grade);
}
select :: ($a: Type, $grade: int) -> Type {
    N_select :: #run () -> int {
        res := 0;
        for a.bases {
            if hamming_weight(it) == grade then res += 1;
        }
        return res;
    }();

    #if N_select == 0 {
        return MultiVector(0);
    } else {
        b :: #run () -> [N_select]B {
            res : [N_select]B = ---;
            indx := 0;
            for a.bases {
                if hamming_weight(it) == grade {
                    res[indx] = it;
                    indx += 1;
                }
            } 
            return res;
        }();

        return MultiVector(..b);
    }
}

select :: inline (a: $_A/MultiVector, $grade: int) -> SELECT(_A, grade).MV {
    res : SELECT(_A, grade).MV = ---;

    #insert #run () -> string {
        builder : String_Builder;
        indx := 0;
        for a.bases {
            if hamming_weight(it) == grade {
                print_to_builder(*builder, "res.values[%] = a.values[%];\n", indx, it_index);
                indx += 1;
            }
        }
        if indx == 0 {  // * no matching grade found
            append(*builder, "res.values[0] = 0;\n");
        }
        return builder_to_string(*builder);
    }();

    return res;
}


select :: (a: $_A/MultiVector, $MV: Type) -> MV {
    #if _A == MV {
        return a;
    } else {
        res : MV;

        #insert #run () -> string {
            builder : String_Builder;

            for res.bases {
                for ab, ai: a.bases {
                    if it == ab {
                        print_to_builder(*builder, "res.values[%] = a.values[%];\n", it_index, ai);
                        break;
                    }
                }
            }

            return builder_to_string(*builder);
        }();

        return res;
    }
}

select :: inline (a: $_A/MultiVector, b: $_B/MultiVector) -> _B {
    return select(a, _B);
}

// *              ##
// *   ##        ##
// *   ##       ##
// * ######    ##    #######
// *   ##     ##
// *   ##    ##
// *        ##

ADD :: struct(a: Type, b: Type) {
    MV :: #run add(a, b);
}
add :: ($atype: Type, $btype: Type) -> Type {
    DEBUG :: false;

    N_combined :: #run () -> int {
        res := atype.bases.count;

        for bb : btype.bases {
            for ab : atype.bases {
                if bb == ab then continue bb;
            }
            res += 1;
        }

        return res;
    }();

    #if DEBUG println("N_combined = % from % and %", N_combined, atype.bases, btype.bases);

    b :: #run () -> [N_combined]B {
        res : [N_combined]B = ---;

        a_indx := 0;
        b_indx := 0;

        for i: 0..N_combined-1 {
            if a_indx < atype.bases.count && b_indx < btype.bases.count {
                ab := atype.bases[a_indx];
                bb := btype.bases[b_indx];
                if ab < bb {
                    res[i] = ab;
                    a_indx += 1;
                } else if ab > bb {
                    res[i] = bb;
                    b_indx += 1;
                } else { // ab == bb
                    res[i] = ab;
                    a_indx += 1;
                    b_indx += 1;
                }
            } else if a_indx < atype.bases.count {
                ab := atype.bases[a_indx];
                res[i] = ab;
                a_indx += 1;
            } else if b_indx < btype.bases.count {
                bb := btype.bases[b_indx];
                res[i] = bb;
                b_indx += 1;
            }
        }
        return res;
    }();
    return MultiVector(..b);
}


operator + :: (a: $_A/MultiVector, b: $_B/MultiVector) -> ADD(_A, _B).MV {
    assert_duality(a,b);
    res : ADD(_A, _B).MV;
    // * assuming the bases are sorted from small to large value
    index_a := 0;
    index_b := 0;
    for res.bases {
        while a.bases[index_a] < it && index_a < a.bases.count-1   index_a += 1;
        while b.bases[index_b] < it && index_b < b.bases.count-1   index_b += 1;
        if a.bases[index_a] == it res.values[it_index] += a.values[index_a];
        if b.bases[index_b] == it res.values[it_index] += b.values[index_b];
    }
    return res;
}
operator - :: (a: $_A/MultiVector, b: $_B/MultiVector) -> ADD(_A, _B).MV {
    assert_duality(a,b);
    res : ADD(_A, _B).MV;
    // * assuming the bases are sorted from small to large value
    index_a := 0;
    index_b := 0;
    for res.bases {
        while a.bases[index_a] < it && index_a < a.bases.count-1   index_a += 1;
        while b.bases[index_b] < it && index_b < b.bases.count-1   index_b += 1;
        if a.bases[index_a] == it res.values[it_index] += a.values[index_a];
        if b.bases[index_b] == it res.values[it_index] -= b.values[index_b];
    }
    return res;
}



operator + :: (a: $_A/MultiVector, b: T) -> ADD(MultiVector(0), _A).MV #symmetric {
    tmp : MultiVector(0, poincare=a.poincare) = ---;
    tmp.values[0] = b;
    return a + tmp;
}
operator - :: (a: $_A/MultiVector, b: T) -> ADD(MultiVector(0), _A).MV {
    tmp : MultiVector(0, poincare=a.poincare) = ---;
    tmp.values[0] = b;
    return a - tmp;
}
operator - :: (b: T, a: $_A/MultiVector) -> ADD(MultiVector(0), _A).MV {
    tmp : MultiVector(0, poincare=a.poincare) = ---;
    tmp.values[0] = b;
    return tmp - a;
}


operator * :: (a: $_A/MultiVector, b: T) -> _A #symmetric {
    res : _A = ---;
    for a.values {
        res.values[it_index] = it * b;
    }
    return res;
}
operator / :: (a: $_A/MultiVector, b: T) -> _A {
    res : _A = ---;
    for a.values {
        res.values[it_index] = it / b;
    }
    return res;
}


// * ########  ######## ##     ## ######## ########   ######  ########
// * ##     ## ##       ##     ## ##       ##     ## ##    ## ##
// * ##     ## ##       ##     ## ##       ##     ## ##       ##
// * ########  ######   ##     ## ######   ########   ######  ######
// * ##   ##   ##        ##   ##  ##       ##   ##         ## ##
// * ##    ##  ##         ## ##   ##       ##    ##  ##    ## ##
// * ##     ## ########    ###    ######## ##     ##  ######  ########


reverse :: (a: $_A/MultiVector) -> _A {
    res := a;
    #insert #run () -> string {
        builder : String_Builder;
        for ab, ai: a.bases {
            hw := hamming_weight(ab);
            if hw < 2 continue;
            binom := binomial(hw, 2);
            if binom % 2 == 1 {
                print_to_builder(*builder, "res.values[%] = -res.values[%];\n", ai, ai);
            }
        }
        return builder_to_string(*builder);
    }();
    return res;
}



// *  #######  ##     ## ######## ######## ########
// * ##     ## ##     ##    ##    ##       ##     ##
// * ##     ## ##     ##    ##    ##       ##     ##
// * ##     ## ##     ##    ##    ######   ########
// * ##     ## ##     ##    ##    ##       ##   ##
// * ##     ## ##     ##    ##    ##       ##    ##
// *  #######   #######     ##    ######## ##     ##

OUTER :: struct(a: Type, b: Type) {
    MV :: #run outer(a,b);
}
outer :: ($a: Type, $b: Type) -> Type {
    o :: GEO(a, b).MV;
    grade :: #run hamming_weight(o.bases[o.bases.count-1]);
    return select(o, grade);
}


meet :: outer;
outer :: inline (a: $_A/MultiVector, b: $_B/MultiVector) -> OUTER(_A, _B).MV {
    assert_duality(a,b);
    DEBUG :: false;
    grade :: #run hamming_weight(OUTER(_A, _B).MV.bases[0]);

    #if DEBUG {
              println("meet --------------------------");
        defer println("meet ^^^^^^^^^^^^^^^^^^^^^^^^^^");

        println("a = %", fmt(a));
        println("b = %", fmt(b));
        println("grade = %", grade);
        println("a * b = %", fmt(a*b));
    }

    return select(a * b, grade);
}



// *  ######   ########  #######     ########  ########   #######  ########
// * ##    ##  ##       ##     ##    ##     ## ##     ## ##     ## ##     ##
// * ##        ##       ##     ##    ##     ## ##     ## ##     ## ##     ##
// * ##   #### ######   ##     ##    ########  ########  ##     ## ##     ##
// * ##    ##  ##       ##     ##    ##        ##   ##   ##     ## ##     ##
// * ##    ##  ##       ##     ##    ##        ##    ##  ##     ## ##     ##
// *  ######   ########  #######     ##        ##     ##  #######  ########

GEO :: struct(a: Type, b: Type) {
    MV :: #run geo(a, b);
}

geo :: ($atype: Type, $btype: Type) -> Type {
    DEBUG :: false;

    N_combined :: #run () -> int {
        res : [..]B;

        for bb : btype.bases {
            for ab : atype.bases {
                combined := ab ^ bb;

                overlap := ab & bb;
                if overlap >> (POS+NEG) continue; // * basis squares to zero

                array_add_if_unique(*res, combined);
            }
        }

        #if DEBUG println("combined bases %", res);

        return res.count;
    }();

    #if DEBUG println("N_combined = % from % and %", N_combined, atype.bases, btype.bases);

    if N_combined == 0 return MultiVector(0);

    b :: #run () -> [N_combined]B {
        res : [N_combined]B = ---;

        tmp : [..]B;
        for bb : btype.bases {
            for ab : atype.bases {
                combined := ab ^ bb;

                overlap := ab & bb;
                if overlap >> (POS+NEG) continue; // * basis squares to zero

                array_add_if_unique(*tmp, combined);
            }
        }
        cmp :: (a: B, b: B) -> s64 {
            return cast(s64)(a > b);
        }
        bubble_sort(tmp, cmp);

        for tmp {
            res[it_index] = it;
        }

        return res;
    }();
    #if DEBUG println("combined bases %", b);
    return MultiVector(..b);
}

operator * :: inline (a: $_A/MultiVector, b: $_B/MultiVector) -> GEO(_A,_B).MV {
    return geo(a, b);
}



geo :: (a: $_A/MultiVector, b: $_B/MultiVector) -> GEO(_A,_B).MV {
    assert_duality(a,b);
    res : GEO(_A,_B).MV;
    // * assuming the bases are sorted from small to large value

    #insert #run ()->string {
        DEBUG :: false;
        builder : String_Builder;

        #if DEBUG println("mv % and % , res % ------------------------", a.bases, b.bases, res.bases);
        for ab, ai: a.bases {
            for bb, bi: b.bases {
                #if DEBUG println("checking bases % and %", ab, bb);
                if ab == bb {
                    #if DEBUG println("equal bases %", ab);
                    if ab >> (NEG+POS) & ((1 << ZER) - 1) > 0 then continue;

                    hw_negs := hamming_weight(ab >> POS & ((1 << NEG)-1));

                    if hw_negs % 2 == 0 {
                        print_to_builder(*builder, "res.values[0] += a.values[%] * b.values[%];\n", ai, bi);
                    } else {
                        print_to_builder(*builder, "res.values[0] -= a.values[%] * b.values[%];\n", ai, bi);
                    }
                } else {
                    #if DEBUG println("non-equal bases % and %", ab, bb);
                    overlap  := ab & bb;
                    combined := ab ^ bb;
                    #if DEBUG println("overlap %, combined %", overlap, combined);

                    // ? An algorithm to find the number of neighbor swaps necessary to sort an array:
                    // ? https://www.geeksforgeeks.org/number-swaps-sort-adjacent-swapping-allowed/
                    
                    if overlap == 0 { // * no basis elements square up
                        for res.bases {
                            if it == combined {

                                n_swaps := count_swaps(ab, bb);

                                if n_swaps % 2 == 0 {
                                    print_to_builder(*builder, "res.values[%] += a.values[%] * b.values[%];\n", it_index, ai, bi);
                                } else {
                                    print_to_builder(*builder, "res.values[%] -= a.values[%] * b.values[%];\n", it_index, ai, bi);
                                }
                                
                                break;
                            }
                        }
                        continue;
                    }

                    // * overlap != 0

                    if overlap >> (NEG+POS) > 0 then continue;
                    #if DEBUG println("no zeros");

                    // * the negative signs from squaring bases that square to -1
                    hw_negs := hamming_weight(overlap >> POS & ((1 << NEG)-1));
                    #if DEBUG println("hamming weight of negs: %", hw_negs);

                    for res.bases {
                        #if DEBUG println("checking res basis %", it);
                        if it == combined {

                            n_swaps := count_swaps(ab, bb);

                            #if DEBUG println("found basis % in res at index %", it, it_index);
                            if (hw_negs+n_swaps) % 2 == 0
                                print_to_builder(*builder, "res.values[%] += a.values[%] * b.values[%];\n", it_index, ai, bi);
                            else
                                print_to_builder(*builder, "res.values[%] -= a.values[%] * b.values[%];\n", it_index, ai, bi);
                            
                            break;
                        }
                    }
                }
            }
        }

        return builder_to_string(*builder);
    }();

    return res;
}

// * #### ##    ## ##    ## ######## ########
// *  ##  ###   ## ###   ## ##       ##     ##
// *  ##  ####  ## ####  ## ##       ##     ##
// *  ##  ## ## ## ## ## ## ######   ########
// *  ##  ##  #### ##  #### ##       ##   ##
// *  ##  ##   ### ##   ### ##       ##    ##
// * #### ##    ## ##    ## ######## ##     ##

INNER :: struct(a: Type, b: Type) {
    MV :: #run inner(a,b);
}
inner :: ($a: Type, $b: Type) -> Type {
    grade_a :: #run hamming_weight(a.bases[a.bases.count-1]);
    grade_b :: #run hamming_weight(b.bases[b.bases.count-1]);
    return select(#run geo(a, b), #run abs(grade_a - grade_b));
}

dot :: inner;
inner :: (a: $_A/MultiVector, b: $_B/MultiVector) -> INNER(_A, _B).MV {
    assert_duality(a,b);
    DEBUG :: false;
    grade_a :: #run hamming_weight(a.bases[a.bases.count-1]);
    grade_b :: #run hamming_weight(b.bases[b.bases.count-1]);
    #if DEBUG       println("inner() ----------------------------");
    #if DEBUG defer println("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^");
    #if DEBUG println("a = %", fmt(a));
    #if DEBUG println("b = %", fmt(b));
    #if DEBUG println("grade a = %, grade b = %", grade_a, grade_b);
    #if DEBUG println("return grade %", abs(grade_a - grade_b));
    #if DEBUG println("a * b = %", fmt(a * b));
    return select(a * b, #run abs(grade_a - grade_b));
}


// *  ######   #######  ##     ## ##     ## ##     ## ########    ###    ########  #######  ########
// * ##    ## ##     ## ###   ### ###   ### ##     ##    ##      ## ##      ##    ##     ## ##     ##
// * ##       ##     ## #### #### #### #### ##     ##    ##     ##   ##     ##    ##     ## ##     ##
// * ##       ##     ## ## ### ## ## ### ## ##     ##    ##    ##     ##    ##    ##     ## ########
// * ##       ##     ## ##     ## ##     ## ##     ##    ##    #########    ##    ##     ## ##   ##
// * ##    ## ##     ## ##     ## ##     ## ##     ##    ##    ##     ##    ##    ##     ## ##    ##
// *  ######   #######  ##     ## ##     ##  #######     ##    ##     ##    ##     #######  ##     ##

COMM :: struct(a: Type, b: Type) {
    MV :: #run comm(a,b);
}

comm :: ($a: Type, $b: Type) ->  Type {
    return geo(a, b);
}

cross :: comm;
comm :: (a: $_A/MultiVector, b: $_B/MultiVector) -> COMM(_A,_B).MV {
    assert_duality(a,b);
    return cast(T)0.5 * (a*b - b*a);
}



// * ########  ##     ##    ###    ##
// * ##     ## ##     ##   ## ##   ##
// * ##     ## ##     ##  ##   ##  ##
// * ##     ## ##     ## ##     ## ##
// * ##     ## ##     ## ######### ##
// * ##     ## ##     ## ##     ## ##
// * ########   #######  ##     ## ########




DualType :: enum {
    POLAR;
    POINCARE;
    HODGE;
}

DUAL :: struct(a: Type, type: DualType) {
    MV :: #run dual(a, type);
}

dual :: ($a: Type, $type:=DEFAULT_DUAL) -> Type {
    #if type == .POLAR    then return dual_polar(a);
    #if type == .POINCARE then return dual_neg_type(a, true);
    #if type == .HODGE    then return dual_neg_type(a, false);
}

dual :: inline (a: $_A/MultiVector, $type:=DEFAULT_DUAL) -> DUAL(_A, type).MV {
    #if type == .POLAR    then return dual_polar(a);
    #if type == .POINCARE then return dual_poincare(a);
    #if type == .HODGE    then return dual_hodge(a);
}

// TODO inverse dual


INV_DUAL :: struct(a: Type, type: DualType) {
    MV :: #run inv_dual(a, type);
}

inv_dual :: ($a: Type, $type:=DEFAULT_DUAL) -> Type {
    #if type == .POLAR    then return inv_dual_polar(a);
    #if type == .POINCARE then return dual_neg_type(a, true);
    #if type == .HODGE    then return dual_neg_type(a, false);
}

inv_dual :: inline (a: $_A/MultiVector, $type:=DEFAULT_DUAL) -> INV_DUAL(_A, type).MV {
    #if type == .POLAR    then return inv_dual_polar(a);
    #if type == .POINCARE then return inv_dual_poincare(a);
    #if type == .HODGE    then return inv_dual_hodge(a);
}


// * ########   #######  ##          ###    ########  #### ######## ##    ##
// * ##     ## ##     ## ##         ## ##   ##     ##  ##     ##     ##  ##
// * ##     ## ##     ## ##        ##   ##  ##     ##  ##     ##      ####
// * ########  ##     ## ##       ##     ## ########   ##     ##       ##
// * ##        ##     ## ##       ######### ##   ##    ##     ##       ##
// * ##        ##     ## ##       ##     ## ##    ##   ##     ##       ##
// * ##         #######  ######## ##     ## ##     ## ####    ##       ##


POLAR :: struct(a: Type) {
    MV :: #run dual_polar(a);
}
dual_polar :: ($a: Type) -> Type {
    return geo(type_of(eI), a);
}
dual_polar :: inline (a: $_A/MultiVector) -> POLAR(_A).MV {
    return eI * a;
}

INV_POLAR :: struct(a: Type) {
    MV :: #run inv_dual_polar(a);
}
inv_dual_polar :: ($a: Type) -> Type {
    return geo(type_of(EI), a);
}

inv_dual_polar :: inline (a: $_A/MultiVector) -> INV_POLAR(_A).MV {
    return EI * a;
}



// * ########  ##     ##    ###    ##
// * ##     ## ##     ##   ## ##   ##
// * ##     ## ##     ##  ##   ##  ##
// * ##     ## ##     ## ##     ## ##
// * ##     ## ##     ## ######### ##
// * ##     ## ##     ## ##     ## ##
// * ########   #######  ##     ## ########

dual_neg_type :: ($a: Type, $poincare:=false) -> Type {
    b :: #run () -> [a.bases.count]B {
        res : [a.bases.count]B = ---;
        max :: (1 << N) - 1;
        for ab, ai: a.bases {
            res[ai] = ~ab & max;
        }

        cmp :: (a: B, b: B) -> s64 {
            return cast(s64)(a > b);
        }
        bubble_sort(res, cmp);

        return res;
    }();

    #if poincare
        return MultiVector(..b, poincare=!a.poincare);
    else
        return MultiVector(..b, poincare=a.poincare);
}

dual_construct_assignments :: ($a: Type, $res: Type) -> string {
    builder : String_Builder;

    max :: (1 << N)-1;
    for ab, ai: a.bases {
        dual := ~ab & max;
        for rb, ri: res.bases {
            if dual != rb then continue;

            print_to_builder(*builder, "res.values[%] = a.values[%];\n", ri, ai);
        }
    }

    return builder_to_string(*builder); // * flip duality bit
}



// * ########   #######  #### ##    ##  ######     ###    ########  ########
// * ##     ## ##     ##  ##  ###   ## ##    ##   ## ##   ##     ## ##      
// * ##     ## ##     ##  ##  ####  ## ##        ##   ##  ##     ## ##      
// * ########  ##     ##  ##  ## ## ## ##       ##     ## ########  ######  
// * ##        ##     ##  ##  ##  #### ##       ######### ##   ##   ##      
// * ##        ##     ##  ##  ##   ### ##    ## ##     ## ##    ##  ##      
// * ##         #######  #### ##    ##  ######  ##     ## ##     ## ########


POINCARE :: struct(a: Type) {
    // MV :: #run poincare(a);
    MV :: #run dual_neg_type(a, poincare=true);
}

// poincare :: ($a: Type) -> Type {
//     b :: #run () -> [a.bases.count]B {
//         res : [a.bases.count]B = ---;
//         max :: (1 << N) - 1;
//         for ab, ai: a.bases {
//             res[ai] = ~ab & max;
//         }

//         cmp :: (a: B, b: B) -> s64 {
//             return cast(s64)(a > b);
//         }
//         bubble_sort(res, cmp);

//         return res;
//     }();

//     return MultiVector(..b, poincare=!a.poincare);
// }

dual_poincare :: (a: $_A/MultiVector) -> POINCARE(_A).MV {
    res : POINCARE(_A).MV = ---;

    #insert #run dual_construct_assignments(_A, POINCARE(_A).MV);

    // #insert #run () -> string {
    //     builder : String_Builder;

    //     max :: (1 << N)-1;
    //     for ab, ai: a.bases {
    //         dual := ~ab & max;
    //         for rb, ri: res.bases {
    //             if dual != rb then continue;

    //             print_to_builder(*builder, "res.values[%] = a.values[%];\n", ri, ai);
    //         }
    //     }

    //     return builder_to_string(*builder); // * flip duality bit
    // }();

    return res;
}

INV_POINCARE :: struct(a: Type) {
    MV :: #run dual_neg_type(a, poincare=true);
}

inv_dual_poincare :: (a: $_A/MultiVector) -> INV_POINCARE(_A).MV {
    res : INV_POINCARE(_A).MV = ---;
    #insert #run dual_construct_assignments(_A, INV_POINCARE(_A).MV);
    return res;
}

// * ########  ##     ##    ###    ##
// * ##     ## ##     ##   ## ##   ##
// * ##     ## ##     ##  ##   ##  ##
// * ##     ## ##     ## ##     ## ##
// * ##     ## ##     ## ######### ##
// * ##     ## ##     ## ##     ## ##
// * ########   #######  ##     ## ########

// * ##     ##  #######  ########   ######   ########
// * ##     ## ##     ## ##     ## ##    ##  ##
// * ##     ## ##     ## ##     ## ##        ##
// * ######### ##     ## ##     ## ##   #### ######
// * ##     ## ##     ## ##     ## ##    ##  ##
// * ##     ## ##     ## ##     ## ##    ##  ##
// * ##     ##  #######  ########   ######   ########

HODGE :: struct(a: Type) {
    // MV :: #run hodge(a);
    MV :: #run dual_neg_type(a, poincare=false);
}

// hodge :: ($a: Type) -> Type {
//     b :: #run () -> [a.bases.count]B {
//         res : [a.bases.count]B = ---;
//         max :: (1 << N) - 1;
//         for ab, ai: a.bases {
//             res[ai] = ~ab & max;
//         }

//         cmp :: (a: B, b: B) -> s64 {
//             return cast(s64)(a > b);
//         }
//         bubble_sort(res, cmp);

//         return res;
//     }();

//     return MultiVector(..b, poincare=a.poincare); // * same duality bit
// }

dual_hodge :: (a: $_A/MultiVector) -> HODGE(_A).MV {
    res : HODGE(_A).MV = ---;

    #insert #run dual_construct_assignments(_A, HODGE(_A).MV);

    // #insert #run () -> string {
    //     builder : String_Builder;

    //     max :: (1 << N)-1;
    //     for ab, ai: a.bases {
    //         dual := ~ab & max;
    //         for rb, ri: res.bases {
    //             if dual != rb then continue;

    //             print_to_builder(*builder, "res.values[%] = a.values[%];\n", ri, ai);
    //         }
    //     }

    //     return builder_to_string(*builder);
    // }();

    return res;
}

INV_HODGE :: struct(a: Type) {
    MV :: #run dual_neg_type(a, poincare=false);
}

inv_dual_hodge :: (a: $_A/MultiVector) -> INV_HODGE(_A).MV {
    res : INV_HODGE(_A).MV = ---;
    #insert #run dual_construct_assignments(_A, INV_HODGE(_A).MV);
    return res;
}

// * ########  ########  ######   ########  ########  ######   ######  #### ##     ## ########
// * ##     ## ##       ##    ##  ##     ## ##       ##    ## ##    ##  ##  ##     ## ##
// * ##     ## ##       ##        ##     ## ##       ##       ##        ##  ##     ## ##
// * ########  ######   ##   #### ########  ######    ######   ######   ##  ##     ## ######
// * ##   ##   ##       ##    ##  ##   ##   ##             ##       ##  ##   ##   ##  ##
// * ##    ##  ##       ##    ##  ##    ##  ##       ##    ## ##    ##  ##    ## ##   ##
// * ##     ## ########  ######   ##     ## ########  ######   ######  ####    ###    ########

REGRESSIVE :: struct(a: Type, b: Type, $dual_type:=DEFAULT_DUAL) {
    MV :: #run regressive(a, b, dual_type);
}
regressive :: ($a: Type, $b: Type, $dual_type: DualType) -> Type {
    // return #run hodge(#run outer(#run hodge(a), #run hodge(b)));
    return #run inv_dual(#run outer(#run dual(a, dual_type), #run dual(b, dual_type)), dual_type);
}

join :: regressive;
regressive :: inline (a: $_A/MultiVector, b: $_B/MultiVector, $dual_type:=DEFAULT_DUAL) -> REGRESSIVE(_A, _B, dual_type).MV {
    assert_duality(a,b);
    DEBUG :: false;
    #if DEBUG {
        println(      "join/regressive ---------------------");
        defer println("join/regressive ^^^^^^^^^^^^^^^^^^^^^");
        println("a = %", fmt(a));
        println("b = %", fmt(b));
        println("dual a = %", fmt(dual_hodge(a)));
        println("dual b = %", fmt(dual_hodge(b)));
        r := meet(dual_hodge(a), dual_hodge(b));
        println("Point == type_of(r): %", Point == type_of(r));
        println("Point : %", Point);
        println("meet : % = %", type_of(r), r);
        // println("     = %", fmt_euclid(r));
    }
    // return dual_hodge(meet(dual_hodge(a), dual_hodge(b)));
    return inv_dual(meet(dual(a, dual_type), dual(b, dual_type)), dual_type);
}
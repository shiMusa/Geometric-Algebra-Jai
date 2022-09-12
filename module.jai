#import "Basic";
#import "Math";
#load "utils.jai";
#module_parameters(POS: int, NEG: int, ZER: int, T := float32, STACK_ALLOC := false);

one  :: #run cast(T)1;
zero :: #run cast(T)0;

N :: POS + NEG + ZER;
#if N <= 8 {
    B :: u8;
} else #if N <= 16 {
    B :: u16;
} else #if N <= 32 {
    B :: u32; 
} else {
    B :: u64;
    #run assert(N < 64, "Only algebras of size up to 64 are supported, but chosen %.", N);
}

#run print("Generating algebra with basis (%,%,%). Basis elements represented as %\n", POS, NEG, ZER, B);

BasisNumber :: struct(basis: B) {
    value: T;
}
BN0 :: BasisNumber(0);
bn :: (value: T, $basis: B = 0) -> BasisNumber(basis) {
    res : BasisNumber(B) = ---;
    res.value = value;
}
fmt :: (n: BasisNumber, $nosign := false) -> string {
    builder: String_Builder;

    #if nosign {
        v := ifx n.value >= 0 then n.value else -n.value;
        print_to_builder(*builder, "% ", v);
    } else 
        print_to_builder(*builder, "% ", n.value);

    for 0..N-1 {
        // append(*builder, "e");
        if n.basis & (1 << it) {
            print_to_builder(*builder, "e%", it+1);
        }
    }

    return builder_to_string(*builder);
}
fmt :: inline (n: T) -> string {
    return sprint("%", n);
}

e0 :: BasisNumber(0).{one};



basis_to_T :: inline (a: B, b: B) -> bool {
    return (a & b) >> (POS+NEG) > 0 || a ^ b == 0;
}

// if the multiplication results in a real number (T)
operator * :: (a: $_A/BasisNumber, b: $_B/BasisNumber) -> T #modify {
    ba := get_const_value(_A, "basis", B);
    bb := get_const_value(_B, "basis", B);
    return basis_to_T(ba, bb);
} {
    DEBUG :: false;
    #if DEBUG print("% * % -> %\n", fmt(a), fmt(b), T);
    
    overlap :: _A.basis & _B.basis;
    #if DEBUG print("overlap %\n", formatInt(overlap,2));
    
    #if overlap >> (POS+NEG) > 0 {
        // a basis element squares to zero
        return zero;
    } else {
        // all basis elements square to +-1
        #if DEBUG println("xor: %", _A.basis ^ _B.basis);
        parity :: #run one - cast(T)2 * cast(T)(_A.basis > _B.basis);
        negs :: #run (overlap >> POS) & ((1 << NEG)-1);
        negs_squared :: #run 1 - 2 * ((hamming_weight(negs)) % 2);
        return a.value * b.value * parity * negs_squared;
    }
}

// if the multiplication results in a BasisNumber
operator * :: (a: $_A/BasisNumber, b: $_B/BasisNumber) -> BasisNumber(_A.basis ^ _B.basis)
#modify {
    ba := get_const_value(_A, "basis", B);
    bb := get_const_value(_B, "basis", B);
    return !basis_to_T(ba, bb);
} {
    // #run assert(type_info(C).type == .STRUCT, "ERROR : no compatible multiplication found! C = %", C);
    res : BasisNumber(_A.basis ^ _B.basis) = ---;
    
    DEBUG :: false;
    #if DEBUG print("% * % -> %, %\n", fmt(a), fmt(b), res, formatInt(res.basis,2));
    
    overlap :: _A.basis & _B.basis;
    #if DEBUG print("overlap %\n", formatInt(overlap,2));
    #if DEBUG println("xor: %", _A.basis ^ _B.basis);


    // default ordering from small basis number to large
    parity :: #run one - cast(T)2 * cast(T)(_A.basis > _B.basis);
    res.value = a.value * b.value * parity;
    #if DEBUG println("res : % = %", type_of(res), res);

    // no basis elements square up
    #if overlap == 0 {
        return res;
    } else {
        // basis elements that square to -1
        negs :: #run (overlap >> POS) & ((1 << NEG)-1);
        #if DEBUG print("negs %\n", formatInt(negs,2));
        #if negs {
            negs_squared :: #run 1 - 2 * ((hamming_weight(negs)) % 2);
            res.value *= negs_squared;
        }

        // the rest squares to +1
        return res;
    }
}

#insert -> string {
    DEBUG :: false;

    builder : String_Builder;

    for grade: 1..N {
        #if DEBUG print("grade %:\n", grade);

        NN := powi(N, grade);
        #if DEBUG print("    NN = %\n", NN);

        indices := NewArray(grade, int); defer array_free(indices);

        // iterate over all possible index combinations
        for i: 0..NN-1 {
            #if DEBUG print("    %: ", i);

            // increase indices
            indices[0] = i;
            for g: 0..grade-1-1 {
                if indices[g] != 0 && indices[g] % N == 0 then
                    indices[g+1] += 1;
            }
            for *indices { <<it %= N; }
            #if DEBUG print("%\n", indices);

            // skip if the same element appears twice
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
                    print_to_builder(*builder, "e% :: BasisNumber(%).{one};\n", it+1, 1 << it);
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

operator - :: (a: $_A/BasisNumber) -> _A {
    res : _A = ---;
    res.value = -a.value;
    return res;
}
operator * :: (a: $_A/BasisNumber, b: T) -> _A #symmetric {
    res : _A = ---;
    res.value = a.value * b;
    return res;
}

operator / :: (a: $_A/BasisNumber, b: T) -> _A {
    res : _A = ---;
    res.value = a.value / b;
    return res;
}


N_ELEMENTS :: #run powi(2, N+1);


// ##     ## ##     ## ##       ######## #### ##     ## ########  ######  ########  #######  ########
// ###   ### ##     ## ##          ##     ##  ##     ## ##       ##    ##    ##    ##     ## ##     ##
// #### #### ##     ## ##          ##     ##  ##     ## ##       ##          ##    ##     ## ##     ##
// ## ### ## ##     ## ##          ##     ##  ##     ## ######   ##          ##    ##     ## ########
// ##     ## ##     ## ##          ##     ##   ##   ##  ##       ##          ##    ##     ## ##   ##
// ##     ## ##     ## ##          ##     ##    ## ##   ##       ##    ##    ##    ##     ## ##    ##
// ##     ##  #######  ########    ##    ####    ###    ########  ######     ##     #######  ##     ##

MultiVector :: struct(bases: ..B) {
    values: [bases.count]T;
}



BasisHelper :: struct(type: Type) {
    BN :: BasisNumber(type.basis);// #run get_basis(type);
}
// get_basis :: ($type: Type) -> Type {
//     return MultiVector(type.basis);
// }

operator + :: (a: $_A/BasisNumber, b: $_B/BasisNumber) -> BasisHelper(_A).BN #modify {
    ba := get_const_value(_A, "basis", B);
    bb := get_const_value(_B, "basis", B);
    return ba == bb;
} {
    res : BasisNumber(_A.basis) = ---;
    res.value = a.value + b.value;
    return res;
}
operator - :: (a: $_A/BasisNumber, b: $_B/BasisNumber) -> BasisHelper(_A).BN #modify {
    ba := get_const_value(_A, "basis", B);
    bb := get_const_value(_B, "basis", B);
    return ba == bb;
} {
    res : BasisNumber(_A.basis) = ---;
    res.value = a.value - b.value;
    return res;
}

MultiVectorHelper :: struct(a: Type, b: Type) {
    MV :: MultiVector(a.basis, b.basis);
}

operator + :: (a: $_A/BasisNumber, b: $_B/BasisNumber) -> MultiVectorHelper(_A,_B).MV #modify {
    ba := get_const_value(_A, "basis", B);
    bb := get_const_value(_B, "basis", B);
    return ba != bb;
} {
    #if _A.basis < _B.basis {
        res : MultiVector(_A.basis, _B.basis) = ---;
        res.values[0] = a.value;
        res.values[1] = b.value;
    } else {
        res : MultiVector(_B.basis, _A.basis) = ---;
        res.values[0] = b.value;
        res.values[1] = a.value;
    }
    return res;
}
operator - :: (a: $_A/BasisNumber, b: $_B/BasisNumber) -> MultiVectorHelper(_A,_B).MV #modify {
    ba := get_const_value(_A, "basis", B);
    bb := get_const_value(_B, "basis", B);
    return ba != bb;
} {
    #if _A.basis < _B.basis {
        res : MultiVector(_A.basis, _B.basis) = ---;
        res.values[0] = a.value;
        res.values[1] = -b.value;
    } else {
        res : MultiVector(_B.basis, _A.basis) = ---;
        res.values[0] = -b.value;
        res.values[1] = a.value;
    }
    return res;
}

operator + :: (a: $_A/BasisNumber, v: T) -> MultiVectorHelper(BN0, _A).MV #symmetric {
    res : MultiVector(0, _A.basis) = ---;
    res.values[0] = v;
    res.values[1] = a.value;
    return res;
}
operator - :: (a: $_A/BasisNumber, v: T) -> MultiVectorHelper(BN0, _A).MV {
    res : MultiVector(0, _A.basis) = ---;
    res.values[0] = -v;
    res.values[1] = a.value;
    return res;
}
operator - :: (v: T, a: $_A/BasisNumber) -> MultiVectorHelper(BN0, _A).MV {
    res : MultiVector(0, _A.basis) = ---;
    res.values[0] = v;
    res.values[1] = -a.value;
    return res;
}


BN_MV_Helper :: struct(bn: Type, mv: Type) {
    MV :: #run combine_bases_bn_mv(bn, mv);
}
combine_bases_bn_mv :: ($bn: Type, $mv: Type) -> Type {
    DEBUG :: false;

    for mv.bases {
        #if DEBUG println("checking % == %? %", it, bn.basis, it == bn.basis);
        if it == bn.basis then return MultiVector(..mv.bases);
    }

    #if DEBUG println("no same basis found.");

    b :: #run () -> [mv.bases.count+1]B {
        res : [mv.bases.count+1]B = ---;

        mv_indx := 0;
        for res {
            if mv_indx < mv.bases.count && mv.bases[mv_indx] < bn.basis {
                res[it_index] = mv.bases[mv_indx];
            } else {
                res[it_index] = bn.basis;
            }
            mv_indx += 1;
        }


        return res;
    }();
    
    #if DEBUG println("combined bases: % from % and %", b, mv.bases, bn.basis);
    return MultiVector(..b);
}

operator + :: (a: $_A/BasisNumber, b: $_B/MultiVector) -> BN_MV_Helper(_A, _B).MV #symmetric {
    res : BN_MV_Helper(_A, _B).MV = ---;
    #if res.bases.count == _B.bases.count {
        for res.bases {
            res.values[it_index] = b.values[it_index];
            if it == a.basis
                res.values[it_index] += a.value;
        }
    } else {
        for b.values {
            res.values[it_index] = it;
        }
        res.values[res.values.count-1] = a.value;
    }
    return res;
}
operator - :: (a: $_A/BasisNumber, b: $_B/MultiVector) -> BN_MV_Helper(_A, _B).MV {
    res : BN_MV_Helper(_A, _B).MV = ---;
    #if res.bases.count == _B.bases.count {
        for res.bases {
            res.values[it_index] = -b.values[it_index];
            if it == a.basis
                res.values[it_index] += a.value;
        }
    } else {
        for b.values {
            res.values[it_index] = -it;
        }
        res.values[res.values.count-1] = a.value;
    }
    return res;
}
operator - :: (b: $_B/MultiVector, a: $_A/BasisNumber) -> BN_MV_Helper(_A, _B).MV {
    res : BN_MV_Helper(_A, _B).MV = ---;
    #if res.bases.count == _B.bases.count {
        for res.bases {
            res.values[it_index] = b.values[it_index];
            if it == a.basis
                res.values[it_index] -= a.value;
        }
    } else {
        for b.values {
            res.values[it_index] = it;
        }
        res.values[res.values.count-1] = -a.value;
    }
    return res;
}




MV_MV_Helper :: struct(a: Type, b: Type) {
    MV :: #run combine_bases_mv_mv(a, b);
}
combine_bases_mv_mv :: ($atype: Type, $btype: Type) -> Type {
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

        for atype.bases {
            res[it_index] = it;
        }
        indx := atype.bases.count;
        for bb: btype.bases {
            for ab: atype.bases {
                if bb == ab then continue bb;
            }
            res[indx] = bb;
            indx += 1;
        }
        return res;
    }();
    return MultiVector(..b);
}


operator + :: (a: $_A/MultiVector, b: $_B/MultiVector) -> MV_MV_Helper(_A, _B).MV {
    res : MV_MV_Helper(_A, _B).MV;
    // assuming the bases are sorted from small to large value
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
operator - :: (a: $_A/MultiVector, b: $_B/MultiVector) -> MV_MV_Helper(_A, _B).MV {
    res : MV_MV_Helper(_A, _B).MV;
    // assuming the bases are sorted from small to large value
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



operator + :: (a: $_A/MultiVector, b: T) -> BN_MV_Helper(BN0, _A).MV #symmetric {
    tmp : BN0 = ---;
    tmp.value = b;
    return a + tmp;
}
operator - :: (a: $_A/MultiVector, b: T) -> BN_MV_Helper(BN0, _A).MV {
    tmp : BN0 = ---;
    tmp.value = -b;
    return a + tmp;
}
operator - :: (b: T, a: $_A/MultiVector) -> BN_MV_Helper(BN0, _A).MV {
    tmp : BN0 = ---;
    tmp.value = b;
    return tmp - a;
}


//  ######   ########  #######     ########  ########   #######  ########
// ##    ##  ##       ##     ##    ##     ## ##     ## ##     ## ##     ##
// ##        ##       ##     ##    ##     ## ##     ## ##     ## ##     ##
// ##   #### ######   ##     ##    ########  ########  ##     ## ##     ##
// ##    ##  ##       ##     ##    ##        ##   ##   ##     ## ##     ##
// ##    ##  ##       ##     ##    ##        ##    ##  ##     ## ##     ##
//  ######   ########  #######     ##        ##     ##  #######  ########

MV_MV_GEO_PROD_Helper :: struct(a: Type, b: Type) {
    MV :: #run combine_bases_mv_mv_geo_prod(a, b);
}
combine_bases_mv_mv_geo_prod :: ($atype: Type, $btype: Type) -> Type {
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

        for atype.bases {
            res[it_index] = it;
        }
        indx := atype.bases.count;
        for bb: btype.bases {
            for ab: atype.bases {
                if bb == ab then continue bb;
            }
            res[indx] = bb;
            indx += 1;
        }
        return res;
    }();
    return MultiVector(..b);
}

operator * :: (a: $_A/MultiVector, b: $_B/MultiVector) -> MV_MV_GEO_PROD_Helper(_A,_B).MV {
    res : MV_MV_GEO_PROD_Helper(_A,_B).MV = ---;
    return res;
}
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
        if n.basis & (1 << it) {
            print_to_builder(*builder, "e%", it+1);
        }
    }

    return builder_to_string(*builder);
}

// e :: inline ($indices: ..int) -> BasisNumber($bn) #modify {
//     bn = cast(B)0;
//     for indices
//         bn |= cast(B)1 << cast(B)(it-1);
//     return true;
// } {
//     return BasisNumber(bn).{one};
// }

operator * :: (a: $_A/BasisNumber, b: $_B/BasisNumber) -> $C #modify {
    DEBUG :: false;
    ba := get_const_value(_A, "basis", B);
    bb := get_const_value(_B, "basis", B);
    #if DEBUG println("ba = %, bb = %", ba, bb);
    
    overlap := ba & bb;
    if overlap >> (POS+NEG) != 0 {
        C = T;
        return true;
    }
    if ba ^ bb == 0 { // all basis elements square to +-1
        C = T;
        return true;
    }

    C = _A;
    set_const_value(cast(*Type_Info_Struct)C, "basis", ba ^ bb);
    #if DEBUG print("C = %\n", C);
    return true;
} {
    DEBUG :: true;
    #if DEBUG print("% * %\n", fmt(a), fmt(b));
    
    // a basis element squares to zero
    overlap :: _A.basis & _B.basis;
    #if DEBUG print("overlap %\n", formatInt(overlap,2));
    #if overlap >> (POS+NEG) != 0 then return zero;

    #if _A.basis ^ _B.basis == 0 { // all basis elements square to +-1
        parity :: #run one - cast(T)2 * cast(T)(_A.basis > _B.basis);
        negs :: #run (overlap >> POS) & ((1 << NEG)-1);
        negs_squared :: #run 1 - 2 * ((hamming_weight(negs)) % 2);
        return a.value * b.value * parity * negs_squared;
    }

    // not squaring to zero, so calculate actual square
    res : C = ---;
    #if DEBUG println("res : % = %", type_of(res), res);

    // // default ordering from small basis number to large
    // parity :: #run one - cast(T)2 * cast(T)(_A.basis > _B.basis);
    // res.value = a.value * b.value * parity;

    // // no basis elements square up
    // #if overlap == 0 then return res;

    // // basis elements that square to -1
    // negs :: #run (overlap >> POS) & ((1 << NEG)-1);
    // #if DEBUG print("negs %\n", formatInt(negs,2));
    // #if negs {
    //     negs_squared :: #run 1 - 2 * ((hamming_weight(negs)) % 2);
    //     res.value *= negs_squared;
    // }

    // // the rest squares to +1
    return res;
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

// operator - :: (a: BasisNumber) -> BasisNumber {
//     res : BasisNumber = ---;
//     res.basis = a.basis;
//     res.value = -a.value;
//     return res;
// }
operator * :: (a: $_A/BasisNumber, b: T) -> _A #symmetric {
    res : _A = ---;
    // res.basis = a.basis;
    res.value = a.value * b;
    return res;
}

// operator / :: (a: BasisNumber, b: T) -> BasisNumber {
//     res : BasisNumber = ---;
//     res.basis = a.basis;
//     res.value = a.value / b;
//     return res;
// }



// N_ELEMENTS :: #run powi(2, N+1);
// #if STACK_ALLOC {
//     ELEMENT_INDEX_MAP :: #run ()->[N_ELEMENTS]int {
//         res : [N_ELEMENTS]int;
//         return res;
//     }();    
// }


// ##     ## ##     ## ##       ######## #### ##     ## ########  ######  ########  #######  ########
// ###   ### ##     ## ##          ##     ##  ##     ## ##       ##    ##    ##    ##     ## ##     ##
// #### #### ##     ## ##          ##     ##  ##     ## ##       ##          ##    ##     ## ##     ##
// ## ### ## ##     ## ##          ##     ##  ##     ## ######   ##          ##    ##     ## ########
// ##     ## ##     ## ##          ##     ##   ##   ##  ##       ##          ##    ##     ## ##   ##
// ##     ## ##     ## ##          ##     ##    ## ##   ##       ##    ##    ##    ##     ## ##    ##
// ##     ##  #######  ########    ##    ####    ###    ########  ######     ##     #######  ##     ##

// MultiVector :: struct {
//     val: T;
//     #if STACK_ALLOC {
//         terms: [N_ELEMENTS]BasisNumber;
//     } else {
//         terms: [..]BasisNumber;
//     }
// }
// mv :: inline (nums: ..BasisNumber) -> MultiVector {
//     return mv(zero, ..nums);
// }
// mv :: (a: T, nums: ..BasisNumber) -> MultiVector {
//     res : MultiVector;

//     res.val = a;

//     #if STACK_ALLOC {
//         for nums {
//             res.terms[ELEMENT_INDEX_MAP[it.basis]] = it;
//         }
//     } else {
//         for nums {
//             indx := -1;
//             for t, i: res.terms {
//                 if t.basis != it.basis continue;
//                 indx = i;
//                 break;
//             }
            
//             if indx != -1 {
//                 res.terms[indx].value += it.value;
//                 continue;
//             }

//             array_add(*res.terms, it);
//         }
//     }
    
//     return res;
// }
// fmt :: (mv: MultiVector) -> string {
//     builder: String_Builder;

//     print_to_builder(*builder, "%", mv.val);

//     if mv.terms.count == 0 then
//         return builder_to_string(*builder);

//     for mv.terms {
//         if it.value >= 0
//             print_to_builder(*builder, " + %",  fmt(it, nosign=true));
//         else 
//             print_to_builder(*builder, " - %",  fmt(it, nosign=true));
//     }

//     return builder_to_string(*builder);
// }

// operator - :: (a: MultiVector) -> MultiVector {
//     res : MultiVector;
//     res.val = -a.val;
//     for a.terms {
//         array_add(*res.terms, -it);
//     }
//     return res;
// }

// operator + :: (a: BasisNumber, b: T) -> MultiVector #symmetric {
//     res : MultiVector;
//     res.val = b;
//     array_add(*res.terms, a);
//     return res;
// }
// operator - :: (a: BasisNumber, b: T) -> MultiVector {
//     res : MultiVector;
//     res.val = -b;
//     array_add(*res.terms, a);
//     return res;
// }
// operator - :: (b: T, a: BasisNumber) -> MultiVector {
//     res : MultiVector;
//     res.val = b;
//     array_add(*res.terms, -a);
//     return res;
// }

// operator + :: (a: MultiVector, b: MultiVector) -> MultiVector {
//     res : MultiVector = ---;
//     res.val = a.val + b.val;
//     res.terms = NewResizableArray(a.terms.count, BasisNumber,  false);

//     for a.terms {
//         array_add(*res.terms, it);
//     }

//     for bt, bi: b.terms {
//         for at, ai: a.terms {
//             if at.basis == bt.basis {
//                 res.terms[ai].value += bt.value;
//                 continue bt;
//             }
//         }
//         // no matching term found
//         array_add(*res.terms, bt);
//     }

//     return res;
// }

// operator + :: (a: MultiVector, b: BasisNumber) -> MultiVector #symmetric {
//     res : MultiVector;
//     array_copy(*res.terms, a.terms);
//     res.val = a.val;

//     found := false;
//     for a.terms {
//         if it.basis == b.basis {
//             res.terms[it_index].value += b.value;
//             found = true;
//             break;
//         }
//     }
//     if !found then array_add(*res.terms, b);

//     return res;
// }
// operator - :: (a: MultiVector, b: BasisNumber) -> MultiVector {
//     res : MultiVector;
//     array_copy(*res.terms, a.terms);
//     res.val = a.val;

//     found := false;
//     for a.terms {
//         if it.basis == b.basis {
//             res.terms[it_index].value -= b.value;
//             found = true;
//             break;
//         }
//     }
//     if !found then array_add(*res.terms, -b);

//     return res;
// }
// operator - :: (b: BasisNumber, a: MultiVector) -> MultiVector {
//     res : MultiVector;

//     for a.terms {
//         array_add(*res.terms, -it.value);
//     }
//     res.val = -a.val;

//     found := false;
//     for a.terms {
//         if it.basis == b.basis {
//             res.terms[it_index].value += b.value;
//             found = true;
//             break;
//         }
//     }
//     if !found then array_add(*res.terms, b);

//     return res;
// }



// operator * :: (a: MultiVector, b: MultiVector) -> MultiVector {
//     res : MultiVector;
//     res. val = a.val * b.val;

//     for at: a.terms {
//         for bt: b.terms {
//             res += at * bt;
//         }
//     }

//     return res;
// }
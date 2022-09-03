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

Number :: struct {
    basis: B;
    value: T;
}
number :: (val: T, basis: B) -> Number {
    res : Number;
    res.value = val;
    res.basis = basis;
    return res;
}
fmt :: (n: Number, $nosign := false) -> string {
    builder: String_Builder;

    #if nosign {
        v := ifx n.value >= 0 then n.value else -n.value;
        print_to_builder(*builder, "% ", v);
    } else 
        print_to_builder(*builder, "% ", n.value);


    // if n.basis then append(*builder, " e");
    for 0..N-1 {
        if n.basis & (1 << it) {
            print_to_builder(*builder, "e%", it);
        }
    }

    return builder_to_string(*builder);
}

e :: (indices: ..int) -> Number {
    b : B;
    for indices
        b |= cast(B)1 << cast(B)it;
    return number(one, b);
}

e0 :: #run number(one, cast(B)0);

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
                    print_to_builder(*builder, "e% :: #run e(%);\n", it+1, it+1);
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

operator - :: (a: Number) -> Number {
    res : Number = ---;
    res.basis = a.basis;
    res.value = -a.value;
    return res;
}
operator * :: (a: Number, b: T) -> Number #symmetric {
    res : Number = ---;
    res.basis = a.basis;
    res.value = a.value * b;
    return res;
}

operator / :: (a: Number, b: T) -> Number {
    res : Number = ---;
    res.basis = a.basis;
    res.value = a.value / b;
    return res;
}

operator * :: (a: Number, b: Number) -> Number {
    DEBUG :: true;
    #if DEBUG print("% * %\n", fmt(a), fmt(b));

    res : Number = ---;
    res.basis = a.basis ^ b.basis;
    #if DEBUG println("basis %", formatInt(res.basis, 2));
    res.value = a.value * b.value * (1.0 - 2.0 * cast(float32)(a.basis > b.basis));

    overlap := a.basis & b.basis;

    // no basis elements square up
    if overlap == 0 then return res;

    // a basis element squares to zero
    #if DEBUG print("overlap %\n", formatInt(overlap,2));
    if overlap >> (POS+NEG) != 0 {
        res.value = zero;
        res.basis = 0;
        return res;
    }

    // basis elements that square to -1
    negs := (overlap >> POS) & ((1 << NEG)-1);
    #if DEBUG print("negs %\n", formatInt(negs,2));
    if negs {
        res.value *= 1 - 2 * ((hamming_weight(negs)) % 2);
    }

    // the rest squares to +1
    return res;
}



N_ELEMENTS :: #run powi(2, N+1);
#if STACK_ALLOC {
    ELEMENT_INDEX_MAP :: #run ()->[N_ELEMENTS]int {
        res : [N_ELEMENTS]int;
        return res;
    }();    
}

MultiVector :: struct {
    val: T;
    #if STACK_ALLOC {
        terms: [N_ELEMENTS]Number;
    } else {
        terms: [..]Number;
    }
}
mv :: inline (nums: ..Number) -> MultiVector {
    return mv(zero, ..nums);
}
mv :: (a: T, nums: ..Number) -> MultiVector {
    res : MultiVector;

    res.val = a;

    #if STACK_ALLOC {
        for nums {
            res.terms[ELEMENT_INDEX_MAP[it.basis]] = it;
        }
    } else {
        for nums {
            indx := -1;
            for t, i: res.terms {
                if t.basis != it.basis continue;
                indx = i;
                break;
            }
            
            if indx != -1 {
                res.terms[indx].value += it.value;
                continue;
            }

            array_add(*res.terms, it);
        }
    }
    
    return res;
}
fmt :: (mv: MultiVector) -> string {
    builder: String_Builder;

    print_to_builder(*builder, "%", mv.val);

    if mv.terms.count == 0 then
        return builder_to_string(*builder);

    for mv.terms {
        if it.value >= 0
            print_to_builder(*builder, " + %",  fmt(it, nosign=true));
        else 
            print_to_builder(*builder, " - %",  fmt(it, nosign=true));
    }

    return builder_to_string(*builder);
}

operator - :: (a: MultiVector) -> MultiVector {
    res : MultiVector;
    res.val = -a.val;
    for a.terms {
        array_add(*res.terms, -it);
    }
    return res;
}

operator + :: (a: Number, b: T) -> MultiVector #symmetric {
    res : MultiVector;
    res.val = b;
    array_add(*res.terms, a);
    return res;
}
operator - :: (a: Number, b: T) -> MultiVector {
    res : MultiVector;
    res.val = -b;
    array_add(*res.terms, a);
    return res;
}
operator - :: (b: T, a: Number) -> MultiVector {
    res : MultiVector;
    res.val = b;
    array_add(*res.terms, -a);
    return res;
}

operator + :: (a: MultiVector, b: MultiVector) -> MultiVector {
    res : MultiVector = ---;
    res.val = a.val + b.val;
    res.terms = NewResizableArray(a.terms.count, Number,  false);

    for a.terms {
        array_add(*res.terms, it);
    }

    for bt, bi: b.terms {
        for at, ai: a.terms {
            if at.basis == bt.basis {
                res.terms[ai].value += bt.value;
                continue bt;
            }
        }
        // no matching term found
        array_add(*res.terms, bt);
    }

    return res;
}

operator + :: (a: MultiVector, b: Number) -> MultiVector #symmetric {
    res : MultiVector;
    array_copy(*res.terms, a.terms);
    res.val = a.val;

    found := false;
    for a.terms {
        if it.basis == b.basis {
            res.terms[it_index].value += b.value;
            found = true;
            break;
        }
    }
    if !found then array_add(*res.terms, b);

    return res;
}
operator - :: (a: MultiVector, b: Number) -> MultiVector {
    res : MultiVector;
    array_copy(*res.terms, a.terms);
    res.val = a.val;

    found := false;
    for a.terms {
        if it.basis == b.basis {
            res.terms[it_index].value -= b.value;
            found = true;
            break;
        }
    }
    if !found then array_add(*res.terms, -b);

    return res;
}
operator - :: (b: Number, a: MultiVector) -> MultiVector {
    res : MultiVector;

    for a.terms {
        array_add(*res.terms, -it.value);
    }
    res.val = -a.val;

    found := false;
    for a.terms {
        if it.basis == b.basis {
            res.terms[it_index].value += b.value;
            found = true;
            break;
        }
    }
    if !found then array_add(*res.terms, b);

    return res;
}



operator * :: (a: MultiVector, b: MultiVector) -> MultiVector {
    res : MultiVector;
    res. val = a.val * b.val;

    for at: a.terms {
        for bt: b.terms {
            res += at * bt;
        }
    }

    return res;
}
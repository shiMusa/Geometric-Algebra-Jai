#import "Basic";
#import "Math";
#load "utils.jai";
#module_parameters(POS: int, NEG: int, ZER: int, T: Type = float32);

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
    #run assert(N < 64, "Only basis sizes of up to 64 supported, but chosen %.", N);
}

#run print("Generating algebra with basis (%,%,%). Basis elements represented as %\n", POS, NEG, ZER, B);

GA_Number :: struct {
    basis: B;
    value: T;
}
ga_number :: (val: T, basis: B) -> GA_Number {
    res : GA_Number;
    res.basis = basis;
    res.value = val;
    return res;
}
ga_el :: (indices: ..int) -> GA_Number {
    b : B;
    for indices
        b |= cast(B)1 << cast(B)it;
    return ga_number(one, b);
}

#insert -> string {
    builder : String_Builder;

    for grade: 1..N {
        print("grade %:\n", grade);

        NN := powi(N, grade);
        print("    NN = %\n", NN);
        
        // grade_dims := NewArray(grade, int); defer array_free(grade_dims);
        // for g: 0..grade-1 {
        //     grade_dims[g] = powi(N,g+1);
        // }
        // print("    grade_dims = %\n", grade_dims);

        indices := NewArray(grade, int); defer array_free(indices);

        // iterate over all possible index combinations
        for i: 0..NN-1 {
            print("    %: ", i);

            // increase indices
            indices[0] = i;
            for g: 0..grade-1-1 {
                if indices[g] % N == 0 then
                    indices[g+1] += 1;
            }
            for *indices { <<it %= N; }
            print("%\n", indices);

            // skip if the same element appears twice
            for i0: 0..grade-1 {
                for j0: i0+1..grade-1 {
                    if indices[i0] == indices[j0] {
                        print("        same elements % & %, continue\n", i0, j0);
                        continue i;
                    }
                }
            }

            append(*builder, "e");
            for indices {
                print_to_builder(*builder, "%", it+1);
            }
            append(*builder, " :: #run ga_el(");
            for indices {
                print_to_builder(*builder, "%", it+1);
                if it_index != indices.count-1 then
                    append(*builder, ",");
            }
            append(*builder, ");\n");
        }
    }

    // for 0..N-1 {
    //     print_to_builder(*builder, "e%1 :: #run ga_el(%1);\n", it);
    // }

    // for i: 0..N-1 {
    //     for j: 0..N-1 {
    //         if i == j then continue;
    //         print_to_builder(*builder, "e%1%2 :: #run ga_el(%1,%2);\n", i, j);
    //     }
    // }

    return builder_to_string(*builder);
}